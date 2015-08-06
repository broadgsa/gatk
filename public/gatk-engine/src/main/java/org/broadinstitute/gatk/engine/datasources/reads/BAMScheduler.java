/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.engine.datasources.reads;

import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.GATKBAMFileSpan;
import htsjdk.samtools.GATKChunk;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.GenomeLocSortedSet;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.interval.IntervalMergingRule;
import org.broadinstitute.gatk.utils.sam.ReadUtils;
import org.broadinstitute.gatk.utils.sam.SAMReaderID;

import java.util.*;

/**
 * Assign intervals to the most appropriate blocks, keeping as little as possible in memory at once.
 */
public class BAMScheduler implements Iterator<FilePointer> {
    private final SAMDataSource dataSource;

    private final Map<SAMReaderID,GATKBAMIndex> indexFiles = new HashMap<SAMReaderID,GATKBAMIndex>();

    private FilePointer nextFilePointer = null;

    private GenomeLocSortedSet loci;
    private PeekableIterator<GenomeLoc> locusIterator;
    private GenomeLoc currentLocus;
    private IntervalMergingRule intervalMergingRule;

    /*
     * Creates BAMScheduler using contigs from the given BAM data source.
     *
     * @param dataSource    BAM source
     * @return non-null BAM scheduler
     */
    public static BAMScheduler createOverMappedReads(final SAMDataSource dataSource) {
        final BAMScheduler scheduler = new BAMScheduler(dataSource, IntervalMergingRule.ALL);
        final GenomeLocSortedSet intervals = GenomeLocSortedSet.createSetFromSequenceDictionary(dataSource.getHeader().getSequenceDictionary());
        scheduler.populateFilteredIntervalList(intervals);
        return scheduler;
    }

    public static BAMScheduler createOverAllReads(final SAMDataSource dataSource, final GenomeLocParser parser) {
        BAMScheduler scheduler = new BAMScheduler(dataSource, IntervalMergingRule.ALL);
        scheduler.populateUnfilteredIntervalList(parser);
        return scheduler;
    }

    public static BAMScheduler createOverIntervals(final SAMDataSource dataSource, final IntervalMergingRule mergeRule, final GenomeLocSortedSet loci) {
        BAMScheduler scheduler = new BAMScheduler(dataSource, mergeRule);
        scheduler.populateFilteredIntervalList(loci);
        return scheduler;
    }


    private BAMScheduler(final SAMDataSource dataSource, final IntervalMergingRule mergeRule) {
        this.dataSource = dataSource;
        this.intervalMergingRule = mergeRule;
        for(SAMReaderID reader: dataSource.getReaderIDs()) {
            GATKBAMIndex index = dataSource.getIndex(reader);
            if(index != null)
                indexFiles.put(reader,dataSource.getIndex(reader));
        }
    }

    /**
     * The consumer has asked for a bounded set of locations.  Prepare an iterator over those locations.
     * @param loci The list of locations to search and iterate over.
     */
    private void populateFilteredIntervalList(final GenomeLocSortedSet loci) {
        this.loci = loci;
        if(!indexFiles.isEmpty()) {
            // If index data is available, start up the iterator.
            locusIterator = new PeekableIterator<GenomeLoc>(loci.iterator());
            if(locusIterator.hasNext())
                currentLocus = locusIterator.next();
            advance();
        }
        else {
            // Otherwise, seed the iterator with a single file pointer over the entire region.
            nextFilePointer = generatePointerOverEntireFileset();
            for(GenomeLoc locus: loci)
                nextFilePointer.addLocation(locus);
            locusIterator = new PeekableIterator<GenomeLoc>(Collections.<GenomeLoc>emptyList().iterator());
        }
    }

    /**
     * The consumer has provided null, meaning to iterate over all available data.  Create a file pointer stretching
     * from just before the start of the region to the end of the region.
     */
    private void populateUnfilteredIntervalList(final GenomeLocParser parser) {
        this.loci = new GenomeLocSortedSet(parser);
        locusIterator = new PeekableIterator<GenomeLoc>(Collections.<GenomeLoc>emptyList().iterator());
        nextFilePointer = generatePointerOverEntireFileset();
    }

    /**
     * Generate a span that runs from the end of the BAM header to the end of the fle.
     * @return A file pointer over the specified region.
     */
    private FilePointer generatePointerOverEntireFileset() {
        FilePointer filePointer = new FilePointer(intervalMergingRule);

        // This is a "monolithic" FilePointer representing all regions in all files we will ever visit, and is
        // the only FilePointer we will create. This allows us to have this FilePointer represent regions from
        // multiple contigs
        filePointer.setIsMonolithic(true);

        Map<SAMReaderID,GATKBAMFileSpan> currentPosition;

        currentPosition = dataSource.getInitialReaderPositions();

        for(SAMReaderID reader: dataSource.getReaderIDs())
            filePointer.addFileSpans(reader,createSpanToEndOfFile(currentPosition.get(reader).getGATKChunks().get(0).getChunkStart()));
        return filePointer;
    }

    public boolean hasNext() {
        return nextFilePointer != null;
    }

    public FilePointer next() {
        if(!hasNext())
            throw new NoSuchElementException("No next element available in interval sharder");
        FilePointer currentFilePointer = nextFilePointer;
        nextFilePointer = null;
        advance();

        return currentFilePointer;
    }

    public void remove() {
        throw new UnsupportedOperationException("Unable to remove FilePointers from an IntervalSharder");
    }

    private void advance() {
        if(loci.isEmpty())
            return;

        while(nextFilePointer == null && currentLocus != null) {
            // special case handling of the unmapped shard.
            if(currentLocus == GenomeLoc.UNMAPPED) {
                nextFilePointer = new FilePointer(intervalMergingRule, GenomeLoc.UNMAPPED);
                for(SAMReaderID id: dataSource.getReaderIDs())
                    nextFilePointer.addFileSpans(id,createSpanToEndOfFile(indexFiles.get(id).getStartOfLastLinearBin()));
                currentLocus = null;
                continue;
            }

            nextFilePointer = new FilePointer(intervalMergingRule);

            int coveredRegionStart = 1;
            int coveredRegionStop = Integer.MAX_VALUE;
            GenomeLoc coveredRegion = null;

            BAMScheduleEntry scheduleEntry = getNextOverlappingBAMScheduleEntry(currentLocus);

            // No overlapping data at all.
            if(scheduleEntry != null) {
                coveredRegionStart = Math.max(coveredRegionStart,scheduleEntry.start);
                coveredRegionStop = Math.min(coveredRegionStop,scheduleEntry.stop);
                coveredRegion = loci.getGenomeLocParser().createGenomeLoc(currentLocus.getContig(),coveredRegionStart,coveredRegionStop);

                nextFilePointer.addFileSpans(scheduleEntry.fileSpans);
            }
            else {
                // Always create a file span, whether there was covered data or not.  If there was no covered data, then the binTree is empty.
                for(SAMReaderID reader: indexFiles.keySet())
                    nextFilePointer.addFileSpans(reader,new GATKBAMFileSpan());
            }

            // Early exit if no bins were found.
            if(coveredRegion == null) {
                // for debugging only: maximum split is 16384.                
                nextFilePointer.addLocation(currentLocus);
                currentLocus = locusIterator.hasNext() ? locusIterator.next() : null;
                continue;
            }

            // Early exit if only part of the first interval was found.
            if(currentLocus.startsBefore(coveredRegion)) {
                int splitPoint = Math.min(coveredRegion.getStart()-currentLocus.getStart(),16384)+currentLocus.getStart();
                GenomeLoc[] splitContigs = currentLocus.split(splitPoint);
                nextFilePointer.addLocation(splitContigs[0]);
                currentLocus = splitContigs[1];
                continue;
            }

            // Define the initial range of the file pointer, aka the region where the locus currently being processed intersects the BAM list.
            GenomeLoc initialLocation = currentLocus.intersect(coveredRegion);
            nextFilePointer.addLocation(initialLocation);

            // See whether the BAM regions discovered overlap the next set of intervals in the interval list.  If so, include every overlapping interval.
            if(!nextFilePointer.locations.isEmpty()) {
                while(locusIterator.hasNext() && locusIterator.peek().overlapsP(coveredRegion)) {
                    currentLocus = locusIterator.next();
                    nextFilePointer.addLocation(currentLocus.intersect(coveredRegion));
                }

                // Chop off the uncovered portion of the locus.  Since we know that the covered region overlaps the current locus,
                  // we can simplify the interval creation process to the end of the covered region to the stop of the given interval.
                if(coveredRegionStop < currentLocus.getStop())
                    currentLocus = loci.getGenomeLocParser().createGenomeLoc(currentLocus.getContig(),coveredRegionStop+1,currentLocus.getStop());
                else if(locusIterator.hasNext())
                    currentLocus = locusIterator.next();
                else
                    currentLocus = null;
            }

        }
    }

    
    /**
     * The last reference sequence processed by this iterator.
     */
    private Integer lastReferenceSequenceLoaded = null;

    /**
     * The stateful iterator used to progress through the genoem.
     */
    private PeekableIterator<BAMScheduleEntry> bamScheduleIterator = null;

    /**
     * Clean up underlying BAMSchedule file handles.
     */
    public void close() {
        if(bamScheduleIterator != null)
            bamScheduleIterator.close();
    }

    /**
     * Get the next overlapping tree of bins associated with the given BAM file.
     * @param currentLocus The actual locus for which to check overlap.
     * @return The next schedule entry overlapping with the given list of loci.
     */
    private BAMScheduleEntry getNextOverlappingBAMScheduleEntry(final GenomeLoc currentLocus) {
        // Make sure that we consult the BAM header to ensure that we're using the correct contig index for this contig name.
        // This will ensure that if the two sets of contigs don't quite match (b36 male vs female ref, hg19 Epstein-Barr), then
        // we'll be using the correct contig index for the BAMs.
        // TODO: Warning: assumes all BAMs use the same sequence dictionary!  Get around this with contig aliasing.
        SAMSequenceRecord currentContigSequenceRecord = dataSource.getHeader().getSequence(currentLocus.getContig());
        if ( currentContigSequenceRecord == null ) {
            throw new UserException(String.format("Contig %s not present in sequence dictionary for merged BAM header: %s",
                                                  currentLocus.getContig(),
                                                  ReadUtils.prettyPrintSequenceRecords(dataSource.getHeader().getSequenceDictionary())));
        }

        final int currentContigIndex = currentContigSequenceRecord.getSequenceIndex();

        // Stale reference sequence or first invocation.  (Re)create the binTreeIterator.
        if(lastReferenceSequenceLoaded == null || lastReferenceSequenceLoaded != currentContigIndex) {
            if(bamScheduleIterator != null)
                bamScheduleIterator.close();
            lastReferenceSequenceLoaded = currentContigIndex;

            // Naive algorithm: find all elements in current contig for proper schedule creation.
            List<GenomeLoc> lociInContig = new LinkedList<GenomeLoc>();
            for(GenomeLoc locus: loci) {
                if (!GenomeLoc.isUnmapped(locus) && dataSource.getHeader().getSequence(locus.getContig()) == null)
                    throw new ReviewedGATKException("BAM file(s) do not have the contig: " + locus.getContig() + ". You are probably using a different reference than the one this file was aligned with");

                if (!GenomeLoc.isUnmapped(locus) && dataSource.getHeader().getSequence(locus.getContig()).getSequenceIndex() == lastReferenceSequenceLoaded)
                    lociInContig.add(locus);
            }

            bamScheduleIterator = new PeekableIterator<BAMScheduleEntry>(new BAMSchedule(dataSource,lociInContig));
        }

        if(!bamScheduleIterator.hasNext())
            return null;

        // Peek the iterator along until finding the first binTree at or following the current locus.
        BAMScheduleEntry bamScheduleEntry = bamScheduleIterator.peek();
        while(bamScheduleEntry != null && bamScheduleEntry.isBefore(currentLocus)) {
            bamScheduleIterator.next();
            bamScheduleEntry = bamScheduleIterator.hasNext() ? bamScheduleIterator.peek() : null;
        }                                   

        return (bamScheduleEntry != null && bamScheduleEntry.overlaps(currentLocus)) ? bamScheduleEntry : null;
    }

    /**
     * Create a span from the given start point to the end of the file.
     * @param startOfRegion Start of the region, in encoded coordinates (block start << 16 & block offset).
     * @return A file span from the given point to the end of the file.
     */
    private GATKBAMFileSpan createSpanToEndOfFile(final long startOfRegion) {
      return new GATKBAMFileSpan(new GATKChunk(startOfRegion,Long.MAX_VALUE));
    }

}

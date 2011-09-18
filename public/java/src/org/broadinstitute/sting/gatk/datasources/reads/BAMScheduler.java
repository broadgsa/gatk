/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.datasources.reads;

import net.sf.picard.util.PeekableIterator;
import net.sf.samtools.GATKBAMFileSpan;
import net.sf.samtools.GATKChunk;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.*;

/**
 * Assign intervals to the most appropriate blocks, keeping as little as possible in memory at once.
 */
public class BAMScheduler implements Iterator<FilePointer> {
    private final SAMDataSource dataSource;

    private final Map<SAMReaderID,GATKBAMIndex> indexFiles = new HashMap<SAMReaderID,GATKBAMIndex>();

    private FilePointer nextFilePointer = null;

    private final GenomeLocSortedSet loci;

    private final PeekableIterator<GenomeLoc> locusIterator;

    private GenomeLoc currentLocus;    

    public BAMScheduler(final SAMDataSource dataSource, final GenomeLocSortedSet loci) {
        this.dataSource = dataSource;
        for(SAMReaderID reader: dataSource.getReaderIDs())
            indexFiles.put(reader,(GATKBAMIndex)dataSource.getIndex(reader));
        this.loci = loci;
        locusIterator = new PeekableIterator<GenomeLoc>(loci.iterator());
        if(locusIterator.hasNext())
            currentLocus = locusIterator.next();
        advance();
    }

    public boolean hasNext() {
        return nextFilePointer != null;
    }

    public FilePointer next() {
        if(!hasNext())
            throw new NoSuchElementException("No next element available in interval sharder");
        FilePointer currentFilePointer = nextFilePointer;
        advance();
        return currentFilePointer;
    }

    public void remove() {
        throw new UnsupportedOperationException("Unable to remove FilePointers from an IntervalSharder");
    }

    private void advance() {
        if(loci.isEmpty())
            return;

        nextFilePointer = null;
        while(nextFilePointer == null && currentLocus != null) {
            // special case handling of the unmapped shard.
            if(currentLocus == GenomeLoc.UNMAPPED) {
                nextFilePointer = new FilePointer(GenomeLoc.UNMAPPED);
                for(SAMReaderID id: dataSource.getReaderIDs())
                    nextFilePointer.addFileSpans(id,new GATKBAMFileSpan(new GATKChunk(indexFiles.get(id).getStartOfLastLinearBin(),Long.MAX_VALUE)));
                currentLocus = null;
                continue;
            }

            nextFilePointer = new FilePointer();

            int coveredRegionStart = 1;
            int coveredRegionStop = Integer.MAX_VALUE;
            GenomeLoc coveredRegion = null;

            BAMScheduleEntry scheduleEntry = getNextOverlappingBAMScheduleEntry(indexFiles,currentLocus);

            // No overlapping data at all.
            if(scheduleEntry != null) {
                coveredRegionStart = Math.max(coveredRegionStart,scheduleEntry.start);
                coveredRegionStop = Math.min(coveredRegionStop,scheduleEntry.stop);
                coveredRegion = loci.getGenomeLocParser().createGenomeLoc(currentLocus.getContig(),coveredRegionStart,coveredRegionStop);

                nextFilePointer.addFileSpans(scheduleEntry.fileSpans);
            }
            else {
                // Always create a file span, whether there was covered data or not.  If there was no covered data, then the binTree is empty.
                //System.out.printf("Shard: index file = %s; reference sequence = %d; ",index.getIndexFile(),currentLocus.getContigIndex());
                for(SAMReaderID reader: indexFiles.keySet())
                    nextFilePointer.addFileSpans(reader,new GATKBAMFileSpan());
            }

            // Early exit if no bins were found.
            if(coveredRegion == null) {
                // for debugging only: maximum split is 16384.                
                if(currentLocus.size() > 16384) {
                    GenomeLoc[] splitContigs = currentLocus.split(currentLocus.getStart()+16384);
                    nextFilePointer.addLocation(splitContigs[0]);
                    currentLocus = splitContigs[1];
                }
                else {
                    nextFilePointer.addLocation(currentLocus);
                    currentLocus = locusIterator.hasNext() ? locusIterator.next() : null;
                }
                continue;
            }

            // Early exit if only part of the first interval was found.
            if(currentLocus.startsBefore(coveredRegion)) {
                // for debugging only: maximum split is 16384.
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
     * Get the next overlapping tree of bins associated with the given BAM file.
     * @param indices BAM indices.
     * @param currentLocus The actual locus for which to check overlap.
     * @return The next schedule entry overlapping with the given list of loci.
     */
    private BAMScheduleEntry getNextOverlappingBAMScheduleEntry(final Map<SAMReaderID,GATKBAMIndex> indices, final GenomeLoc currentLocus) {
        // Stale reference sequence or first invocation.  (Re)create the binTreeIterator.
        if(lastReferenceSequenceLoaded == null || lastReferenceSequenceLoaded != currentLocus.getContigIndex()) {
            if(bamScheduleIterator != null)
                bamScheduleIterator.close();
            lastReferenceSequenceLoaded = currentLocus.getContigIndex();

            // Naive algorithm: find all elements in current contig for proper schedule creation.
            List<GenomeLoc> lociInContig = new LinkedList<GenomeLoc>();
            for(GenomeLoc locus: loci) {
                if(locus.getContigIndex() == lastReferenceSequenceLoaded)
                    lociInContig.add(locus);
            }

            bamScheduleIterator = new PeekableIterator<BAMScheduleEntry>(new BAMSchedule(indices,lociInContig));
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

}

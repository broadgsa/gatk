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
import net.sf.samtools.AbstractBAMFileIndex;
import net.sf.samtools.Bin;
import net.sf.samtools.BrowseableBAMIndex;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.*;

/**
 * Shard intervals based on position within the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalSharder {
    private static Logger logger = Logger.getLogger(IntervalSharder.class);

    public static Iterator<FilePointer> shardIntervals(final SAMDataSource dataSource, final GenomeLocSortedSet loci) {
        return new IntervalSharder.FilePointerIterator(dataSource,loci);
    }

    /**
     * A lazy-loading iterator over file pointers.
     */
    private static class FilePointerIterator implements Iterator<FilePointer> {
        final SAMDataSource dataSource;
        final GenomeLocSortedSet loci;
        final PeekableIterator<GenomeLoc> locusIterator;
        final Queue<FilePointer> cachedFilePointers = new LinkedList<FilePointer>();

        public FilePointerIterator(final SAMDataSource dataSource, final GenomeLocSortedSet loci) {
            this.dataSource = dataSource;
            this.loci = loci;
            locusIterator = new PeekableIterator<GenomeLoc>(loci.iterator());
            advance();
        }

        public boolean hasNext() {
            return !cachedFilePointers.isEmpty();
        }

        public FilePointer next() {
            if(!hasNext())
                throw new NoSuchElementException("FilePointerIterator iteration is complete");
            FilePointer filePointer = cachedFilePointers.remove();
            if(cachedFilePointers.isEmpty())
                advance();
            return filePointer;
        }

        public void remove() {
            throw new UnsupportedOperationException("Cannot remove from a FilePointerIterator");
        }

        private void advance() {
            GenomeLocSortedSet nextBatch = new GenomeLocSortedSet(loci.getGenomeLocParser());
            String contig = null;

            // If the next section of the BAM to be processed is unmapped, handle this region separately.
            while(locusIterator.hasNext() && nextBatch.isEmpty()) {
                contig = null;
                while(locusIterator.hasNext() && (contig == null || (!GenomeLoc.isUnmapped(locusIterator.peek()) && locusIterator.peek().getContig().equals(contig)))) {
                    GenomeLoc nextLocus = locusIterator.next();
                    contig = nextLocus.getContig();
                    nextBatch.add(nextLocus);
                }
            }

            if(nextBatch.size() > 0) {
                cachedFilePointers.addAll(shardIntervalsOnContig(dataSource,contig,nextBatch));
            }
        }
    }

    /**
     * Merge / split intervals based on an awareness of the structure of the BAM file.
     * @param dataSource
     * @param contig Contig against which to align the intervals.  If null, create a file pointer across unmapped reads.
     * @param loci
     * @return
     */
    private static List<FilePointer> shardIntervalsOnContig(final SAMDataSource dataSource, final String contig, final GenomeLocSortedSet loci) {
        // If the contig is null, eliminate the chopping process and build out a file pointer consisting of the unmapped region of all BAMs.
        if(contig == null) {
            FilePointer filePointer = new FilePointer(GenomeLoc.UNMAPPED);
            for(SAMReaderID id: dataSource.getReaderIDs())
                filePointer.addFileSpans(id,null);
            return Collections.singletonList(filePointer);
        }

        // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.
        List<FilePointer> filePointers = new ArrayList<FilePointer>();
        FilePointer lastFilePointer = null;
        BAMOverlap lastBAMOverlap = null;

        Map<SAMReaderID,BrowseableBAMIndex> readerToIndexMap = new HashMap<SAMReaderID,BrowseableBAMIndex>();
        IntervalSharder.BinMergingIterator binMerger = new IntervalSharder.BinMergingIterator();
        for(SAMReaderID id: dataSource.getReaderIDs()) {
            final SAMSequenceRecord referenceSequence = dataSource.getHeader(id).getSequence(contig);
            // If this contig can't be found in the reference, skip over it.
            if(referenceSequence == null && contig != null)
                continue;
            final BrowseableBAMIndex index = (BrowseableBAMIndex)dataSource.getIndex(id);
            binMerger.addReader(id,
                                index,
                                referenceSequence.getSequenceIndex(),
                                index.getBinsOverlapping(referenceSequence.getSequenceIndex(),1,referenceSequence.getSequenceLength()).iterator());
            // Cache the reader for later data lookup.
            readerToIndexMap.put(id,index);
        }

        PeekableIterator<BAMOverlap> binIterator = new PeekableIterator<BAMOverlap>(binMerger);

        for(GenomeLoc location: loci) {
            if(!location.getContig().equals(contig))
                throw new ReviewedStingException("Location outside bounds of contig");

            if(!binIterator.hasNext())
                break;

            int locationStart = location.getStart();
            final int locationStop = location.getStop();

            // Advance to first bin.
            while(binIterator.peek().stop < locationStart)
                binIterator.next();

            // Add all relevant bins to a list.  If the given bin extends beyond the end of the current interval, make
            // sure the extending bin is not pruned from the list.
            List<BAMOverlap> bamOverlaps = new ArrayList<BAMOverlap>();
            while(binIterator.hasNext() && binIterator.peek().stop <= locationStop)
                bamOverlaps.add(binIterator.next());
            if(binIterator.hasNext() && binIterator.peek().start <= locationStop)
                bamOverlaps.add(binIterator.peek());

            // Bins found; try to match bins with locations.
            Iterator<BAMOverlap> bamOverlapIterator = bamOverlaps.iterator();

            while(locationStop >= locationStart) {
                int binStart = lastFilePointer!=null ? lastFilePointer.overlap.start : 0;
                int binStop =  lastFilePointer!=null ? lastFilePointer.overlap.stop : 0;

                while(binStop < locationStart && bamOverlapIterator.hasNext()) {
                    if(lastFilePointer != null && lastFilePointer.locations.size() > 0)
                        filePointers.add(lastFilePointer);

                    lastBAMOverlap = bamOverlapIterator.next();
                    lastFilePointer = new FilePointer(lastBAMOverlap);
                    binStart = lastFilePointer.overlap.start;
                    binStop = lastFilePointer.overlap.stop;
                }

                if(locationStart < binStart) {
                    // The region starts before the first bin in the sequence.  Add the region occurring before the sequence.
                    if(lastFilePointer != null && lastFilePointer.locations.size() > 0) {
                        filePointers.add(lastFilePointer);
                        lastFilePointer = null;
                        lastBAMOverlap = null;
                    }

                    final int regionStop = Math.min(locationStop,binStart-1);

                    GenomeLoc subset = loci.getGenomeLocParser().createGenomeLoc(location.getContig(),locationStart,regionStop);
                    lastFilePointer = new FilePointer(subset);

                    locationStart = regionStop + 1;
                }
                else if(locationStart > binStop) {
                    // The region starts after the last bin in the sequence.  Add the region occurring after the sequence.
                    if(lastFilePointer != null && lastFilePointer.locations.size() > 0) {
                        filePointers.add(lastFilePointer);
                        lastFilePointer = null;
                        lastBAMOverlap = null;
                    }

                    GenomeLoc subset = loci.getGenomeLocParser().createGenomeLoc(location.getContig(),locationStart,locationStop);
                    filePointers.add(new FilePointer(subset));

                    locationStart = locationStop + 1;
                }
                else {
                    if(lastFilePointer == null)
                        throw new ReviewedStingException("Illegal state: initializer failed to create cached file pointer.");

                    // The start of the region overlaps the bin.  Add the overlapping subset.
                    final int regionStop = Math.min(locationStop,binStop);
                    lastFilePointer.addLocation(loci.getGenomeLocParser().createGenomeLoc(location.getContig(),locationStart,regionStop));
                    locationStart = regionStop + 1;
                }
            }
        }

        if(lastFilePointer != null && lastFilePointer.locations.size() > 0)
            filePointers.add(lastFilePointer);

        // Lookup the locations for every file pointer in the index.
        for(SAMReaderID id: readerToIndexMap.keySet()) {
            BrowseableBAMIndex index = readerToIndexMap.get(id);
            for(FilePointer filePointer: filePointers)
                filePointer.addFileSpans(id,index.getSpanOverlapping(filePointer.overlap.getBin(id)));
        }

        return filePointers;
    }

    private static class BinMergingIterator implements Iterator<BAMOverlap> {
        private PriorityQueue<BinQueueState> binQueue = new PriorityQueue<BinQueueState>();
        private Queue<BAMOverlap> pendingOverlaps = new LinkedList<BAMOverlap>();

        public void addReader(final SAMReaderID id, final BrowseableBAMIndex index, final int referenceSequence, Iterator<Bin> bins) {
            binQueue.add(new BinQueueState(id,index,referenceSequence,new IntervalSharder.LowestLevelBinFilteringIterator(index,bins)));
        }

        public boolean hasNext() {
            return pendingOverlaps.size() > 0 || !binQueue.isEmpty();
        }

        public BAMOverlap next() {
            if(!hasNext())
                throw new NoSuchElementException("No elements left in merging iterator");
            if(pendingOverlaps.isEmpty())
                advance();
            return pendingOverlaps.remove();
        }

        public void advance() {
            List<ReaderBin> bins = new ArrayList<ReaderBin>();
            int boundsStart, boundsStop;

            // Prime the pump
            if(binQueue.isEmpty())
                return;
            bins.add(getNextBin());
            boundsStart = bins.get(0).getStart();
            boundsStop  = bins.get(0).getStop();

            // Accumulate all the bins that overlap the current bin, in sorted order.
            while(!binQueue.isEmpty() && peekNextBin().getStart() <= boundsStop) {
                ReaderBin bin = getNextBin();
                bins.add(bin);
                boundsStart = Math.min(boundsStart,bin.getStart());
                boundsStop = Math.max(boundsStop,bin.getStop());
            }

            List<Pair<Integer,Integer>> range = new ArrayList<Pair<Integer,Integer>>();
            int start = bins.get(0).getStart();
            int stop = bins.get(0).getStop();
            while(start <= boundsStop) {
                // Find the next stopping point.
                for(ReaderBin bin: bins) {
                    stop = Math.min(stop,bin.getStop());
                    if(start < bin.getStart())
                        stop = Math.min(stop,bin.getStart()-1);
                }

                range.add(new Pair<Integer,Integer>(start,stop));
                // If the last entry added included the last element, stop.
                if(stop >= boundsStop)
                    break;

                // Find the next start.
                start = stop + 1;
                for(ReaderBin bin: bins) {
                    if(start >= bin.getStart() && start <= bin.getStop())
                        break;
                    else if(start < bin.getStart()) {
                        start = bin.getStart();
                        break;
                    }
                }
            }

            // Add the next series of BAM overlaps to the window.
            for(Pair<Integer,Integer> window: range) {
                BAMOverlap bamOverlap = new BAMOverlap(window.first,window.second);
                for(ReaderBin bin: bins)
                    bamOverlap.addBin(bin.id,bin.bin);
                pendingOverlaps.add(bamOverlap);
            }
        }

        public void remove() { throw new UnsupportedOperationException("Cannot remove from a merging iterator."); }

        private ReaderBin peekNextBin() {
            if(binQueue.isEmpty())
                throw new NoSuchElementException("No more bins are available");
            BinQueueState current = binQueue.peek();
            return new ReaderBin(current.getReaderID(),current.getIndex(),current.getReferenceSequence(),current.peekNextBin());
        }

        private ReaderBin getNextBin() {
            if(binQueue.isEmpty())
                throw new NoSuchElementException("No more bins are available");
            BinQueueState current = binQueue.remove();
            ReaderBin readerBin = new ReaderBin(current.getReaderID(),current.getIndex(),current.getReferenceSequence(),current.nextBin());
            if(current.hasNextBin())
                binQueue.add(current);
            return readerBin;
        }

    }

    /**
     * Filters out bins not at the lowest level in the tree.
     */
    private static class LowestLevelBinFilteringIterator implements Iterator<Bin> {
        private BrowseableBAMIndex index;
        private Iterator<Bin> wrappedIterator;

        private Bin nextBin;

        public LowestLevelBinFilteringIterator(final BrowseableBAMIndex index, Iterator<Bin> iterator) {
            this.index = index;
            this.wrappedIterator = iterator;
            advance();
        }

        public boolean hasNext() {
            return nextBin != null;
        }

        public Bin next() {
            Bin bin = nextBin;
            advance();
            return bin;
        }

        public void remove() { throw new UnsupportedOperationException("Remove operation is not supported"); }

        private void advance() {
            nextBin = null;
            while(wrappedIterator.hasNext() && nextBin == null) {
                Bin bin = wrappedIterator.next();
                if(index.getLevelForBin(bin) == AbstractBAMFileIndex.getNumIndexLevels()-1)
                    nextBin = bin;
            }
        }
    }
}

class BinQueueState implements Comparable<org.broadinstitute.sting.gatk.datasources.reads.BinQueueState> {
    private final SAMReaderID id;
    private final BrowseableBAMIndex index;
    private final int referenceSequence;
    private final PeekableIterator<Bin> bins;

    private int firstLocusInCurrentBin;
    private int lastLocusInCurrentBin;

    public BinQueueState(final SAMReaderID id, final BrowseableBAMIndex index, final int referenceSequence, final Iterator<Bin> bins) {
        this.id = id;
        this.index = index;
        this.referenceSequence = referenceSequence;
        this.bins = new PeekableIterator<Bin>(bins);
        refreshLocusInBinCache();
    }

    public SAMReaderID getReaderID() {
        return id;
    }

    public BrowseableBAMIndex getIndex() {
        return index;
    }

    public int getReferenceSequence() {
        return referenceSequence;
    }

    public boolean hasNextBin() {
        return bins.hasNext();
    }

    public Bin peekNextBin() {
        return bins.peek();
    }

    public Bin nextBin() {
        Bin nextBin = bins.next();
        refreshLocusInBinCache();
        return nextBin;
    }

    public int compareTo(org.broadinstitute.sting.gatk.datasources.reads.BinQueueState other) {
        if(!this.bins.hasNext() && !other.bins.hasNext()) return 0;
        if(!this.bins.hasNext()) return -1;
        if(!this.bins.hasNext()) return 1;

        // Both BinQueueStates have next bins.  Before proceeding, make sure the bin cache is valid.
        if(this.firstLocusInCurrentBin <= 0 || this.lastLocusInCurrentBin <= 0 ||
           other.firstLocusInCurrentBin <= 0 || other.lastLocusInCurrentBin <= 0) {
            throw new ReviewedStingException("Sharding mechanism error - bin->locus cache is invalid.");
        }

        // Straight integer subtraction works here because lhsStart, rhsStart always positive.
        if(this.firstLocusInCurrentBin != other.firstLocusInCurrentBin)
            return this.firstLocusInCurrentBin - other.firstLocusInCurrentBin;

        // Straight integer subtraction works here because lhsStop, rhsStop always positive.
        return this.lastLocusInCurrentBin - other.lastLocusInCurrentBin;
    }

    private void refreshLocusInBinCache() {
        firstLocusInCurrentBin = -1;
        lastLocusInCurrentBin = -1;
        if(bins.hasNext()) {
            Bin bin = bins.peek();
            firstLocusInCurrentBin = index.getFirstLocusInBin(bin);
            lastLocusInCurrentBin = index.getLastLocusInBin(bin);
        }
    }
}
/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.datasources.shards;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.BlockDrivenSAMDataSource;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.SAMReaderID;

import java.util.*;

import net.sf.samtools.*;
import net.sf.picard.util.PeekableIterator;

/**
 * Shard intervals based on position within the BAM file.
 *
 * @author mhanna
 * @version 0.1
 */
public class IntervalSharder {
    public static Iterator<FilePointer> shardIntervals(final BlockDrivenSAMDataSource dataSource, final List<GenomeLoc> loci) {
        return new FilePointerIterator(dataSource,loci);
    }

    /**
     * A lazy-loading iterator over file pointers.
     */
    private static class FilePointerIterator implements Iterator<FilePointer> {
        final BlockDrivenSAMDataSource dataSource;
        final PeekableIterator<GenomeLoc> locusIterator;
        final Queue<FilePointer> cachedFilePointers = new LinkedList<FilePointer>();

        public FilePointerIterator(final BlockDrivenSAMDataSource dataSource, final List<GenomeLoc> loci) {
            this.dataSource = dataSource;
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
            List<GenomeLoc> nextBatch = new ArrayList<GenomeLoc>();
            String contig = null;

            while(locusIterator.hasNext() && nextBatch.isEmpty()) {
                contig = null;
                while(locusIterator.hasNext() && (contig == null || locusIterator.peek().getContig().equals(contig))) {
                    GenomeLoc nextLocus = locusIterator.next();
                    contig = nextLocus.getContig();
                    nextBatch.add(nextLocus);
                }
            }

            if(nextBatch.size() > 0)
                cachedFilePointers.addAll(shardIntervalsOnContig(dataSource,contig,nextBatch));
        }
    }
    
    private static List<FilePointer> shardIntervalsOnContig(final BlockDrivenSAMDataSource dataSource, final String contig, final List<GenomeLoc> loci) {
        // Gather bins for the given loci, splitting loci as necessary so that each falls into exactly one lowest-level bin.
        List<FilePointer> filePointers = new ArrayList<FilePointer>();
        FilePointer lastFilePointer = null;
        BAMOverlap lastBAMOverlap = null;

        Map<SAMReaderID,CachingBAMFileIndex> readerToIndexMap = new HashMap<SAMReaderID,CachingBAMFileIndex>();
        BinMergingIterator binMerger = new BinMergingIterator();
        for(SAMReaderID id: dataSource.getReaderIDs()) {
            final SAMSequenceRecord referenceSequence = dataSource.getHeader(id).getSequence(contig);
            if(referenceSequence == null)
                continue;
            final CachingBAMFileIndex index = dataSource.getIndex(id);
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
                throw new StingException("Location outside bounds of contig");

            if(!binIterator.hasNext())
                break;

            int locationStart = (int)location.getStart();
            final int locationStop = (int)location.getStop();

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
                    lastFilePointer = new FilePointer(contig,lastBAMOverlap);
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

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,regionStop);
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

                    GenomeLoc subset = GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,locationStop);
                    filePointers.add(new FilePointer(subset));

                    locationStart = locationStop + 1;
                }
                else {
                    if(lastFilePointer == null)
                        throw new StingException("Illegal state: initializer failed to create cached file pointer.");

                    // The start of the region overlaps the bin.  Add the overlapping subset.
                    final int regionStop = Math.min(locationStop,binStop);
                    lastFilePointer.addLocation(GenomeLocParser.createGenomeLoc(location.getContig(),locationStart,regionStop));
                    locationStart = regionStop + 1;
                }
            }
        }

        if(lastFilePointer != null && lastFilePointer.locations.size() > 0)
            filePointers.add(lastFilePointer);

        // Lookup the locations for every file pointer in the index.
        for(SAMReaderID id: readerToIndexMap.keySet()) {
            CachingBAMFileIndex index = readerToIndexMap.get(id);
            for(FilePointer filePointer: filePointers)
                filePointer.addFileSpans(id,index.getChunksOverlapping(filePointer.overlap.getBin(id)));
            index.close();
        }
        
        return filePointers;
    }

    private static class BinMergingIterator implements Iterator<BAMOverlap> {
        private PriorityQueue<BinQueueState> binQueue = new PriorityQueue<BinQueueState>();
        private Queue<BAMOverlap> pendingOverlaps = new LinkedList<BAMOverlap>();

        public void addReader(final SAMReaderID id, final CachingBAMFileIndex index, final int referenceSequence, Iterator<Bin> bins) {
            binQueue.add(new BinQueueState(id,index,referenceSequence,new LowestLevelBinFilteringIterator(index,bins)));
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
            return new ReaderBin(current.id,current.index,current.referenceSequence,current.bins.peek());
        }

        private ReaderBin getNextBin() {
            if(binQueue.isEmpty())
                throw new NoSuchElementException("No more bins are available");
            BinQueueState current = binQueue.remove();
            ReaderBin readerBin = new ReaderBin(current.id,current.index,current.referenceSequence,current.bins.next());
            if(current.bins.hasNext())
                binQueue.add(current);
            return readerBin;
        }

        private class ReaderBin {
            public final SAMReaderID id;
            public final CachingBAMFileIndex index;
            public final int referenceSequence;
            public final Bin bin;

            public ReaderBin(final SAMReaderID id, final CachingBAMFileIndex index, final int referenceSequence, final Bin bin) {
                this.id = id;
                this.index = index;
                this.referenceSequence = referenceSequence;
                this.bin = bin;
            }

            public int getStart() {
                return index.getFirstLocusInBin(bin);
            }

            public int getStop() {
                return index.getLastLocusInBin(bin);
            }
        }

        private class BinQueueState implements Comparable<BinQueueState> {
            public final SAMReaderID id;
            public final CachingBAMFileIndex index;
            public final int referenceSequence;
            public final PeekableIterator<Bin> bins;

            public BinQueueState(final SAMReaderID id, final CachingBAMFileIndex index, final int referenceSequence, final Iterator<Bin> bins) {
                this.id = id;
                this.index = index;
                this.referenceSequence = referenceSequence;
                this.bins = new PeekableIterator<Bin>(bins);
            }

            public int compareTo(BinQueueState other) {
                if(!this.bins.hasNext() && !other.bins.hasNext()) return 0;
                if(!this.bins.hasNext()) return -1;
                if(!this.bins.hasNext()) return 1;

                int thisStart = this.index.getFirstLocusInBin(this.bins.peek());
                int otherStart = other.index.getFirstLocusInBin(other.bins.peek());

                // Straight integer subtraction works here because lhsStart, rhsStart always positive.
                if(thisStart != otherStart)
                    return thisStart - otherStart;

                int thisStop = this.index.getLastLocusInBin(this.bins.peek());
                int otherStop = other.index.getLastLocusInBin(other.bins.peek());

                // Straight integer subtraction works here because lhsStop, rhsStop always positive.
                return thisStop - otherStop;
            }
        }
    }

    /**
     * Filters out bins not at the lowest level in the tree.
     */
    private static class LowestLevelBinFilteringIterator implements Iterator<Bin> {
        private CachingBAMFileIndex index;
        private Iterator<Bin> wrappedIterator;

        private Bin nextBin;

        public LowestLevelBinFilteringIterator(final CachingBAMFileIndex index, Iterator<Bin> iterator) {
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
                if(index.getLevelForBin(bin) == index.getNumIndexLevels()-1)
                    nextBin = bin;
            }
        }
    }    
}

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
class FilePointer {
    protected final Map<SAMReaderID,SAMFileSpan> fileSpans = new HashMap<SAMReaderID,SAMFileSpan>();
    protected final String referenceSequence;
    protected final BAMOverlap overlap;
    protected final List<GenomeLoc> locations;

    public FilePointer(final GenomeLoc location) {
        referenceSequence = location.getContig();
        overlap = null;
        locations = Collections.singletonList(location);
    }

    public FilePointer(final String referenceSequence,final BAMOverlap overlap) {
        this.referenceSequence = referenceSequence;
        this.overlap = overlap;
        this.locations = new ArrayList<GenomeLoc>();
    }

    public void addLocation(GenomeLoc location) {
        locations.add(location);
    }

    public void addFileSpans(SAMReaderID id, SAMFileSpan fileSpan) {
        this.fileSpans.put(id,fileSpan);
    }
}

/**
 * Models a bin at which all BAM files in the merged input stream overlap.
 */
class BAMOverlap {
    public final int start;
    public final int stop;

    private final Map<SAMReaderID,Bin> bins = new HashMap<SAMReaderID,Bin>();

    public BAMOverlap(final int start, final int stop) {
        this.start = start;
        this.stop = stop;
    }

    public void addBin(final SAMReaderID id, final Bin bin) {
        bins.put(id,bin);
    }

    public Bin getBin(final SAMReaderID id) {
        return bins.get(id);
    }
}




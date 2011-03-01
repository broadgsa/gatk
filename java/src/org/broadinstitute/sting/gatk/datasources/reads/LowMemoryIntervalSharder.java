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
import net.sf.samtools.GATKBin;
import net.sf.samtools.GATKChunk;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Assign intervals to the most appropriate blocks, keeping as little as possible in memory at once.
 */
public class LowMemoryIntervalSharder implements Iterator<FilePointer> {
    private static Logger logger = Logger.getLogger(IntervalSharder.class);

    private final SAMDataSource dataSource;

    private final GenomeLocSortedSet loci;

    private final PeekableIterator<GenomeLoc> locusIterator;

    private GenomeLoc currentLocus;

    private FilePointer nextFilePointer = null;

    public LowMemoryIntervalSharder(final SAMDataSource dataSource, final GenomeLocSortedSet loci) {
        this.dataSource = dataSource;
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
            nextFilePointer = new FilePointer(currentLocus.getContig());

            int coveredRegionStart = 1;
            int coveredRegionStop = Integer.MAX_VALUE;
            GenomeLoc coveredRegion = null;

            for(SAMReaderID reader: dataSource.getReaderIDs()) {
                GATKBAMIndex index = (GATKBAMIndex)dataSource.getIndex(reader);
                BinTree binTree = getNextOverlappingBinTree((GATKBAMIndex)dataSource.getIndex(reader),currentLocus);
                if(binTree != null) {
                    coveredRegionStart = Math.max(coveredRegionStart,binTree.getStart());
                    coveredRegionStop = Math.min(coveredRegionStop,binTree.getStop());
                    coveredRegion = loci.getGenomeLocParser().createGenomeLoc(currentLocus.getContig(),coveredRegionStart,coveredRegionStop);

                    GATKBAMFileSpan fileSpan = generateFileSpan(index,binTree,currentLocus);
                    nextFilePointer.addFileSpans(reader,fileSpan);
                }
            }

            // Early exit if no bins were found.
            if(coveredRegion == null) {
                nextFilePointer.addLocation(currentLocus);
                currentLocus = locusIterator.next();
                continue;
            }

            // Define the initial range of the file pointer, aka the region where the locus currently being processed intersects the BAM list.
            GenomeLoc initialLocation = currentLocus.intersect(coveredRegion);
            nextFilePointer.addLocation(initialLocation);

            // See whether the BAM regions discovered overlap the next set of intervals in the interval list.  If so, include every overlapping interval.
            if(!nextFilePointer.locations.isEmpty()) {
                while(locusIterator.hasNext() && locusIterator.peek().overlapsP(coveredRegion)) {
                    currentLocus = locusIterator.next();
                    nextFilePointer.addLocation(locusIterator.next());
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
    private int lastReferenceSequenceLoaded = -1;

    /**
     * The stateful iterator used to progress through the genoem.
     */
    private PeekableIterator<BinTree> binTreeIterator = null;

    /**
     * Get the next overlapping tree of bins associated with the given BAM file.
     * @param index BAM index representation.
     * @param locus Locus for which to grab the bin tree, if available.
     * @return The BinTree overlapping the given locus.
     */
    private BinTree getNextOverlappingBinTree(final GATKBAMIndex index, final GenomeLoc locus) {
        // Stale reference sequence or first invocation.  (Re)create the binTreeIterator.
        if(locus.getContigIndex() != lastReferenceSequenceLoaded) {
            if(binTreeIterator != null)
                binTreeIterator.close();
            lastReferenceSequenceLoaded = locus.getContigIndex();
            binTreeIterator = new PeekableIterator<BinTree>(new BinTreeIterator(index,index.getIndexFile(),locus.getContigIndex()));
        }

        if(!binTreeIterator.hasNext())
            return null;

        BinTree binTree = binTreeIterator.peek();
        while(binTree.isBefore(locus)) {
            binTreeIterator.next(); // Before the point of interest.  Consume this one.
            binTree = binTreeIterator.peek();
        }

        if(binTree.overlaps(locus)) {
            binTreeIterator.next();
            return binTree;
        }

        return null;
    }

    /**
     * Converts a bin list to a file span, trimmed based on the linear index and with overlapping regions removed.
     * @param index BAM index.
     * @param binTree Tree of data found to overlap the region.  binTree.overlaps(initialRegion) must return true.
     * @param initialRegion The region to employ when trimming the linear index.
     * @return File span mapping to given region.
     */
    private GATKBAMFileSpan generateFileSpan(final GATKBAMIndex index, final BinTree binTree, final GenomeLoc initialRegion) {
        List<GATKChunk> chunks = new ArrayList<GATKChunk>(binTree.size());
        for(GATKBin bin: binTree.getBins()) {
            if(bin == null)
                continue;
            // The optimizer below will mutate the chunk list.  Make sure each element is a clone of the reference sequence.
            for(GATKChunk chunk: bin.getChunkList())
                chunks.add(chunk.clone());
        }

        // Optimize the chunk list with a linear index optimization
        chunks = index.optimizeChunkList(chunks,index.getLinearIndex(initialRegion.getContigIndex()).getMinimumOffset(initialRegion.getStart()));
        
        GATKBAMFileSpan fileSpan = new GATKBAMFileSpan(chunks.toArray(new GATKChunk[chunks.size()]));

        return fileSpan;
    }
}

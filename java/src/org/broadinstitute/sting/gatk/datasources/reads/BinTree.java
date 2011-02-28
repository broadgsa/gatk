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
import net.sf.samtools.Bin;
import net.sf.samtools.GATKBin;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Represents a tree of overlapping bins in a single
 * BAM index.
 */
public class BinTree {
    /**
     * The BAM index from which this bin data is sourced.
     */
    private final GATKBAMIndex index;

    /**
     * The bins in this tree, organized by level.
     */
    private final GATKBin[] bins;

    public BinTree(GATKBAMIndex index,final GATKBin[] bins) {
        this.index = index;
        this.bins = bins;            
    }

    public GATKBin getLowestLevelBin() {
        for(int i = bins.length-1; i >= 0; i--) {
            if(bins[i] != null)
                return bins[i];
        }
        return null;
    }


    /**
     * Retrieve the bins from the bin tree.
     * @return list of bins.
     */
    public GATKBin[] getBins() {
        return bins;
    }

    /**
     * Gets the number of bins in a given bin list.
     * @return Number of bins in the list.
     */
    public int size() {
        return bins.length;
    }

    /**
     * Checks overlap between this bin tree and other bin trees.
     * @param position the position over which to detect overlap.
     * @return True if the segment overlaps.  False otherwise.
     */
    public boolean overlaps(final GenomeLoc position) {
        for(GATKBin gatkBin: bins) {
            if(gatkBin == null)
                continue;
            Bin bin = new Bin(gatkBin.getReferenceSequence(),gatkBin.getBinNumber());
            // Overlap occurs when the position is not disjoint with the bin boundaries.
            if(!(position.getStop() < index.getFirstLocusInBin(bin) || position.getStart() > index.getLastLocusInBin(bin)))
                return true;
        }
        return false;
    }
}

/**
 * Iterate through all bin trees in sequence, from those covering base 1 to those covering MAX_BINS.
 */
class BinTreeIterator implements Iterator<BinTree> {
    /**
     * The index over which to iterate.
     */
    private final GATKBAMIndex index;

    /**
     * Iterators over each individual level.
     */
    private final PeekableIterator<GATKBin>[] levelIterators;

    /**
     * The next bin tree to be returned.
     */
    private BinTree nextBinTree;

    /**
     * Each iteration through the bin tree has a corresponding lowest level.  Make sure
     * every lowest-level bin is covered, whether that bin is present or not.
     */
    private int currentBinInLowestLevel;

    public BinTreeIterator(final GATKBAMIndex index, final File indexFile, final int referenceSequence) {
        this.index = index;
        
        BAMIndexBinIterator binIterator = new BAMIndexBinIterator(index,indexFile,referenceSequence);
        levelIterators = new PeekableIterator[GATKBAMIndex.getNumIndexLevels()];
        for(int level = 0; level < GATKBAMIndex.getNumIndexLevels(); level++)
            levelIterators[level] = new PeekableIterator<GATKBin>(binIterator.getIteratorOverLevel(level));

        // Set the current bin to one less that the first bin in the sequence.  advance() will push it
        // ahead to the first bin in the lowest level.
        currentBinInLowestLevel = GATKBAMIndex.getFirstBinInLevel(GATKBAMIndex.getNumIndexLevels()-1) - 1;

        advance();
    }

    public boolean hasNext() {
        return nextBinTree != null;    
    }

    /**
     * Return the next BinTree in the level.
     * @return Next BinTree in sequence.
     */
    public BinTree next() {
        if(!hasNext())
            throw new NoSuchElementException("BinTreeIterator is out of elements");
        BinTree currentBinTree = nextBinTree;
        advance();
        return currentBinTree;
    }

    /**
     * Bring the bin tree ahead to the next overlapping structure.
     */
    private void advance() {
        final int lowestLevel = GATKBAMIndex.getNumIndexLevels()-1;
        final int firstBinInLowestLevel = GATKBAMIndex.getFirstBinInLevel(lowestLevel);
        final int binsInLowestLevel = index.getLevelSize(lowestLevel);

        GATKBin[] bins = new GATKBin[GATKBAMIndex.getNumIndexLevels()];
        nextBinTree = null;
        while(nextBinTree == null) {
            currentBinInLowestLevel++;
            boolean levelIteratorsExhausted = true;

            for(int level = lowestLevel; level >= 0; level--) {
                if(!levelIterators[level].hasNext())
                    continue;
                levelIteratorsExhausted = false;

                final int firstBinInThisLevel = GATKBAMIndex.getFirstBinInLevel(level);
                final int binsInThisLevel = index.getLevelSize(level);
                final int currentBinInThisLevel = ((currentBinInLowestLevel-firstBinInLowestLevel)*binsInThisLevel/binsInLowestLevel) + firstBinInThisLevel;

                while(levelIterators[level].hasNext() && levelIterators[level].peek().getBinNumber() < currentBinInThisLevel)
                    levelIterators[level].next();

                if(levelIterators[level].hasNext() && levelIterators[level].peek().getBinNumber() == currentBinInThisLevel)
                    bins[level] = levelIterators[level].peek();
            }

            // No more bins available for this reference sequence?  Break out of the loop.
            if(levelIteratorsExhausted)
                break;

            // Found a compelling bin tree?  Break out of the loop.
            for(int level = 0; level <= lowestLevel; level++) {
                if(bins[level] != null) {
                    nextBinTree = new BinTree(index,bins);
                    break;
                }
            }
        }
    }

    /**
     * Remove unsupported.
     */
    public void remove() {
       throw new UnsupportedOperationException("Cannot remove elements from a BinTreeIterator");
    }
}

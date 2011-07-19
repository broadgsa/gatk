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
import net.sf.samtools.SAMFileSpan;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

import java.util.*;

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
public class FilePointer {
    protected final SortedMap<SAMReaderID,SAMFileSpan> fileSpans = new TreeMap<SAMReaderID,SAMFileSpan>();
    protected final BAMOverlap overlap;
    protected final List<GenomeLoc> locations;

    /**
     * Does this file pointer point into an unmapped region?
     */
    protected final boolean isRegionUnmapped;

    public FilePointer() {
        this((BAMOverlap)null);
    }

    public FilePointer(final GenomeLoc location) {
        this.overlap = null;
        this.locations = Collections.singletonList(location);
        this.isRegionUnmapped = GenomeLoc.isUnmapped(location);
    }

    public FilePointer(final BAMOverlap overlap) {
        this.overlap = overlap;
        this.locations = new ArrayList<GenomeLoc>();
        this.isRegionUnmapped = false;
    }

    /**
     * Returns an immutable variant of the list of locations.
     * @return
     */
    public List<GenomeLoc> getLocations() {
        return Collections.unmodifiableList(locations);
    }

    @Override
    public boolean equals(final Object other) {
        if(!(other instanceof FilePointer))
            return false;
        FilePointer otherFilePointer = (FilePointer)other;

        // intervals
        if(this.locations.size() != otherFilePointer.locations.size())
            return false;
        for(int i = 0; i < locations.size(); i++) {
            if(!this.locations.get(i).equals(otherFilePointer.locations.get(i)))
                return false;
        }

        // fileSpans
        if(this.fileSpans.size() != otherFilePointer.fileSpans.size())
            return false;
        Iterator<Map.Entry<SAMReaderID,SAMFileSpan>> thisEntries = this.fileSpans.entrySet().iterator();
        Iterator<Map.Entry<SAMReaderID,SAMFileSpan>> otherEntries = otherFilePointer.fileSpans.entrySet().iterator();
        while(thisEntries.hasNext() || otherEntries.hasNext()) {
            if(!thisEntries.next().equals(otherEntries.next()))
                return false;
        }
        
        return true;
    }

    public void addLocation(final GenomeLoc location) {
        locations.add(location);
    }

    public void addFileSpans(final SAMReaderID id, final SAMFileSpan fileSpan) {
        this.fileSpans.put(id,fileSpan);
    }

    public void addFileSpans(final Map<SAMReaderID, GATKBAMFileSpan> fileSpans) {
        this.fileSpans.putAll(fileSpans);
    }


    /**
     * Computes the size of this file span, in uncompressed bytes.
     * @return Size of the file span.
     */
    public long size() {
        long size = 0L;
        for(SAMFileSpan fileSpan: fileSpans.values())
            size += ((GATKBAMFileSpan)fileSpan).size();
        return size;
    }

    /**
     * Returns the difference in size between two filespans.
     * @param other Other filespan against which to measure.
     * @return The difference in size between the two file pointers.
     */
    public long minus(final FilePointer other) {
        long difference = 0;
        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> thisIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(this.fileSpans.entrySet().iterator());
        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> otherIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(other.fileSpans.entrySet().iterator());

        while(thisIterator.hasNext()) {
            // If there are no elements left in the 'other' iterator, spin out this iterator.
            if(!otherIterator.hasNext()) {
                GATKBAMFileSpan nextSpan = (GATKBAMFileSpan)thisIterator.next().getValue();
                difference += nextSpan.size();
                continue;
            }

            // Otherwise, compare the latest value.
            int compareValue = thisIterator.peek().getKey().compareTo(otherIterator.peek().getKey());

            if(compareValue < 0) {
                // This before other.
                difference += ((GATKBAMFileSpan)thisIterator.next().getValue()).size();
            }
            else if(compareValue > 0) {
                // Other before this.
                difference += ((GATKBAMFileSpan)otherIterator.next().getValue()).size();
            }
            else {
                // equality; difference the values.
                GATKBAMFileSpan thisRegion = (GATKBAMFileSpan)thisIterator.next().getValue();
                GATKBAMFileSpan otherRegion = (GATKBAMFileSpan)otherIterator.next().getValue();
                difference += Math.abs(thisRegion.minus(otherRegion).size());
            }
        }
        return difference;
    }

    /**
     * Combines two file pointers into one.
     * @param parser The genomelocparser to use when manipulating intervals.
     * @param other File pointer to combine into this one.
     * @return A completely new file pointer that is the combination of the two.
     */
    public FilePointer combine(final GenomeLocParser parser, final FilePointer other) {
        FilePointer combined = new FilePointer();

        List<GenomeLoc> intervals = new ArrayList<GenomeLoc>();
        intervals.addAll(locations);
        intervals.addAll(other.locations);
        for(GenomeLoc interval: IntervalUtils.sortAndMergeIntervals(parser,intervals,IntervalMergingRule.ALL))
            combined.addLocation(interval);

        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> thisIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(this.fileSpans.entrySet().iterator());
        PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>> otherIterator = new PeekableIterator<Map.Entry<SAMReaderID,SAMFileSpan>>(other.fileSpans.entrySet().iterator());

        while(thisIterator.hasNext() || otherIterator.hasNext()) {
            int compareValue;
            if(!otherIterator.hasNext()) {
                compareValue = -1;
            }
            else if(!thisIterator.hasNext())
                compareValue = 1;
            else
                compareValue = thisIterator.peek().getKey().compareTo(otherIterator.peek().getKey());

            // This before other.
            if(compareValue < 0)
                mergeElementsInto(combined,thisIterator);
            // Other before this.
            else if(compareValue > 0)
                mergeElementsInto(combined,otherIterator);
            // equality; union the values.
            else
                mergeElementsInto(combined,thisIterator,otherIterator);
        }
        return combined;
    }

    /**
     * Roll the next element in the iterator into the combined entry.
     * @param combined Entry into which to roll the next element.
     * @param iterators Sources of next elements.
     */
    private void mergeElementsInto(final FilePointer combined, Iterator<Map.Entry<SAMReaderID,SAMFileSpan>>... iterators) {
        if(iterators.length == 0)
            throw new ReviewedStingException("Tried to add zero elements to an existing file pointer.");
        Map.Entry<SAMReaderID,SAMFileSpan> initialElement = iterators[0].next();
        GATKBAMFileSpan fileSpan = (GATKBAMFileSpan)initialElement.getValue();
        for(int i = 1; i < iterators.length; i++)
            fileSpan = fileSpan.union((GATKBAMFileSpan)iterators[i].next().getValue());
        combined.addFileSpans(initialElement.getKey(),fileSpan);
    }
}

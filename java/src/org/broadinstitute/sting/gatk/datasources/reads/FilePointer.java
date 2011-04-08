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
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.sting.utils.interval.IntervalMergingRule;
import org.broadinstitute.sting.utils.interval.IntervalUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Represents a small section of a BAM file, and every associated interval.
 */
class FilePointer {
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

        while(thisIterator.hasNext() || otherIterator.hasNext()) {
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
            int compareValue = thisIterator.peek().getKey().compareTo(otherIterator.peek().getKey());

            if(compareValue < 0) {
                // This before other.
                Map.Entry<SAMReaderID,SAMFileSpan> entry = thisIterator.next();
                combined.addFileSpans(entry.getKey(),entry.getValue());
            }
            else if(compareValue > 0) {
                // Other before this.
                Map.Entry<SAMReaderID,SAMFileSpan> entry = otherIterator.next();
                combined.addFileSpans(entry.getKey(),entry.getValue());
            }
            else {
                // equality; union the values.
                SAMReaderID reader = thisIterator.peek().getKey();
                GATKBAMFileSpan thisRegion = (GATKBAMFileSpan)thisIterator.next().getValue();
                GATKBAMFileSpan otherRegion = (GATKBAMFileSpan)otherIterator.next().getValue();
                combined.addFileSpans(reader,thisRegion.union(otherRegion));
            }
        }
        return combined;
    }
}

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

package org.broadinstitute.gatk.utils.downsampling;

import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;
import org.broadinstitute.gatk.utils.sam.ReadUtils;

import java.util.*;

/**
 * Utilities for using the downsamplers for common tasks
 *
 * User: depristo
 * Date: 3/6/13
 * Time: 4:26 PM
 */
public class DownsamplingUtils {
    private DownsamplingUtils() { }

    /**
     * Level the coverage of the reads in each sample to no more than downsampleTo reads, no reducing
     * coverage at any read start to less than minReadsPerAlignmentStart
     *
     * This algorithm can be used to handle the situation where you have lots of coverage in some interval, and
     * want to reduce the coverage of the big peak down without removing the many reads at the edge of this
     * interval that are in fact good
     *
     * This algorithm separately operates on the reads for each sample independently.
     *
     * @param reads a sorted list of reads
     * @param downsampleTo the targeted number of reads we want from reads per sample
     * @param minReadsPerAlignmentStart don't reduce the number of reads starting at a specific alignment start
     *                                  to below this.  That is, if this value is 2, we'll never reduce the number
     *                                  of reads starting at a specific start site to less than 2
     * @return a sorted list of reads
     */
    public static List<GATKSAMRecord> levelCoverageByPosition(final List<GATKSAMRecord> reads, final int downsampleTo, final int minReadsPerAlignmentStart) {
        if ( reads == null ) throw new IllegalArgumentException("reads must not be null");

        final List<GATKSAMRecord> downsampled = new ArrayList<GATKSAMRecord>(reads.size());

        final Map<String, Map<Integer, List<GATKSAMRecord>>> readsBySampleByStart = partitionReadsBySampleAndStart(reads);
        for ( final Map<Integer, List<GATKSAMRecord>> readsByPosMap : readsBySampleByStart.values() ) {
            final LevelingDownsampler<List<GATKSAMRecord>, GATKSAMRecord> downsampler = new LevelingDownsampler<List<GATKSAMRecord>, GATKSAMRecord>(downsampleTo, minReadsPerAlignmentStart);
            downsampler.submit(readsByPosMap.values());
            downsampler.signalEndOfInput();
            for ( final List<GATKSAMRecord> downsampledReads : downsampler.consumeFinalizedItems())
                downsampled.addAll(downsampledReads);
        }

        return ReadUtils.sortReadsByCoordinate(downsampled);
    }

    /**
     * Build the data structure mapping for each sample -> (position -> reads at position)
     *
     * Note that the map position -> reads isn't ordered in any meaningful way
     *
     * @param reads a list of sorted reads
     * @return a map containing the list of reads at each start location, for each sample independently
     */
    private static Map<String, Map<Integer, List<GATKSAMRecord>>> partitionReadsBySampleAndStart(final List<GATKSAMRecord> reads) {
        final Map<String, Map<Integer, List<GATKSAMRecord>>> readsBySampleByStart = new LinkedHashMap<String, Map<Integer, List<GATKSAMRecord>>>();

        for ( final GATKSAMRecord read : reads ) {
            Map<Integer, List<GATKSAMRecord>> readsByStart = readsBySampleByStart.get(read.getReadGroup().getSample());

            if ( readsByStart == null ) {
                readsByStart = new LinkedHashMap<Integer, List<GATKSAMRecord>>();
                readsBySampleByStart.put(read.getReadGroup().getSample(), readsByStart);
            }

            List<GATKSAMRecord> readsAtStart = readsByStart.get(read.getAlignmentStart());
            if ( readsAtStart == null ) {
                readsAtStart = new LinkedList<GATKSAMRecord>();
                readsByStart.put(read.getAlignmentStart(), readsAtStart);
            }

            readsAtStart.add(read);
        }

        return readsBySampleByStart;
    }
}

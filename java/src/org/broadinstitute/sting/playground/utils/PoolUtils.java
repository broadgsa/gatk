/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.utils;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.playground.gatk.walkers.poolseq.ReadOffsetQuad;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 27, 2009
 * Time: 12:31:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class PoolUtils {

    private PoolUtils() {
    }

    public static ReadOffsetQuad splitReadsByReadDirection(List<SAMRecord> reads, List<Integer> offsets) {
        ArrayList<SAMRecord> forwardReads;
        ArrayList<SAMRecord> reverseReads;
        ArrayList<Integer> forwardOffsets;
        ArrayList<Integer> reverseOffsets;

        if (reads == null) {
            forwardReads = null;
            reverseReads = null;
            forwardOffsets = null;
            reverseOffsets = null;
        } else {
            forwardReads = new ArrayList<SAMRecord>();
            reverseReads = new ArrayList<SAMRecord>();
            forwardOffsets = new ArrayList<Integer>();
            reverseOffsets = new ArrayList<Integer>();

            for (int readNo = 0; readNo < reads.size(); readNo++) {
                if (reads.get(readNo).getReadNegativeStrandFlag()) {
                    forwardReads.add(reads.get(readNo));
                    forwardOffsets.add(offsets.get(readNo));
                } else {
                    reverseReads.add(reads.get(readNo));
                    reverseOffsets.add(offsets.get(readNo));
                }
            }
        }

        return new ReadOffsetQuad(forwardReads,forwardOffsets,reverseReads,reverseOffsets);
    }


    public static Pair<List<SAMRecord>, List<Integer>> thresholdReadsByQuality(List<SAMRecord> reads, List<Integer> offsets, byte qThresh) {
        List<SAMRecord> threshReads;
        List<Integer> threshOffsets;
        if (reads == null) {
            threshReads = null;
            threshOffsets = null;
        } else if (qThresh <= 0) {
            threshReads = reads;
            threshOffsets = offsets;
        } else {
            threshReads = new ArrayList<SAMRecord>();
            threshOffsets = new ArrayList<Integer>();

            for (int readNo = 0; readNo < reads.size(); readNo++) {
                if (reads.get(readNo).getBaseQualities()[offsets.get(readNo)] >= qThresh) {
                    threshReads.add(reads.get(readNo));
                    threshOffsets.add(offsets.get(readNo));
                } // else do nothing
            }
        }

        return new Pair<List<SAMRecord>,List<Integer>>(threshReads, threshOffsets);
    }

    public static Pair<List<SAMRecord>,List<Integer>> thresholdReadsByMappingQuality( List<SAMRecord> reads, List<Integer> offsets, int mapQual ) {
        List<SAMRecord> goodMapReads;
        List<Integer> goodMapOffsets;
        if ( reads == null ) {
            goodMapReads = null;
            goodMapOffsets = null;
        } else if ( mapQual < 0 ) {
            goodMapReads = reads;
            goodMapOffsets = offsets;
        } else {
            goodMapReads = new ArrayList<SAMRecord>();
            goodMapOffsets = new ArrayList<Integer>();

            for ( int readNo = 0; readNo < reads.size(); readNo ++ ) {
                if ( reads.get(readNo).getMappingQuality() > mapQual ) {
                    goodMapReads.add(reads.get(readNo));
                    goodMapOffsets.add(offsets.get(readNo));
                }
            }
        }

        return new Pair<List<SAMRecord>,List<Integer>>(goodMapReads,goodMapOffsets);
    }


}

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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import net.sf.samtools.SAMRecord;

import java.io.PrintStream;

/**
 * This walker prints out quality score counts for first and second reads of a pair aggregated over all reads
 * in the interval.
 *
 * @Author: Chris Hartl
 */
public class PairedQualityScoreCountsWalker extends ReadWalker<Pair<byte[],Boolean>,Pair<CycleQualCounts,CycleQualCounts>> {
    @Output
    public PrintStream out;

    @Argument(fullName="readLength", shortName="rl", doc="Length of reads in the bam file", required=true)
    public int readLength = -1;

    public void initialize() { return; }

    public Pair<CycleQualCounts,CycleQualCounts> reduceInit() {
        return new Pair<CycleQualCounts,CycleQualCounts>( new CycleQualCounts(readLength), new CycleQualCounts(readLength) );
    }

    public Pair<CycleQualCounts,CycleQualCounts> reduce( Pair<byte[],Boolean> mapCounts, Pair<CycleQualCounts,CycleQualCounts> reduceCounts ) {
        if ( mapCounts != null ) {
            if ( mapCounts.second ) {
                reduceCounts.first.update(mapCounts.first);
            } else {
                reduceCounts.second.update(mapCounts.first);
            }
        }

        return reduceCounts;
    }

    public Pair<byte[],Boolean> map( ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        if ( canUseRead(read) ) {
            return getCorrectlyOrientedBaseQualities(read);
        } else {
            return null;
        }
    }

    private boolean canUseRead(SAMRecord read) {
        return ( ! read.getMateUnmappedFlag() && ! read.getReadUnmappedFlag() ) && ( read.getReadPairedFlag() && read.getReadLength() == readLength );
    }

    private Pair<byte[],Boolean> getCorrectlyOrientedBaseQualities(SAMRecord read) {
        byte[] quals = read.getReadNegativeStrandFlag() ? Utils.reverse(read.getBaseQualities()) : read.getBaseQualities();
        return new Pair<byte[], Boolean>(quals, read.getFirstOfPairFlag());
    }

    public void onTraversalDone(Pair<CycleQualCounts,CycleQualCounts> finalCounts) {
        StringBuilder output = new StringBuilder();
        output.append(String.format("%s\t%s\t%s%n","Cycle","First_read_counts","Second_read_counts"));
        for ( int offset = 0; offset < readLength; offset++ ) {
            output.append(String.format("%d\t%s\t%s%n",offset,finalCounts.first.getCountDistribution(offset),finalCounts.second.getCountDistribution(offset)));
        }
        out.printf("%s",output.toString());
    }

}

class CycleQualCounts {
    private long[][] qualityCountsByCycle;
    private int cycleLength;
    private int qualMax = QualityUtils.MAX_REASONABLE_Q_SCORE + 1;

    public CycleQualCounts(int cycleLength) {
        this.cycleLength = cycleLength;
        qualityCountsByCycle = new long[cycleLength][qualMax];
        for ( int cycle = 0; cycle < cycleLength; cycle++ ) {
            for ( int qual = 0; qual < qualMax; qual++) {
                qualityCountsByCycle[cycle][qual] = 0;
            }
        }
    }

    public void update(int offset, byte quality) {
        qualityCountsByCycle[offset][qualityToQualityIndex(quality)]++;
    }


    public void update(byte[] qualArray) {
        for ( int o = 0; o < cycleLength; o++ ) {
            update(o,qualArray[o]);
        }
    }

    private int qualityToQualityIndex(byte qual) {
        return qual < 0 ? 0 : qual > qualMax ? qualMax : qual;
    }

    public long[][] getCounts() { return qualityCountsByCycle; }

    public String getCountDistribution(int offset) {
        StringBuilder b = new StringBuilder();
        for ( int qual = 0; qual < qualMax-1; qual++ ) {
            b.append(String.format("%d;",qualityCountsByCycle[offset][qual]));    
        }
        b.append(String.format("%d",qualityCountsByCycle[offset][qualMax-1]));

        return b.toString();
    }
}
package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.NGSPlatform;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.EnumSet;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 30, 2009
 *
 * The Cycle covariate.
 * For Solexa the cycle is simply the position in the read (counting backwards if it is a negative strand read)
 * For 454 the cycle is the TACG flow cycle, that is, each flow grabs all the TACG's in order in a single cycle
 * For example, for the read: AAACCCCGAAATTTTTACTG
 * the cycle would be 11111111222333333344
 * For SOLiD the cycle is a more complicated mixture of ligation cycle and primer round
 */

public class CycleCovariate implements StandardCovariate {

    private static final int MAXIMUM_CYCLE_VALUE = 1000;

    private static final EnumSet<NGSPlatform> DISCRETE_CYCLE_PLATFORMS = EnumSet.of(NGSPlatform.ILLUMINA, NGSPlatform.SOLID, NGSPlatform.PACBIO, NGSPlatform.COMPLETE_GENOMICS);
    private static final EnumSet<NGSPlatform> FLOW_CYCLE_PLATFORMS = EnumSet.of(NGSPlatform.LS454, NGSPlatform.ION_TORRENT);

    // Initialize any member variables using the command-line arguments passed to the walkers
    @Override
    public void initialize(final RecalibrationArgumentCollection RAC) {
        if (RAC.DEFAULT_PLATFORM != null && !NGSPlatform.isKnown(RAC.DEFAULT_PLATFORM))
            throw new UserException.CommandLineException("The requested default platform (" + RAC.DEFAULT_PLATFORM + ") is not a recognized platform.");
    }

    // Used to pick out the covariate's value from attributes of the read
    @Override
    public void recordValues(final GATKSAMRecord read, final ReadCovariates values) {
        final int readLength = read.getReadLength();
        final NGSPlatform ngsPlatform = read.getNGSPlatform();

        // Discrete cycle platforms
        if (DISCRETE_CYCLE_PLATFORMS.contains(ngsPlatform)) {
            final int readOrderFactor = read.getReadPairedFlag() && read.getSecondOfPairFlag() ? -1 : 1;
            final int increment;
            int cycle;
            if (read.getReadNegativeStrandFlag()) {
                cycle = readLength * readOrderFactor;
                increment = -1 * readOrderFactor;
            }
            else {
                cycle = readOrderFactor;
                increment = readOrderFactor;
            }

            final int CUSHION = 4;
            final int MAX_CYCLE = readLength - CUSHION - 1;
            for (int i = 0; i < readLength; i++) {
                final int key = (i<CUSHION || i>MAX_CYCLE) ? -1 : keyFromCycle(cycle);
                values.addCovariate(key, key, key, i);
                cycle += increment;
            }
        }

        // Flow cycle platforms
        else if (FLOW_CYCLE_PLATFORMS.contains(ngsPlatform)) {

            final byte[] bases = read.getReadBases();

            // Differentiate between first and second of pair.
            // The sequencing machine cycle keeps incrementing for the second read in a pair. So it is possible for a read group
            // to have an error affecting quality at a particular cycle on the first of pair which carries over to the second of pair.
            // Therefore the cycle covariate must differentiate between first and second of pair reads.
            // This effect can not be corrected by pulling out the first of pair and second of pair flags into a separate covariate because
            //   the current sequential model would consider the effects independently instead of jointly.
            final boolean multiplyByNegative1 = read.getReadPairedFlag() && read.getSecondOfPairFlag();

            int cycle = multiplyByNegative1 ? -1 : 1;     // todo -- check if this is the right behavior for mate paired reads in flow cycle platforms.

            // BUGBUG: Consider looking at degradation of base quality scores in homopolymer runs to detect when the cycle incremented even though the nucleotide didn't change
            // For example, AAAAAAA was probably read in two flow cycles but here we count it as one
            if (!read.getReadNegativeStrandFlag()) { // Forward direction
                int iii = 0;
                while (iii < readLength) {
                    while (iii < readLength && bases[iii] == (byte) 'T') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    while (iii < readLength && bases[iii] == (byte) 'A') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    while (iii < readLength && bases[iii] == (byte) 'C') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    while (iii < readLength && bases[iii] == (byte) 'G') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }
                    if (iii < readLength) {
                        if (multiplyByNegative1)
                            cycle--;
                        else
                            cycle++;
                    }
                    if (iii < readLength && !BaseUtils.isRegularBase(bases[iii])) {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii++;
                    }

                }
            }
            else { // Negative direction
                int iii = readLength - 1;
                while (iii >= 0) {
                    while (iii >= 0 && bases[iii] == (byte) 'T') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    while (iii >= 0 && bases[iii] == (byte) 'A') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    while (iii >= 0 && bases[iii] == (byte) 'C') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    while (iii >= 0 && bases[iii] == (byte) 'G') {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                    if (iii >= 0) {
                        if (multiplyByNegative1)
                            cycle--;
                        else
                            cycle++;
                    }
                    if (iii >= 0 && !BaseUtils.isRegularBase(bases[iii])) {
                        final int key = keyFromCycle(cycle);
                        values.addCovariate(key, key, key, iii);
                        iii--;
                    }
                }
            }
        }

        // Unknown platforms
        else {
            throw new UserException("The platform (" + read.getReadGroup().getPlatform() + ") associated with read group " + read.getReadGroup() + " is not a recognized platform. Implemented options are e.g. illumina, 454, and solid");
        }
    }

    // Used to get the covariate's value from input csv file during on-the-fly recalibration
    @Override
    public final Object getValue(final String str) {
        return Integer.parseInt(str);
    }

    @Override
    public String formatKey(final int key) {
        int cycle = key >> 1;  // shift so we can remove the "sign" bit
        if ( (key & 1) != 0 )   // is the last bit set?
            cycle *= -1;        // then the cycle is negative
        return String.format("%d", cycle);
    }

    @Override
    public int keyFromValue(final Object value) {
        return (value instanceof String) ? keyFromCycle(Integer.parseInt((String) value)) : keyFromCycle((Integer) value);
    }

    @Override
    public int maximumKeyValue() {
        return (MAXIMUM_CYCLE_VALUE << 1) + 1;
    }

    private static int keyFromCycle(final int cycle) {
        // no negative values because values must fit into the first few bits of the long
        int result = Math.abs(cycle);
        result = result << 1; // shift so we can add the "sign" bit
        if ( cycle < 0 )
            result++;    // negative cycles get the lower-most bit set
        return result;
    }
}
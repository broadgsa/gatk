/*
 * Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.locusiterator;

import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.locusiterator.old.SAMRecordAlignmentState;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;

/**
 * Caliper microbenchmark of fragment pileup
 */
public class AlignmentStateMachinePerformance {
    final static int readLength = 101;
    final static int nReads = 10000;
    final static int locus = 1;

    public static void main(String[] args) {
        final int rep = Integer.valueOf(args[0]);
        final boolean useNew = Boolean.valueOf(args[1]);
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);

        int nIterations = 0;
        for ( final String cigar : Arrays.asList("101M", "50M10I40M", "50M10D40M") ) {
            for ( int j = 0; j < nReads; j++ ) {
                GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, readLength);
                read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
                final byte[] quals = new byte[readLength];
                for ( int i = 0; i < readLength; i++ )
                    quals[i] = (byte)(i % QualityUtils.MAX_QUAL_SCORE);
                read.setBaseQualities(quals);
                read.setCigarString(cigar);

                for ( int i = 0; i < rep; i++ ) {
                    if ( useNew ) {
                        final AlignmentStateMachine alignmentStateMachine = new AlignmentStateMachine(read);
                        while ( alignmentStateMachine.stepForwardOnGenome() != null ) {
                            nIterations++;
                        }
                    } else {
                        final SAMRecordAlignmentState alignmentStateMachine = new SAMRecordAlignmentState(read);
                        while ( alignmentStateMachine.stepForwardOnGenome() != null ) {
                            alignmentStateMachine.getRead();
                            nIterations++;
                        }
                    }
                }
            }
        }

        System.out.printf("iterations %d%n", nIterations);
    }
}

package org.broadinstitute.sting.gatk.walkers.recalibration;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.NGSPlatform;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.EnumSet;
import java.util.List;

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
 *  For Solexa the cycle is simply the position in the read (counting backwards if it is a negative strand read)
 *  For 454 the cycle is the TACG flow cycle, that is, each flow grabs all the TACG's in order in a single cycle
 *     For example, for the read: AAACCCCGAAATTTTTACTG
 *             the cycle would be 11111111222333333344
 *  For SOLiD the cycle is a more complicated mixture of ligation cycle and primer round
 */

public class CycleCovariate implements StandardCovariate {
    private final static EnumSet<NGSPlatform> DISCRETE_CYCLE_PLATFORMS = EnumSet.of(NGSPlatform.ILLUMINA, NGSPlatform.SOLID, NGSPlatform.PACBIO, NGSPlatform.COMPLETE_GENOMICS);
    private final static EnumSet<NGSPlatform> FLOW_CYCLE_PLATFORMS = EnumSet.of(NGSPlatform.LS454,  NGSPlatform.ION_TORRENT);

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        if( RAC.DEFAULT_PLATFORM != null ) {
            if( RAC.DEFAULT_PLATFORM.equalsIgnoreCase( "SLX" ) || RAC.DEFAULT_PLATFORM.equalsIgnoreCase( "ILLUMINA" ) ||
                RAC.DEFAULT_PLATFORM.contains( "454" ) || RAC.DEFAULT_PLATFORM.equalsIgnoreCase( "SOLID" ) || RAC.DEFAULT_PLATFORM.equalsIgnoreCase( "ABI_SOLID" ) ) {
                // nothing to do
            } else {
                throw new UserException.CommandLineException("The requested default platform (" + RAC.DEFAULT_PLATFORM +") is not a recognized platform. Implemented options are illumina, 454, and solid");
            }
        }
    }

    // Used to pick out the covariate's value from attributes of the read
    public void getValues(SAMRecord read, Comparable[] comparable) {

        //-----------------------------
        // Illumina, Solid, PacBio, and Complete Genomics
        //-----------------------------

        final NGSPlatform ngsPlatform = ((GATKSAMRecord)read).getNGSPlatform();
        if( DISCRETE_CYCLE_PLATFORMS.contains(ngsPlatform) ) {
            final int init;
            final int increment;
            if( !read.getReadNegativeStrandFlag() ) {
                // Differentiate between first and second of pair.
                // The sequencing machine cycle keeps incrementing for the second read in a pair. So it is possible for a read group
                // to have an error affecting quality at a particular cycle on the first of pair which carries over to the second of pair.
                // Therefore the cycle covariate must differentiate between first and second of pair reads.
                // This effect can not be corrected by pulling out the first of pair and second of pair flags into a separate covariate because
                //   the current sequential model would consider the effects independently instead of jointly.
                if( read.getReadPairedFlag() && read.getSecondOfPairFlag() ) {
                    //second of pair, positive strand
                    init = -1;
                    increment = -1;
                }
                else
                {
                    //first of pair, positive strand
                    init = 1;
                    increment = 1;
                }

            } else {
                if( read.getReadPairedFlag() && read.getSecondOfPairFlag() ) {
                    //second of pair, negative strand
                    init = -read.getReadLength();
                    increment = 1;
                }
                else
                {
                    //first of pair, negative strand
                    init = read.getReadLength();
                    increment = -1;
                }
            }

            int cycle = init;
            for(int i = 0; i < read.getReadLength(); i++) {
                comparable[i] = cycle;
                cycle += increment;
            }
        }

        //-----------------------------
        // 454 and Ion Torrent
        //-----------------------------
        else if( FLOW_CYCLE_PLATFORMS.contains(ngsPlatform) ) {

            final int readLength = read.getReadLength();
            final byte[] bases = read.getReadBases();

            // Differentiate between first and second of pair.
            // The sequencing machine cycle keeps incrementing for the second read in a pair. So it is possible for a read group
            // to have an error affecting quality at a particular cycle on the first of pair which carries over to the second of pair.
            // Therefore the cycle covariate must differentiate between first and second of pair reads.
            // This effect can not be corrected by pulling out the first of pair and second of pair flags into a separate covariate because
            //   the current sequential model would consider the effects independently instead of jointly.
            final boolean multiplyByNegative1 = read.getReadPairedFlag() && read.getSecondOfPairFlag();

            int cycle = multiplyByNegative1 ? -1 : 1;

            // BUGBUG: Consider looking at degradation of base quality scores in homopolymer runs to detect when the cycle incremented even though the nucleotide didn't change
            // For example, AAAAAAA was probably read in two flow cycles but here we count it as one
            if( !read.getReadNegativeStrandFlag() ) { // Forward direction
                int iii = 0;
                while( iii < readLength )
                {
                    while( iii < readLength && bases[iii] == (byte)'T' ) { comparable[iii] = cycle; iii++; }
                    while( iii < readLength && bases[iii] == (byte)'A' ) { comparable[iii] = cycle; iii++; }
                    while( iii < readLength && bases[iii] == (byte)'C' ) { comparable[iii] = cycle; iii++; }
                    while( iii < readLength && bases[iii] == (byte)'G' ) { comparable[iii] = cycle; iii++; }
                    if( iii < readLength ) { if (multiplyByNegative1) cycle--; else cycle++; }
                    if( iii < readLength && !BaseUtils.isRegularBase(bases[iii]) ) { comparable[iii] = cycle; iii++; }

                }
            } else { // Negative direction
                int iii = readLength-1;
                while( iii >= 0 )
                {
                    while( iii >= 0 && bases[iii] == (byte)'T' ) { comparable[iii] = cycle; iii--; }
                    while( iii >= 0 && bases[iii] == (byte)'A' ) { comparable[iii] = cycle; iii--; }
                    while( iii >= 0 && bases[iii] == (byte)'C' ) { comparable[iii] = cycle; iii--; }
                    while( iii >= 0 && bases[iii] == (byte)'G' ) { comparable[iii] = cycle; iii--; }
                    if( iii >= 0 ) { if (multiplyByNegative1) cycle--; else cycle++; }
                    if( iii >= 0 && !BaseUtils.isRegularBase(bases[iii]) ) { comparable[iii] = cycle; iii--; }
                }
            }
        }
        else  {
            throw new IllegalStateException("This method hasn't been implemented yet for " + read.getReadGroup().getPlatform());
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }
}
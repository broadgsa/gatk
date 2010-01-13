package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.BaseUtils;
import net.sf.samtools.SAMRecord;

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

    private static boolean warnedUserBadPlatform = false;
    private static String defaultPlatform;

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
        if( RAC.DEFAULT_PLATFORM.equalsIgnoreCase( "ILLUMINA" ) || RAC.DEFAULT_PLATFORM.contains( "454" ) || RAC.DEFAULT_PLATFORM.equalsIgnoreCase( "SOLID" ) ) {
            defaultPlatform = RAC.DEFAULT_PLATFORM;
        } else {
            throw new StingException( "The requested default platform (" + RAC.DEFAULT_PLATFORM +") is not a recognized platform. Implemented options are illumina, 454, and solid");
        }
    }

    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read, final int offset ) {

        int cycle = 1;

        //-----------------------------
        // ILLUMINA and SOLID
        //-----------------------------

        if( read.getReadGroup().getPlatform().equalsIgnoreCase( "ILLUMINA" ) || read.getReadGroup().getPlatform().equalsIgnoreCase( "SOLID" ) ) {
            cycle = offset + 1;
	        if( read.getReadNegativeStrandFlag() ) {
	            cycle = read.getReadLength() - offset;
	        }
        }

        //-----------------------------
        // 454
        //-----------------------------

        else if( read.getReadGroup().getPlatform().contains( "454" ) ) { // Some bams have "LS454" and others have just "454"
            final byte[] bases = read.getReadBases();

            // BUGBUG: Consider looking at degradation of base quality scores in homopolymer runs to detect when the cycle incremented even though the nucleotide didn't change
            // For example, AAAAAAA was probably read in two flow cycles but here we count it as one
            if( !read.getReadNegativeStrandFlag() ) { // Forward direction
                int iii = 0;
                while( iii <= offset )
                {
                    while( iii <= offset && bases[iii] == (byte)'T' ) { iii++; }
                    while( iii <= offset && bases[iii] == (byte)'A' ) { iii++; }
                    while( iii <= offset && bases[iii] == (byte)'C' ) { iii++; }
                    while( iii <= offset && bases[iii] == (byte)'G' ) { iii++; }
                    if( iii <= offset ) { cycle++; }
                    if( iii <= offset && !BaseUtils.isRegularBase(bases[iii]) ) { iii++; }

                }
            } else { // Negative direction
                int iii = bases.length-1;
                while( iii >= offset )
                {
                    while( iii >= offset && bases[iii] == (byte)'T' ) { iii--; }
                    while( iii >= offset && bases[iii] == (byte)'A' ) { iii--; }
                    while( iii >= offset && bases[iii] == (byte)'C' ) { iii--; }
                    while( iii >= offset && bases[iii] == (byte)'G' ) { iii--; }
                    if( iii >= offset ) { cycle++; }
                    if( iii >= offset && !BaseUtils.isRegularBase(bases[iii]) ) { iii--; }
                }
            }
        }

        //-----------------------------
        // SOLID (unused), only to be used in conjunction with PrimerRoundCovariate
        //-----------------------------

        //else if( read.getReadGroup().getPlatform().equalsIgnoreCase( "SOLID" ) ) {
        //    // The ligation cycle according to http://www3.appliedbiosystems.com/cms/groups/mcb_marketing/documents/generaldocuments/cms_057511.pdf
        //    int pos = offset + 1;
	    //    if( read.getReadNegativeStrandFlag() ) {
	    //        pos = read.getReadLength() - offset;
	    //    }
        //	cycle = pos / 5; // integer division
        //}

        //-----------------------------
        // UNRECOGNIZED PLATFORM
        //-----------------------------

        else { // Platform is unrecognized so revert to the default platform but warn the user first
        	if( !warnedUserBadPlatform ) {
                if( defaultPlatform != null) { // The user set a default platform
                    Utils.warnUser( "Platform string (" + read.getReadGroup().getPlatform() + ") unrecognized in CycleCovariate. " +
                            "Reverting to " + defaultPlatform + " definition of machine cycle." );
                } else { // The user did not set a default platform
                    Utils.warnUser( "Platform string (" + read.getReadGroup().getPlatform() + ") unrecognized in CycleCovariate. " +
                            "Reverting to Illumina definition of machine cycle. Users may set the default platform using the --default_platform <String> argument." );
                    defaultPlatform = "Illumina";
                }
                warnedUserBadPlatform = true;
            }
            read.getReadGroup().setPlatform( defaultPlatform );
            return getValue( read, offset ); // A recursive call
        }

        // This functionality was moved to PairedReadOrderCovariate
        // Differentiate between first and second of pair
        //if( read.getReadPairedFlag() && read.getSecondOfPairFlag() ) {
        //    cycle *= -1;
        //}

        return cycle;
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }
}
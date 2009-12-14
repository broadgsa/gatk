package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
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
 *  For 454 the cycle is the number of discontinuous nucleotides seen during the length of the read
 *     For example, for the read: AAACCCCGAAATTTTTTT
 *             the cycle would be 111222234445555555
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

        //-----------------------------
        // ILLUMINA
        //-----------------------------

        if( read.getReadGroup().getPlatform().equalsIgnoreCase( "ILLUMINA" ) ) {
            int cycle = offset;
	        if( read.getReadNegativeStrandFlag() ) {
	            cycle = read.getReadLength() - (offset + 1);
	        }
	        return cycle;
        }

        //-----------------------------
        // 454
        //-----------------------------

        else if( read.getReadGroup().getPlatform().contains( "454" ) ) { // Some bams have "LS454" and others have just "454"
            int cycle = 0;
            byte[] bases = read.getReadBases();
            if( !read.getReadNegativeStrandFlag() ) { // forward direction
                byte prevBase = bases[0];
                for( int iii = 1; iii <= offset; iii++ ) {
                    if( bases[iii] != prevBase ) { // This base doesn't match the previous one so it is a new cycle
                        cycle++;
                        prevBase = bases[iii];
                    }
                }
            } else { // negative direction
                byte prevBase = bases[bases.length-1];
                for( int iii = bases.length-2; iii >= offset; iii-- ) {
                    if( bases[iii] != prevBase ) { // This base doesn't match the previous one so it is a new cycle
                        cycle++;
                        prevBase = bases[iii];
                    }
                }
            }
            return cycle;
        }

        //-----------------------------
        // SOLID
        //-----------------------------

        else if( read.getReadGroup().getPlatform().equalsIgnoreCase( "SOLID" ) ) {
            // The ligation cycle according to http://www3.appliedbiosystems.com/cms/groups/mcb_marketing/documents/generaldocuments/cms_057511.pdf
            int pos = offset;
	        if( read.getReadNegativeStrandFlag() ) {
	            pos = read.getReadLength() - (offset + 1);
	        }
        	return pos / 5; // integer division
        }

        //-----------------------------
        // UNRECOGNIZED PLATFORM
        //-----------------------------

        else { // Platform is unrecognized so revert to the default platform but warn the user first
        	if( !warnedUserBadPlatform ) {
                if( defaultPlatform != null) { // the user set a default platform
                    Utils.warnUser( "Platform string (" + read.getReadGroup().getPlatform() + ") unrecognized in CycleCovariate. " +
                            "Reverting to " + defaultPlatform + " definition of machine cycle." );
                } else { // the user did not set a default platform
                    Utils.warnUser( "Platform string (" + read.getReadGroup().getPlatform() + ") unrecognized in CycleCovariate. " +
                            "Reverting to Illumina definition of machine cycle. Users may set the default platform using the --default_platform <String> argument." );
                    defaultPlatform = "Illumina";
                }
                warnedUserBadPlatform = true;
            }
            read.getReadGroup().setPlatform( defaultPlatform );
            return getValue( read, offset ); // a recursive call
        }
    }

    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }

    // Used to estimate the amount space required for the full data HashMap
    public final int estimatedNumberOfBins() {
        return 100;
    }
}
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

package org.broadinstitute.gatk.tools.walkers.qc;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMReadGroupRecord;
import org.broadinstitute.gatk.utils.commandline.Advanced;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.AlignmentUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Collect read clipping statistics
 *
 * <p>This tool collects statistics about the read length, number of clipping events, and length
 * of the clipping in all reads in the dataset.</p>
 *
 * <h3>Input</h3>
 * One or more BAM files.
 *
 * <h3>Output</h3>
 * A simple tabulated text file with read length and clipping statistics for every read (or every given number of reads
 * if the "skip" option is used).
 *
 * <h3>Caveat</h3>
 * <p>This tool ignores "N" events in the CIGAR string.</p>
 *
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS})
public class ReadClippingStats extends ReadWalker<ReadClippingStats.ReadClippingInfo,Integer> {
    @Output
    protected PrintStream out;

    /**
     * when this flag is set (default), statistics will be collected on unmapped reads as well. The default behavior
     * is to ignore unmapped reads."
     */
    @Argument(fullName="include_unmapped", shortName="u", doc="Include unmapped reads in the analysis", required=false)
    protected boolean INCLUDE_UNMAPPED = false;

    /**
     * print every read whose read number is divisible by SKIP. READ_NUMBER % SKIP == 0. First read in the file has read number = 1,
     * second is 2, third is 3, ... A value of 1 means print every read. A value of 1000 means print every 1000th read.
     */
    @Advanced
    @Argument(fullName="skip", shortName="skip", doc="Do not print all reads, skip some.", required=false)
    protected int SKIP = 1;

    public class ReadClippingInfo {
        SAMReadGroupRecord rg;
        int readLength, nClippingEvents, nClippedBases;
    }

    public ReadClippingInfo map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker metaDataTracker) {
        if ( AlignmentUtils.isReadUnmapped(read) && !INCLUDE_UNMAPPED)
            return null;

        ReadClippingInfo info = new ReadClippingInfo();
        info.rg = read.getReadGroup();

        if ( info.rg == null ) throw new UserException.ReadMissingReadGroup(read);

        for ( CigarElement elt : read.getCigar().getCigarElements() ) {
            switch ( elt.getOperator()) {
                case H : // ignore hard clips
                case S : // soft clip
                    info.nClippingEvents++;
                    info.nClippedBases += elt.getLength();
                    break;
                case EQ : // sequence match
                case M : // alignment match
                case D : // deletion w.r.t. the reference
                case P : // ignore pads
                case I : // insertion w.r.t. the reference
                case N : // reference skip (looks and gets processed just like a "deletion", just different logical meaning)
                case X : // sequence mismatch
                    break;
                default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + elt.getOperator());
            }
            info.readLength = read.getReadLength();
        }

        return info;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        out.println(Utils.join(" \t", Arrays.asList("ReadGroup", "ReadLength", "NClippingEvents", "NClippedBases", "PercentClipped")));
        return 0;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param info  result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(ReadClippingInfo info, Integer sum) {
        if ( info != null ) {
            if ( sum % SKIP == 0 ) {
                String id = info.rg.getReadGroupId();
                out.printf("%s\t %d\t %d\t %d\t %.2f%n",
                        id, info.readLength, info.nClippingEvents, info.nClippedBases,
                        100.0 * MathUtils.ratio(info.nClippedBases, info.readLength));
            }
            return sum + 1;
        } else {
            return sum;
        }
    }

}
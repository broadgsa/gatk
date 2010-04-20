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

package org.broadinstitute.sting.gatk.walkers.fasta;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.WalkerName;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.commandline.Argument;
import net.sf.samtools.SAMRecord;

import java.io.*;

/**
 * Converts reads from the input BAM file into fastq format, stripping away all alignment information.
 * Optionally re-reverses the negative strand alignments using the --re-reverse command-line argument.
 */
@WalkerName("BamToFastq")
public class BamToFastqWalker extends ReadWalker<Integer, Integer> {

    @Argument(fullName="re_reverse", shortName="reverse", doc="re-reverse bases and quals of reads from the negative strand", required=false)
    private Boolean RE_REVERSE = false;

    @Argument(fullName="SQFile", shortName="SQ", doc="Output path for secondary quality map (readName => SAM SQ field)", required=false)
    String SQFile = null;

    private PrintWriter sqbw = null;

    public void initialize() {
        if ( SQFile != null ) {
            try {
                sqbw = new PrintWriter(SQFile);
            } catch (IOException e) {
                throw new RuntimeException("Unable to open sq output file: " + SQFile);
            }
        }
    }

	public Integer map(char[] ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        out.println("@" + read.getReadName());
        if ( !RE_REVERSE || !read.getReadNegativeStrandFlag() ) {
            out.println(read.getReadString());
            out.println("+");
            out.println(read.getBaseQualityString());
        } else {
            out.println(BaseUtils.simpleReverseComplement(read.getReadString()));
            out.println("+");
            out.println(BaseUtils.reverse(read.getBaseQualityString()));
        }

        if ( sqbw != null ) {
            byte[] sqs = (byte[])read.getAttribute("SQ");
            if ( sqs != null ) {
                sqbw.print(read.getReadName() + "\t" + "SQ:H:");
                for ( byte sq : sqs ) {
                    sqbw.printf("%02X", sq);
                }
                sqbw.println();
            }
        }

        return 1;
	}

    public Integer reduceInit() {
        return 0;
    }

	public Integer reduce(Integer value, Integer sum) {
		return value + sum;
	}

    public void onTraversalDone(Integer sum) {
        if ( sqbw != null )
            sqbw.close();
        logger.info("Number of reads converted: " + sum);
    }
}
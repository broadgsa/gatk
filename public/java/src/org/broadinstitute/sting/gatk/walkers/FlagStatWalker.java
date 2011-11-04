package org.broadinstitute.sting.gatk.walkers;

import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;


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
 * A reimplementation of the 'samtools flagstat' subcommand in the GATK.  Walks
 * over all input data, accumulating statistics such as total number of reads,
 * reads with QC failure flag set, number of duplicates, percentage mapped, etc.
 * @author aaron
 */
@Requires({DataSource.READS})
public class FlagStatWalker extends ReadWalker<Integer, Integer> {
    @Output
    PrintStream out;

    // what comes out of the flagstat
    static class FlagStat {
        long readCount = 0L;
        long QC_failure = 0L;
        long duplicates = 0L;
        long mapped = 0L;
        long paired_in_sequencing = 0L;
        long read1 = 0L;
        long read2 = 0L;
        long properly_paired = 0L;
        long with_itself_and_mate_mapped = 0L;
        long singletons = 0L;
        long with_mate_mapped_to_a_different_chr = 0L;
        long with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 = 0L;

        public String toString() {
            String ret = "";
            StringBuilder builder = new StringBuilder(ret);
            NumberFormat percentFormatter = new DecimalFormat("#0.00");
            builder.append(readCount);
            builder.append(" in total\n");

            builder.append(QC_failure);
            builder.append(" QC failure\n");


            builder.append(duplicates);
            builder.append(" duplicates\n");

            builder.append(mapped);
            builder.append(" mapped (");
            builder.append(percentFormatter.format(( (float)mapped / (float)readCount ) * 100.0));
            builder.append("%)\n");

            builder.append(paired_in_sequencing);
            builder.append(" paired in sequencing\n");


            builder.append(read1);
            builder.append(" read1\n");

            builder.append(read2);
            builder.append(" read2\n");

            builder.append(properly_paired);
            builder.append(" properly paired (");
            builder.append(percentFormatter.format(( (float)properly_paired / (float)readCount ) * 100.0));
            builder.append("%)\n");


            builder.append(with_itself_and_mate_mapped);
            builder.append(" with itself and mate mapped\n");

            builder.append(singletons);
            builder.append(" singletons (");
            builder.append(percentFormatter.format(( (float)singletons / (float)readCount ) * 100.0));
                        builder.append("%)\n");


            builder.append(with_mate_mapped_to_a_different_chr);
            builder.append(" with mate mapped to a different chr\n");

            builder.append(with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5);
            builder.append(" with mate mapped to a different chr (mapQ>=5)");

            return builder.toString();
        }

    }


    private FlagStat myStat = new FlagStat();

    public Integer map( ReferenceContext ref, GATKSAMRecord read, ReadMetaDataTracker metaDataTracker ) {
        myStat.readCount++;
        if (read.getReadFailsVendorQualityCheckFlag()) {
            myStat.QC_failure++;
        }
        if (read.getDuplicateReadFlag()) {
            myStat.duplicates++;
        }
        if (read.getReferenceIndex() >= 0) {
            myStat.mapped++;
        }
        if (read.getReadPairedFlag()) {
            myStat.paired_in_sequencing++;

            if (read.getSecondOfPairFlag()) {
                myStat.read2++;
            } else if (read.getReadPairedFlag()) {
                myStat.read1++;
            }
            if (read.getProperPairFlag()) {

                myStat.properly_paired++;
            }
            if (!read.getMateUnmappedFlag() && read.getReferenceIndex() >= 0) {
                myStat.with_itself_and_mate_mapped++;
            }
            if (read.getMateUnmappedFlag()) {
                myStat.singletons++;
            }
        }
        if (read.getReferenceIndex() >= 0 && read.getMateReferenceIndex() >= 0 && ! read.getReferenceIndex().equals(read.getMateReferenceIndex())) {
            myStat.with_mate_mapped_to_a_different_chr++;

            if (read.getMappingQuality() >= 5) {
                myStat.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5++;
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

    public void onTraversalDone(Integer result) {
        //out.println("[REDUCE RESULT] Traversal result is: " + result);
        out.println(myStat.toString());
    }
}
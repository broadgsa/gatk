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

import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.NanoSchedulable;
import org.broadinstitute.gatk.engine.walkers.ReadWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.NumberFormat;

/**
 * Collect statistics about sequence reads based on their SAM flags
 *
 * <p>This tool emulates the behavior of 'samtools flagstat'. It collects statistics such as total number of reads,
 * reads with QC failure flag set, number of duplicates, percentage mapped, etc.</p>
 *
 * <h3>Input</h3>
 * <p>
 * A BAM file containing the sequence data.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Resulting stats are written to file if an output file name is given (with -o), otherwise output to stdout.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T FlagStat \
 *   -R reference.fasta \
 *   -I reads.bam \
 *   [-o output.txt]
 * </pre>
 *
 * @author aaron
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS})
public class FlagStat extends ReadWalker<FlagStat.FlagStatus, FlagStat.FlagStatus> implements NanoSchedulable {
    @Output
    PrintStream out;

    // what comes out of the flagstat
    public final static class FlagStatus {
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

        public FlagStatus add(final FlagStatus other) {
            readCount += other.readCount;
            QC_failure += other.QC_failure;
            duplicates += other.duplicates;
            mapped += other.mapped;
            paired_in_sequencing += other.paired_in_sequencing;
            read1 += other.read1;
            read2 += other.read2;
            properly_paired += other.properly_paired;
            with_itself_and_mate_mapped += other.with_itself_and_mate_mapped;
            singletons += other.singletons;
            with_mate_mapped_to_a_different_chr += other.with_mate_mapped_to_a_different_chr;
            with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5 += other.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5;

            return this;
        }

        public FlagStatus add(final GATKSAMRecord read) {
            this.readCount++;

            if (read.getReadFailsVendorQualityCheckFlag()) {
                this.QC_failure++;
            }
            if (read.getDuplicateReadFlag()) {
                this.duplicates++;
            }
            if (!read.getReadUnmappedFlag()) {
                this.mapped++;
            }
            if (read.getReadPairedFlag()) {
                this.paired_in_sequencing++;

                if (read.getSecondOfPairFlag()) {
                    this.read2++;
                } else if (read.getReadPairedFlag()) {
                    this.read1++;
                }
                if (read.getProperPairFlag()) {
                    this.properly_paired++;
                }
                if (!read.getReadUnmappedFlag() && !read.getMateUnmappedFlag()) {
                    this.with_itself_and_mate_mapped++;

                    if (!read.getReferenceIndex().equals(read.getMateReferenceIndex())) {
                        this.with_mate_mapped_to_a_different_chr++;

                        if (read.getMappingQuality() >= 5) {
                            this.with_mate_mapped_to_a_different_chr_maq_greaterequal_than_5++;
                        }
                    }
                }
                if (!read.getReadUnmappedFlag() && read.getMateUnmappedFlag()) {
                    this.singletons++;
                }
            }

            return this;
        }
    }


    @Override
    public FlagStatus map( final ReferenceContext ref, final GATKSAMRecord read, final RefMetaDataTracker metaDataTracker ) {
        return new FlagStatus().add(read);
   }

    @Override
    public FlagStatus reduceInit() {
        return new FlagStatus();
    }

    @Override
    public FlagStatus reduce(final FlagStatus value, final FlagStatus sum) {
        return sum.add(value);
    }

    @Override
    public void onTraversalDone(final FlagStatus result) {
        out.println(result.toString());
    }
}
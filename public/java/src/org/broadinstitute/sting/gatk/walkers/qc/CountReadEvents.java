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

package org.broadinstitute.sting.gatk.walkers.qc;

import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Walks over the input data set, counting the number of read events (from the CIGAR operator)
 *
 * <h3>Input</h3>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Number of read events for each category, formatted as a GATKReport table.
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -T CountReadEvents \
 *   -R ref.fasta \
 *   -I input.bam \
 *   -o output.grp \
 *   [-L input.intervals]
 * </pre>
 */

@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountReadEvents extends ReadWalker<Map<CigarOperator, ArrayList<Integer>> , Map<Integer, Map<CigarOperator, Long>>> {
    @Output
    PrintStream out;

    public Map<CigarOperator, ArrayList<Integer>> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
        return ReadUtils.getCigarOperatorForAllBases(read);
    }

    public Map<Integer, Map<CigarOperator, Long>> reduceInit() {
        return new HashMap<Integer, Map<CigarOperator, Long>>();
    }

    public Map<Integer, Map<CigarOperator, Long>> reduce(Map<CigarOperator, ArrayList<Integer>> value, Map<Integer, Map<CigarOperator, Long>> sum) {
        for (Map.Entry<CigarOperator, ArrayList<Integer>> entry : value.entrySet()) {
            CigarOperator op = entry.getKey();
            ArrayList<Integer> positions = entry.getValue();

            for (int p : positions) {
                Map<CigarOperator, Long> operatorCount = sum.get(p);
                if (operatorCount == null) {
                    operatorCount = new HashMap<CigarOperator, Long>();
                    sum.put(p, operatorCount);
                }

                Long count = operatorCount.get(op);
                if (count == null)
                    count = 0L;
                count++;
                operatorCount.put(op, count);
            }
        }
        return sum;
    }

    @Override
    public void onTraversalDone(Map<Integer, Map<CigarOperator, Long>> result) {
        GATKReport report = GATKReport.newSimpleReport("Events", "Position", "Event", "Observations");
        for (Map.Entry<Integer, Map<CigarOperator, Long>> entry : result.entrySet()) {
            int position = entry.getKey();
            Map<CigarOperator, Long> operatorCount = entry.getValue();

            for (Map.Entry<CigarOperator, Long> subEntry: operatorCount.entrySet()) {
                String operator = subEntry.getKey().name();
                Long observations = subEntry.getValue();
                report.addRow(position, operator, observations);
            }
        }
        report.print(out);
    }
}

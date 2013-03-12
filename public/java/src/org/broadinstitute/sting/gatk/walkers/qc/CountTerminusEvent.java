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

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

import java.io.PrintStream;
import java.util.List;

/**
 * Walks over the input data set, counting the number of reads ending in insertions/deletions or soft-clips
 *
 * <h3>Input</h3>
 * <p>
 * One or more BAM files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * Number of reads ending in each category.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CountTerminusEvent \
 *   -o output.txt \
 *   -I input.bam \
 *   [-L input.intervals]
 * </pre>
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Requires({DataSource.READS, DataSource.REFERENCE})
public class CountTerminusEvent extends ReadWalker<Pair<Long, Long>, Pair<Long, Long>> {
    @Output
    public PrintStream out;

    public Pair<Long, Long> map(ReferenceContext ref, GATKSAMRecord read, RefMetaDataTracker tracker) {
        List<CigarElement> cigarElements = read.getCigar().getCigarElements();

        CigarElement lastElement = null;
        for (CigarElement element : cigarElements) {
            if (element.getOperator() != CigarOperator.HARD_CLIP)
                lastElement = element;
        }

        if (lastElement == null)
            throw new UserException.MalformedBAM(read, "read does not have any bases, it's all hard clips");

        long endsInIndel = lastElement.getOperator() == CigarOperator.INSERTION || lastElement.getOperator() == CigarOperator.DELETION? 1 : 0;
        long endsInSC = lastElement.getOperator() == CigarOperator.SOFT_CLIP ? 1 : 0;

        return new Pair<Long, Long>(endsInIndel, endsInSC);
    }

    public Pair<Long, Long> reduceInit() { return new Pair<Long, Long>(0L, 0L); }

    public Pair<Long, Long> reduce(Pair<Long, Long> value, Pair<Long, Long> sum) {
        sum.set(sum.getFirst() + value.getFirst(), sum.getSecond() + value.getSecond());
        return sum;
    }

    @Override
    public void onTraversalDone(Pair<Long, Long> result) {
        out.println(String.format("\tReads ending in indels : %d\n\tReads ending in soft-clips: %d\n", result.getFirst(), result.getSecond()));
    }
}

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

package org.broadinstitute.gatk.tools.walkers.diagnostics;


import org.broadinstitute.gatk.engine.walkers.By;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;


import java.io.*;
import java.util.Collection;

/**
 * Report well-covered intervals
 *
 * <p>
 * This tool evaluates whether sites are well-covered or not according to specific coverage quality parameters, and
 * outputs a list of intervals that are considered well-covered, i.e. where most samples have good coverage. This is
 * useful for masking out poorly-covered sites where we cannot expect meaningful results in downstream analyses. See
 * argument defaults for what constitutes "most" samples and "good" coverage.
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant file and optionally, minimum coverage and sample percentage values.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * An list of well-covered intervals.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CoveredByNSamplesSites \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -out output.intervals \
 *   -minCov 15
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@By(DataSource.REFERENCE_ORDERED_DATA)
public class CoveredByNSamplesSites extends RodWalker<GenomeLoc, Integer> implements TreeReducible<Integer> {

    @Output(fullName = "OutputIntervals", shortName = "out", doc = "Name of file for output intervals")
    PrintStream outputStream;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(fullName = "minCoverage", shortName = "minCov",doc = "only samples that have coverage bigger than minCoverage will be counted",required = false)
    int minCoverage = 10;

    @Argument(fullName = "percentageOfSamples", shortName = "percentage", doc = "only sites where at least percentageOfSamples of the samples have good coverage, will be emitted", required = false)
    double percentageOfSamples = 0.9;

    @Override
    public GenomeLoc map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        if ( VCs.size() == 0 )
            return null;

        boolean emitSite = false;
        for(VariantContext vc : VCs){
            int coveredSamples = 0;
            final GenotypesContext genotypes = vc.getGenotypes();
            final int numOfGenotypes = genotypes.size();
            for(Genotype g : genotypes){
                if(g.getDP() >= minCoverage)
                    coveredSamples++;
            }
            if((double)coveredSamples/numOfGenotypes > percentageOfSamples){
                emitSite = true;
            }
        }
        if (emitSite)
            return ref.getLocus();
        else
            return null;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(GenomeLoc value, Integer sum) {
        if ( value != null ) {
            outputStream.println(value);
            sum++;
        }
        return sum;
    }

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    /**
     *
     * @param result the number of sites that passed the filter.
     */
    public void onTraversalDone(Integer result) {
        logger.info(result+" sites that have "+(percentageOfSamples*100)+"% of the samples with at least "+minCoverage+" coverage.\n");
    }



}

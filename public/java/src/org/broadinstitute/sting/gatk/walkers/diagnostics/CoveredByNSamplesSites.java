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

package org.broadinstitute.sting.gatk.walkers.diagnostics;


import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;


import java.io.*;
import java.util.Collection;

/**
 * print intervals file with all the variant sites that have "most" ( >= 90% by default) of the samples with "good" (>= 10 by default)coverage ("most" and "good" can be set in the command line).
 *
 * <p>
 * CoveredByNSamplesSites is a GATK tool for filter out sites based on their coverage.
 * The sites that pass the filter are printed out to an intervals file.
 *
 * <h2>Input</h2>
 * <p>
 * A variant file and optionally min coverage and sample percentage values.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An intervals file.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CoveredByNSamplesSites \
 *   -V input.vcf \
 *   -out output.intervals \
 *   -minCov 15
 * </pre>
 *
 */

@By(DataSource.REFERENCE_ORDERED_DATA)
public class CoveredByNSamplesSites extends RodWalker<GenomeLoc, Integer> implements TreeReducible<Integer> {

    @Output(fullName = "OutputIntervals", shortName = "out", doc = "Name of file for output intervals", required = true)
    PrintStream outputStream;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Argument(fullName = "minCoverage", shortName = "minCov",doc = "only samples that have covarage bigger then minCoverage will be counted",required = false)
    int minCoverage = 10;

    @Argument(fullName = "precentageOfSamples", shortName = "percentage", doc = "only sites where at list percentageOfSamples of the samples have good coverage, will be emited", required = false)
    double percentageOfSamples = 0.9;

    @Override
    public GenomeLoc map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());
        if ( VCs.size() == 0 )
            return null;
        if(VCs.size() != 1)
            throw new RuntimeException("there are more then one vc: "+VCs.size());

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
     * Tell the user the number of sites processed and how many passed. Close out the new intervals file.
     *
     * @param result  pair of *the number of sites seen and number of sites passed the filter.
     */
    public void onTraversalDone(Integer result) {
        logger.info(result+" sites that have "+(percentageOfSamples*100)+"% of the samples with at list "+minCoverage+" coverage.\n");
    }



}

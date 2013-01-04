package org.broadinstitute.sting.gatk.walkers.coverage;

/**
 * Created with IntelliJ IDEA.
 * User: ami
 * Date: 1/3/13
 * Time: 1:51 PM
 * To change this template use File | Settings | File Templates.
 */

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.ArgumentCollection;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContextUtils;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocSortedSet;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypesContext;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Iterator;
import java.util.Map;

@By(DataSource.REFERENCE_ORDERED_DATA)
public class CoveredByNSamplesSites extends RodWalker<Integer,Integer> implements TreeReducible<Integer> {
    @Output
    PrintStream out;

    @Output(fullName = "OutputIntervals", shortName = "intervals", doc = "Name of file for output intervals", required = true)
    File intervalFile;

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();


    @Argument(fullName = "minCoverage", shortName = "minCov",doc = "only samples that have coverage higher then the min are taking into account", required =  false)
    int minCoverage = 10;

    @Argument(fullName = "samplePercentage", shortName = "percentage",doc = "only sites that have more than samplePercentage samples with minCoverage", required = false)
    double samplePercentage = 0.9;

    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context){
        if(tracker == null)
            return 0;

        Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());
        if (vcs.size() == 0)
            return 0;

        for (VariantContext vc : vcs){
            int countCoveredSamples = 0;
            final GenotypesContext genotypes =  vc.getGenotypes();
            final int numOfGenotype = genotypes.size();
            for (Genotype g : genotypes){
                if (g.getDP() > minCoverage)
                    countCoveredSamples++;
            }
            if (countCoveredSamples/numOfGenotype >= samplePercentage)
                out.print(ref.getLocus());

        }

        return 1;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        logger.info("Processed " + result + " loci.\n");
    }
}

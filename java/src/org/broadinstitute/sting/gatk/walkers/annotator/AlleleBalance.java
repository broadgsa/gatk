package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.List;
import java.util.Map;


public class AlleleBalance extends StandardVariantAnnotation {

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {

        if ( !variation.isBiallelic() || !(variation instanceof VariantBackedByGenotype) )
            return null;
        final List<Genotype> genotypes = ((VariantBackedByGenotype)variation).getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        double ratio = 0.0;
        double totalWeights = 0.0;
        for ( Genotype genotype : genotypes ) {
            // we care only about het calls
            if ( !genotype.isHet() )
                continue;

            if ( !(genotype instanceof SampleBacked) )
                continue;

            String sample = ((SampleBacked)genotype).getSampleName();
            StratifiedAlignmentContext context = stratifiedContexts.get(sample);
            if ( context == null )
                continue;

            final String genotypeStr = genotype.getBases().toUpperCase();
            if ( genotypeStr.length() != 2 )
                return null;

            final String bases = new String(context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup().getBases()).toUpperCase();
            if ( bases.length() == 0 )
                return null;

            char a = genotypeStr.charAt(0);
            char b = genotypeStr.charAt(1);
            int aCount = Utils.countOccurrences(a, bases);
            int bCount = Utils.countOccurrences(b, bases);

            int refCount = (a == ref.getBase() ? aCount : bCount);
            int altCount = (a == ref.getBase() ? bCount : aCount);

            // sanity check
            if ( refCount + altCount == 0 )
                continue;

            // weight the allele balance by genotype quality so that e.g. mis-called homs don't affect the ratio too much
            ratio += genotype.getNegLog10PError() * ((double)refCount / (double)(refCount + altCount));
            totalWeights += genotype.getNegLog10PError();
        }

        // make sure we had a het genotype
        if ( MathUtils.compareDoubles(totalWeights, 0.0) == 0 )
            return null;

        return String.format("%.2f", (ratio / totalWeights));
    }

    public String getKeyName() { return "AB"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine("AB", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Allele Balance for hets (ref/(ref+alt))"); }
}

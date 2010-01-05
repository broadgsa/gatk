package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.SampleBacked;
import org.broadinstitute.sting.utils.genotype.vcf.VCFInfoHeaderLine;

import java.util.Map;
import java.util.List;
import java.util.ArrayList;


public class QualByDepth extends StandardVariantAnnotation {

    public String annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, Variation variation) {
        if ( !(variation instanceof VariantBackedByGenotype) )
            return null;
        final List<Genotype> genotypes = ((VariantBackedByGenotype)variation).getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        //double QbyD = genotypeQualByDepth(ref.getBase(), genotypes, stratifiedContexts);
        double QbyD = 10.0 * variation.getNegLog10PError() / (double)variationQualByDepth(ref.getBase(), genotypes, stratifiedContexts);

        return String.format("%.2f", QbyD);
    }

    public String getKeyName() { return "QD"; }

    public VCFInfoHeaderLine getDescription() { return new VCFInfoHeaderLine(getKeyName(), 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Genotype Quality by Depth"); }

    private int variationQualByDepth(char ref, final List<Genotype> genotypes, Map<String, StratifiedAlignmentContext> stratifiedContexts) {
        int depth = 0;
        for ( Genotype genotype : genotypes ) {
            if ( !(genotype instanceof SampleBacked) )
                continue;

            // we care only about variant calls
            if ( !genotype.isVariant(ref) )
                continue;

            String sample = ((SampleBacked)genotype).getSampleName();
            StratifiedAlignmentContext context = stratifiedContexts.get(sample);
            if ( context != null )
                depth += context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size();
        }

        return depth;
    }

    private double genotypeQualByDepth(char ref, final List<Genotype> genotypes, Map<String, StratifiedAlignmentContext> stratifiedContexts) {
        ArrayList<Double> qualsByDepth = new ArrayList<Double>();
        for ( Genotype genotype : genotypes ) {
            if ( !(genotype instanceof SampleBacked) )
                continue;

            // we care only about variant calls
            if ( !genotype.isVariant(ref) )
                continue;

            String sample = ((SampleBacked)genotype).getSampleName();
            StratifiedAlignmentContext context = stratifiedContexts.get(sample);
            if ( context == null )
                continue;

            qualsByDepth.add(genotype.getNegLog10PError() / context.getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).size());
        }

        double mean = 0.0;
        for ( Double value : qualsByDepth )
            mean += value;
        mean /= qualsByDepth.size();

        return mean;
    }
}
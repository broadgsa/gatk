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

package org.broadinstitute.gatk.tools.walkers.annotator;


import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.gatk.engine.contexts.AlignmentContext;
import org.broadinstitute.gatk.engine.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.gatk.utils.genotyper.MostLikelyAllele;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.pileup.PileupElement;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


/**
 * Allele balance per sample
 *
 * <p>The allele balance is the fraction of ref bases over ref + alt bases.</p>
 *
 * <h3>Caveats</h3>
 * <p>Note that this annotation will only work properly for biallelic samples that are called as heterozygous.</p>
 * <h4>This is an experimental annotation. As such, it is unsupported; we do not make any guarantees that it will work properly, and you use it at your own risk.</h4>
 */
public class AlleleBalanceBySample extends GenotypeAnnotation implements ExperimentalAnnotation {

    public void annotate(final RefMetaDataTracker tracker,
                         final AnnotatorCompatible walker,
                         final ReferenceContext ref,
                         final AlignmentContext stratifiedContext,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final PerReadAlleleLikelihoodMap alleleLikelihoodMap){


        // We need a heterozygous genotype and either a context or alleleLikelihoodMap
        if ( g == null || !g.isCalled() || !g.isHet() || ( stratifiedContext == null && alleleLikelihoodMap == null) )
            return;

        // Test for existence of <NON_REF> allele, and manually check isSNP() 
        // and isBiallelic() while ignoring the <NON_REF> allele
        boolean biallelicSNP = vc.isSNP() && vc.isBiallelic();

        if(vc.hasAllele(GVCF_NONREF)){
            // If we have the GVCF <NON_REF> allele, then the SNP is biallelic
            // iff there are 3 alleles and both the reference and first alt
            // allele are length 1.
            biallelicSNP = vc.getAlleles().size() == 3 &&
                    vc.getReference().length() == 1 &&
                    vc.getAlternateAllele(0).length() == 1;
        }

        if ( !biallelicSNP )
            return;

        Double ratio;
        if (alleleLikelihoodMap != null && !alleleLikelihoodMap.isEmpty())
            ratio = annotateWithLikelihoods(alleleLikelihoodMap, vc);
        else if ( stratifiedContext != null )
            ratio = annotateWithPileup(stratifiedContext, vc);
        else
            return;

        if (ratio == null)
            return;

        gb.attribute(getKeyNames().get(0), Double.valueOf(String.format("%.2f", ratio)));
    }

    private static final Allele GVCF_NONREF = Allele.create("<NON_REF>", false);

    private Double annotateWithPileup(final AlignmentContext stratifiedContext, final VariantContext vc) {

        final HashMap<Byte, Integer> alleleCounts = new HashMap<>();
        for ( final Allele allele : vc.getAlleles() )
            alleleCounts.put(allele.getBases()[0], 0);

        final ReadBackedPileup pileup = stratifiedContext.getBasePileup();
        for ( final PileupElement p : pileup ) {
            if ( alleleCounts.containsKey(p.getBase()) )
                alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+1);
        }

        // we need to add counts in the correct order
        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference().getBases()[0]);
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get(vc.getAlternateAllele(i).getBases()[0]);

        // sanity check
        if(counts[0] + counts[1] == 0)
            return null;

        return ((double) counts[0] / (double)(counts[0] + counts[1]));
    }

    private Double annotateWithLikelihoods(final PerReadAlleleLikelihoodMap perReadAlleleLikelihoodMap, final VariantContext vc) {
        final Set<Allele> alleles = new HashSet<>(vc.getAlleles());

        // make sure that there's a meaningful relationship between the alleles in the perReadAlleleLikelihoodMap and our VariantContext
        if ( ! perReadAlleleLikelihoodMap.getAllelesSet().containsAll(alleles) )
            throw new IllegalStateException("VC alleles " + alleles + " not a strict subset of per read allele map alleles " + perReadAlleleLikelihoodMap.getAllelesSet());

        final HashMap<Allele, Integer> alleleCounts = new HashMap<>();
        for ( final Allele allele : vc.getAlleles() ) { alleleCounts.put(allele, 0); }

        for ( final Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : perReadAlleleLikelihoodMap.getLikelihoodReadMap().entrySet()) {
            final MostLikelyAllele a = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue(), alleles);
            if (! a.isInformative() ) continue; // read is non-informative
            final int prevCount = alleleCounts.get(a.getMostLikelyAllele());
            alleleCounts.put(a.getMostLikelyAllele(), prevCount + 1);
        }

        final int[] counts = new int[alleleCounts.size()];
        counts[0] = alleleCounts.get(vc.getReference());
        for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
            counts[i+1] = alleleCounts.get( vc.getAlternateAllele(i) );

        // sanity check
        if(counts[0] + counts[1] == 0)
            return null;

        return ((double) counts[0] / (double)(counts[0] + counts[1]));

    }

    public List<String> getKeyNames() { return Arrays.asList("AB"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Allele balance for each het genotype")); }
}
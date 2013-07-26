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

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ExperimentalAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.variant.vcf.VCFFormatHeaderLine;
import org.broadinstitute.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.GenotypeBuilder;
import org.broadinstitute.variant.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;


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
        if ( stratifiedContext == null )
            return;

        Double ratio = annotateSNP(stratifiedContext, vc, g);
        if (ratio == null)
            return;

        gb.attribute(getKeyNames().get(0), Double.valueOf(String.format("%.2f", ratio.doubleValue())));
    }

    private Double annotateSNP(AlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
        double ratio = -1;

        if ( !vc.isSNP() )
            return null;

        if ( !vc.isBiallelic() )
            return null;

        if ( g == null || !g.isCalled() )
            return null;

        if (!g.isHet())
            return null;

        Collection<Allele> altAlleles = vc.getAlternateAlleles();
        if ( altAlleles.size() == 0 )
            return null;

        final String bases = new String(stratifiedContext.getBasePileup().getBases());
        if ( bases.length() == 0 )
            return null;
        char refChr = vc.getReference().toString().charAt(0);
        char altChr = vc.getAlternateAllele(0).toString().charAt(0);

        int refCount = MathUtils.countOccurrences(refChr, bases);
        int altCount = MathUtils.countOccurrences(altChr, bases);

        // sanity check
        if ( refCount + altCount == 0 )
            return null;

        ratio = ((double)refCount / (double)(refCount + altCount));
        return ratio;
    }

    public List<String> getKeyNames() { return Arrays.asList("AB"); }

    public List<VCFFormatHeaderLine> getDescriptions() { return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0), 1, VCFHeaderLineType.Float, "Allele balance for each het genotype")); }
}
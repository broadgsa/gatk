/*
 * Copyright (c) 2010 The Broad Institute
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
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.GenotypesContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * The allele balance (fraction of ref bases over ref + alt bases) across all bialleleic het-called samples
 */
public class AlleleBalance extends InfoFieldAnnotation {


    char[] BASES = {'A','C','G','T'};
    public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref, Map<String, AlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( stratifiedContexts.size() == 0 )
            return null;

        if ( !vc.isBiallelic() )
            return null;
        final GenotypesContext genotypes = vc.getGenotypes();
        if ( !vc.hasGenotypes() )
            return null;

        double ratioHom = 0.0;
        double ratioHet = 0.0;
        double weightHom = 0.0;
        double weightHet = 0.0;
        double overallNonDiploid = 0.0;
        for ( Genotype genotype : genotypes ) {
            // we care only about het calls

            AlignmentContext context = stratifiedContexts.get(genotype.getSampleName());
            if ( context == null )
                continue;

            final ReadBackedPileup pileup = context.getBasePileup();
            if ( vc.isSNP() ) {
                final String bases = new String(pileup.getBases());
                if ( bases.length() == 0 )
                    return null;

                double pTrue = 1.0 - Math.pow(10.0,genotype.getLog10PError());
                if ( genotype.isHet() ) {
                    final char refChr = vc.getReference().toString().charAt(0);
                    final char altChr = vc.getAlternateAllele(0).toString().charAt(0);

                    final int refCount = MathUtils.countOccurrences(refChr, bases);
                    final int altCount = MathUtils.countOccurrences(altChr, bases);
                    final int otherCount = bases.length()-refCount-altCount;

                    // sanity check
                    if ( refCount + altCount == 0 )
                        continue;

                    // weight the allele balance by genotype quality so that e.g. mis-called homs don't affect the ratio too much
                    ratioHet += pTrue * ((double)refCount / (double)(refCount + altCount));
                    weightHet += pTrue;
                    overallNonDiploid += ( (double) otherCount )/(bases.length()*genotypes.size());
                } else if ( genotype.isHom() ) {
                    char alleleChr;
                    if ( genotype.isHomRef() ) {
                        alleleChr = vc.getReference().toString().charAt(0);
                    } else {
                        alleleChr = vc.getAlternateAllele(0).toString().charAt(0);
                    }
                    final int alleleCount = MathUtils.countOccurrences(alleleChr,bases);
                    int bestOtherCount = 0;
                    for ( char b : BASES ) {
                        if ( b == alleleChr )
                            continue;
                        int count = MathUtils.countOccurrences(b,bases);
                        if ( count > bestOtherCount )
                            bestOtherCount = count;
                    }
                    final int otherCount = bases.length() - alleleCount;
                    ratioHom += pTrue*( (double) alleleCount)/(alleleCount+bestOtherCount);
                    weightHom += pTrue;
                    overallNonDiploid += ((double ) otherCount)/(bases.length()*genotypes.size());
                }
                // Allele Balance for indels was not being computed correctly (since there was no allele matching).  Instead of
                // prolonging the life of imperfect code, I've decided to delete it.  If someone else wants to try again from
                // scratch, be my guest - but make sure it's done correctly!  [EB]
            }
        }

        // make sure we had a het genotype

        Map<String, Object> map = new HashMap<String, Object>();
        if ( weightHet > 0.0 ) {
            map.put("ABHet",ratioHet/weightHet);
        }

        if ( weightHom > 0.0 ) {
            map.put("ABHom",ratioHom/weightHom);
        }

        if ( overallNonDiploid > 0.0 ) {
            map.put("OND",overallNonDiploid);
        }
        return map;
    }


    public List<String> getKeyNames() { return Arrays.asList("ABHet","ABHom","OND"); }

    public List<VCFInfoHeaderLine> getDescriptions() { return Arrays.asList(new VCFInfoHeaderLine("ABHet", 1, VCFHeaderLineType.Float, "Allele Balance for hets (ref/(ref+alt))"),
            new VCFInfoHeaderLine("ABHom", 1, VCFHeaderLineType.Float, "Allele Balance for homs (A/(A+O))"),
            new VCFInfoHeaderLine("OND", 1, VCFHeaderLineType.Float, "Overall non-diploid ratio (alleles/(alleles+non-alleles))")); }
}

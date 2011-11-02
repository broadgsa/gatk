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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatibleWalker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.GenotypeAnnotation;
import org.broadinstitute.sting.utils.codecs.vcf.VCFFormatHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.pileup.ExtendedEventPileupElement;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedExtendedEventPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Unsupported
 */
@Hidden
public class ReadDepthAndAllelicFractionBySample extends GenotypeAnnotation {

        private static String REF_ALLELE = "REF";

        private static String DEL = "DEL"; // constant, for speed: no need to create a key string for deletion allele every time

        public Map<String, Object> annotate(RefMetaDataTracker tracker, AnnotatorCompatibleWalker walker, ReferenceContext ref,
                                            AlignmentContext stratifiedContext, VariantContext vc, Genotype g) {
            if ( g == null || !g.isCalled() )
                return null;

            if ( vc.isSNP() )
                return annotateSNP(stratifiedContext, vc);
            if ( vc.isIndel() )
                return annotateIndel(stratifiedContext, vc);

            return null;
        }

        private Map<String,Object> annotateSNP(AlignmentContext stratifiedContext, VariantContext vc) {

            if ( ! stratifiedContext.hasBasePileup() ) return null;

            HashMap<Byte, Integer> alleleCounts = new HashMap<Byte, Integer>();
            for ( Allele allele : vc.getAlternateAlleles() )
                alleleCounts.put(allele.getBases()[0], 0);

            ReadBackedPileup pileup = stratifiedContext.getBasePileup();
            int totalDepth = pileup.getNumberOfElements();

            Map<String, Object> map = new HashMap<String, Object>();
            map.put(getKeyNames().get(0), totalDepth); // put total depth in right away

            if ( totalDepth == 0 ) return map; // done, can not compute FA at 0 coverage!!

            int mq0 = 0; // number of "ref" reads that are acually mq0
            for ( PileupElement p : pileup ) {
                if ( p.getMappingQual() == 0 ) {
                    mq0++;
                    continue;
                }
                if ( alleleCounts.containsKey(p.getBase()) ) // non-mq0 read and it's an alt
                    alleleCounts.put(p.getBase(), alleleCounts.get(p.getBase())+1);
            }

            if ( mq0 == totalDepth ) return map; // if all reads are mq0, there is nothing left to do

            // we need to add counts in the correct order
            String[] fracs = new String[alleleCounts.size()];
            for (int i = 0; i < vc.getAlternateAlleles().size(); i++) {
                fracs[i] = String.format("%.3f", ((float)alleleCounts.get(vc.getAlternateAllele(i).getBases()[0]))/(totalDepth-mq0));
            }

            map.put(getKeyNames().get(1), fracs);
            return map;
        }

        private Map<String,Object> annotateIndel(AlignmentContext
            stratifiedContext, VariantContext
            vc) {

            if ( ! stratifiedContext.hasExtendedEventPileup() ) {
                return null;
            }

            ReadBackedExtendedEventPileup pileup = stratifiedContext.getExtendedEventPileup();
            if ( pileup == null )
                return null;
            int totalDepth = pileup.getNumberOfElements();

            Map<String, Object> map = new HashMap<String, Object>();
            map.put(getKeyNames().get(0), totalDepth); // put total depth in right away

            if ( totalDepth == 0 ) return map;
            int mq0 = 0; // number of "ref" reads that are acually mq0

            HashMap<String, Integer> alleleCounts = new HashMap<String, Integer>();
            Allele refAllele = vc.getReference();

            for ( Allele allele : vc.getAlternateAlleles() ) {

                if ( allele.isNoCall() ) {
                    continue; // this does not look so good, should we die???
                }

                alleleCounts.put(getAlleleRepresentation(allele), 0);
            }

            for ( ExtendedEventPileupElement e : pileup.toExtendedIterable() ) {

                if ( e.getMappingQual() == 0 ) {
                    mq0++;
                    continue;
                }

                if ( e.isInsertion() ) {

                    final String b =  e.getEventBases();
                    if ( alleleCounts.containsKey(b) ) {
                        alleleCounts.put(b, alleleCounts.get(b)+1);
                    }

                } else {
                    if ( e.isDeletion() ) {
                        if ( e.getEventLength() == refAllele.length() ) {
                            // this is indeed the deletion allele recorded in VC
                            final String b = DEL;
                            if ( alleleCounts.containsKey(b) ) {
                                alleleCounts.put(b, alleleCounts.get(b)+1);
                            }
                        }
//                    else {
//                        System.out.print("   deletion of WRONG length found");
//                    }
                    }
                }
            }

            if ( mq0 == totalDepth ) return map;

            String[] fracs = new String[alleleCounts.size()];
            for (int i = 0; i < vc.getAlternateAlleles().size(); i++)
                fracs[i] = String.format("%.3f",
                        ((float)alleleCounts.get(getAlleleRepresentation(vc.getAlternateAllele(i))))/(totalDepth-mq0));

            map.put(getKeyNames().get(1), fracs);

            //map.put(getKeyNames().get(0), counts);
            return map;
        }

        private String getAlleleRepresentation(Allele allele) {
            if ( allele.isNull() ) { // deletion wrt the ref
                 return DEL;
            } else { // insertion, pass actual bases
                return allele.getBaseString();
            }

        }

     //   public String getIndelBases()
        public List<String> getKeyNames() { return Arrays.asList("DP","FA"); }

        public List<VCFFormatHeaderLine> getDescriptions() {
            return Arrays.asList(new VCFFormatHeaderLine(getKeyNames().get(0),
                            1,
                            VCFHeaderLineType.Integer,
                            "Total read depth per sample, including MQ0"),
                            new VCFFormatHeaderLine(getKeyNames().get(1),
                            VCFHeaderLineCount.UNBOUNDED,
                            VCFHeaderLineType.Float,
                            "Fractions of reads (excluding MQ0 from both ref and alt) supporting each reported alternative allele, per sample"));
        }
}

/*
* Copyright 2012-2016 Broad Institute, Inc.
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
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.gatk.tools.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.gatk.utils.MathUtils;
import org.broadinstitute.gatk.utils.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Allele balance across all samples
 *
 * <p> This is a set of experimental annotations that attempt to estimate whether the data supporting a variant call fits allelic ratio expectations, or whether there might be some bias in the data. ABHom is the proportion of reads from homozygous samples that support the call (REF or ALT depending on whether the call is hom-ref or hom-var). ABHet is the proportion of REF reads from heterozygous samples. OND represents the overall fraction of data that diverges from the diploid hypothesis, based on the number of reads that support something other than the genotyped alleles (called "non-alleles"). Note that each sample's contribution is weighted by its genotype quality so that individual mis-calls don't affect the overall ratio too much.</p>
 * <h3>Calculations</h3>
 * <p> $$ ABHom = \frac{# REF or ALT reads from homozygous samples}{# REF + ALT reads from homozygous samples} $$ <br />
 *     $$ ABHet = \frac{# REF reads from heterozygous samples}{# REF + ALT reads from heterozygous samples} $$ <br />
 *     $$ OND = \frac{# reads from non-alleles}{# all reads} $$
 * </p>
 * <p> For ABHom, the value should be close to 1.00 because ideally, all the reads should support a single allele. For ABHet, the value should be close to 0.5, so half of the reads support the ref allele and half of the reads support the alt allele. Divergence from these expected ratios may indicate that there is some bias in favor of one allele. Note the caveats below regarding cancer and RNAseq analysis. </p>
 * <h3>Caveats</h3>
 * <ul>
 *     <li>This annotation will only work properly for biallelic SNPs in diploid organisms where all samples are either called heterozygous or homozygous.</li>
 *     <li>This annotation cannot currently be calculated for indels.</li>
 *     <li>The reasoning underlying this annotation only applies to germline variants in DNA sequencing data. In somatic/cancer analysis, divergent ratios are expected due to tumor heterogeneity and normal contamination. In RNAseq analysis, divergent ratios may indicate differential allele expression.</li>
 *     <li>As stated above, this annotation is experimental and should be interpreted with caution as we cannot guarantee that it is appropriate. Basically, use it at your own risk.</li>
 * </ul>
 * <h3>Related annotations</h3>
 * <ul>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_AlleleBalanceBySample.php">AlleleBallanceBySample</a></b> calculates allele balance for each individual sample.</li>
 *     <li><b><a href="https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_annotator_DepthPerAlleleBySample.php">DepthPerAlleleBySample</a></b> calculates depth of coverage for each allele per sample.</li>
 * </ul>
 */

public class AlleleBalance extends InfoFieldAnnotation implements ActiveRegionBasedAnnotation {

    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( !(vc.isBiallelic() && vc.hasGenotypes())) {
            return null;
        }

        double refCountInHetSamples = 0.0;
        double altCountInHetSamples = 0.0;
        double correctCountInHomSamples = 0.0;
        double incorrectCountInHomSamples = 0.0;
        double nonDiploidCount = 0.0;
        double totalReadCount = 0.0;

        for ( final Genotype genotype : vc.getGenotypes() ) {
            if ( !vc.isSNP() ) {
                continue;
            }

            final int[] alleleCounts = getCounts(genotype, stratifiedContexts, vc);

            if (alleleCounts == null) continue;

            final long totalReads = MathUtils.sum(alleleCounts);
            if ( genotype.isHet() ) {
                // weight read counts by genotype quality so that e.g. mis-called homs don't affect the ratio too much
                refCountInHetSamples += alleleCounts[0];
                altCountInHetSamples += alleleCounts[1];
                nonDiploidCount += totalReads - (alleleCounts[0] + alleleCounts[1]);
                totalReadCount += totalReads;
            } else if ( genotype.isHom() ) {
                final int alleleIndex = genotype.isHomRef() ?  0 : 1 ;
                final int alleleCount = alleleCounts[alleleIndex];
                int bestOtherCount = 0;
                for(int n = 0; n < alleleCounts.length; n++){
                    if( n != alleleIndex && alleleCounts[n] > bestOtherCount ) {
                        bestOtherCount = alleleCounts[n];
                    }
                }
                correctCountInHomSamples += alleleCount;
                incorrectCountInHomSamples += bestOtherCount;
                nonDiploidCount += totalReads - alleleCount;
                totalReadCount += totalReads;
            }
            // Allele Balance for indels was not being computed correctly (since there was no allele matching).  Instead of
            // prolonging the life of imperfect code, I've decided to delete it.  If someone else wants to try again from
            // scratch, be my guest - but make sure it's done correctly!  [EB]

        }
        final double diploidCountInHetSamples = altCountInHetSamples + refCountInHetSamples;
        final double diploidCountInHomSamples = correctCountInHomSamples + incorrectCountInHomSamples;

        Map<String, Object> map = new HashMap<>();
        if ( diploidCountInHetSamples > 0.0 ) {
            map.put(GATKVCFConstants.ALLELE_BALANCE_HET_KEY, refCountInHetSamples / diploidCountInHetSamples);
        }

        if ( diploidCountInHomSamples > 0.0 ) {
            map.put(GATKVCFConstants.ALLELE_BALANCE_HOM_KEY, correctCountInHomSamples / diploidCountInHomSamples);
        }

        if ( totalReadCount > 0.0 ) {
            map.put(GATKVCFConstants.NON_DIPLOID_RATIO_KEY, nonDiploidCount / totalReadCount);
        }
        return map;
    }

    /**
     * Get the number of reads per allele, using the following (in order of preference):
     * - genotype.getAD()
     * - reads from an AlignmentContext
     * - reads from a PerReadAlleleLikelihoodMap (Not yet implemented) 
     *
     * @param genotype The genotype of interest
     * @param stratifiedContexts A mapping of sample name to read alignments at a location
     * @param vc The Variant Context
     * @return The number of reads per allele
     */
    private int[] getCounts(final Genotype genotype,
                            final Map<String, AlignmentContext> stratifiedContexts,
                            final VariantContext vc){
        if(genotype == null)
            return null;

        if (genotype.hasAD()) {
            return genotype.getAD();
        } else {    // If  getAD() returned no information we count alleles from the pileup
            final AlignmentContext context = stratifiedContexts == null ? null : stratifiedContexts.get(genotype.getSampleName());
            if (context == null) return null;

            final byte[] bases = context.getBasePileup().getBases();
            // Should be able to replace with the following, but this annotation was not found when using -A AlleleBalance
            // return vc.getAlleles().stream().map(a -> MathUtils.countOccurrences(a.getBases()[0], bases)).mapToInt(Integer::intValue).toArray();
            final List<Allele> alleles = vc.getAlleles();
            final int[] result = new int[alleles.size()];
            // Calculate the depth for each allele, assuming that the allele is a single base
            for(int n = 0; n < alleles.size(); n++){
                result[n] = MathUtils.countOccurrences(alleles.get(n).getBases()[0], bases);
            }
            return result;

        }
    }

    @Override
    public List<String> getKeyNames() {
        return Arrays.asList(GATKVCFConstants.ALLELE_BALANCE_HET_KEY,
                             GATKVCFConstants.ALLELE_BALANCE_HOM_KEY,
                             GATKVCFConstants.NON_DIPLOID_RATIO_KEY);
    }
}

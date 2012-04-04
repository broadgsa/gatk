/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

@Analysis(description = "Evaluation summary for indels")
public class IndelSummary extends VariantEvaluator implements StandardEval {
    final protected static Logger logger = Logger.getLogger(IndelSummary.class);

    @DataPoint(description = "Number of SNPs", format = "%d")
    public int n_SNPs = 0;

    @DataPoint(description = "Number of singleton SNPs", format = "%d")
    public int n_singleton_SNPs = 0;

    @DataPoint(description = "Number of Indels", format = "%d")
    public int n_indels = 0;

    // Number of Indels Sites (counts one for any number of alleles at site)
    public int nIndelSites = 0;

    @DataPoint(description = "Number of singleton Indels", format = "%d")
    public int n_singleton_indels = 0;

    // counts 1 for each site where the number of alleles > 2
    public int nMultiIndelSites = 0;

    @DataPoint(description = "Percent of indel sites that are multi-allelic")
    public String percent_of_sites_with_more_than_2_alleles;

    @DataPoint(description = "SNP to indel ratio")
    public String SNP_to_indel_ratio;

    @DataPoint(description = "Singleton SNP to indel ratio")
    public String SNP_to_indel_ratio_for_singletons;

    @DataPoint(description = "Indel novelty rate")
    public String indel_novelty_rate;

    @DataPoint(description = "1 to 2 bp indel ratio")
    public String ratio_of_1_to_2_bp_indels;

    @DataPoint(description = "1 to 3 bp indel ratio")
    public String ratio_of_1_to_3_bp_indels;

    @DataPoint(description = "2 to 3 bp indel ratio")
    public String ratio_of_2_to_3_bp_indels;

    @DataPoint(description = "1 and 2 to 3 bp indel ratio")
    public String ratio_of_1_and_2_to_3_bp_indels;

    @DataPoint(description = "Frameshift percent")
    public String frameshift_rate_for_coding_indels;

    //
    // insertions to deletions
    //
    @DataPoint(description = "Insertion to deletion ratio")
    public String insertion_to_deletion_ratio;

    @DataPoint(description = "Insertion to deletion ratio for 1 bp events")
    public String insertion_to_deletion_ratio_for_1bp_indels;

    //
    // Frameshifts
    //
    @DataPoint(description = "Number of indels in protein-coding regions labeled as frameshift")
    public int n_coding_indels_frameshifting = 0;

    @DataPoint(description = "Number of indels in protein-coding regions not labeled as frameshift")
    public int n_coding_indels_in_frame = 0;

    //
    // Het : hom ratios
    //
    @DataPoint(description = "Het to hom ratio for SNPs")
    public String SNP_het_to_hom_ratio;

    @DataPoint(description = "Het to hom ratio for indels")
    public String indel_het_to_hom_ratio;
    
    int nSNPHets = 0, nSNPHoms = 0, nIndelHets = 0, nIndelHoms = 0;

    int nKnownIndels = 0, nInsertions = 0;
    int n1bpInsertions = 0, n1bpDeletions = 0;
    int[] countByLength = new int[]{0, 0, 0, 0}; // note that the first element isn't used


    public final static int LARGE_INDEL_SIZE_THRESHOLD = 10;

    @DataPoint(description = "Number of large (>10 bp) deletions")
    public int n_large_deletions = 0;

    @DataPoint(description = "Number of large (>10 bp) insertions")
    public int n_large_insertions = 0;

    @DataPoint(description = "Ratio of large (>10 bp) insertions to deletions")
    public String insertion_to_deletion_ratio_for_large_indels;

    @Override public int getComparisonOrder() { return 2; }

    public void update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval == null || (getWalker().ignoreAC0Sites() && eval.isMonomorphicInSamples()) )
            return;

        // update counts
        switch ( eval.getType() ) {
            case SNP:
                n_SNPs += eval.getNAlleles() - 1; // -1 for ref
                if ( variantWasSingleton(eval) ) n_singleton_SNPs++;

                // collect information about het / hom ratio
                for ( final Genotype g : eval.getGenotypes() ) {
                    if ( g.isHet() ) nSNPHets++;
                    if ( g.isHomVar() ) nSNPHoms++;
                }
                break;
            case INDEL:
                if ( eval.isComplexIndel() ) break; // don't count complex substitutions
                
                nIndelSites++;
                if ( ! eval.isBiallelic() ) nMultiIndelSites++;
                if ( variantWasSingleton(eval) ) n_singleton_indels++;

                // collect information about het / hom ratio
                for ( final Genotype g : eval.getGenotypes() ) {
                    if ( g.isHet() ) nIndelHets++;
                    if ( g.isHomVar() ) nIndelHoms++;
                }

                for ( Allele alt : eval.getAlternateAlleles() ) {
                    n_indels++; // +1 for each alt allele

                    if ( comp != null ) nKnownIndels++; // TODO -- make this test allele specific?

                    // ins : del ratios
                    final int alleleSize = alt.length() - eval.getReference().length();
                    if ( alleleSize == 0 ) throw new ReviewedStingException("Allele size not expected to be zero for indel: alt = " + alt + " ref = " + eval.getReference());
                    if ( alleleSize > 0 ) nInsertions++;
                    if ( alleleSize == 1 ) n1bpInsertions++;
                    if ( alleleSize == -1 ) n1bpDeletions++;

                    // requires snpEFF annotations
                    if ( eval.getAttributeAsString("SNPEFF_GENE_BIOTYPE", "missing").equals("protein_coding") ) {
                        final String effect = eval.getAttributeAsString("SNPEFF_EFFECT", "missing");
                        if ( effect.equals("missing") ) 
                            throw new ReviewedStingException("Saw SNPEFF_GENE_BIOTYPE but unexpected no SNPEFF_EFFECT at " + eval);
                        if ( effect.equals("FRAME_SHIFT") )
                            n_coding_indels_frameshifting++;
                        else if ( effect.startsWith("CODON") )
                            n_coding_indels_in_frame++;
                        else
                            ; // lots of protein coding effects that shouldn't be counted, such as INTRON
                    }

                    if ( alleleSize > LARGE_INDEL_SIZE_THRESHOLD )
                        n_large_insertions++;
                    else if ( alleleSize < -LARGE_INDEL_SIZE_THRESHOLD )
                        n_large_deletions++;
                    
                    // update the baby histogram
                    final int absSize = Math.abs(alleleSize);
                    if ( absSize < countByLength.length ) countByLength[absSize]++;

                }

                break;
            default:
                throw new UserException.BadInput("Unexpected variant context type: " + eval);
        }

        return;
    }

    public void finalizeEvaluation() {
        percent_of_sites_with_more_than_2_alleles = Utils.formattedRatio(nMultiIndelSites, nIndelSites);
        SNP_to_indel_ratio = Utils.formattedRatio(n_SNPs, n_indels);
        SNP_to_indel_ratio_for_singletons = Utils.formattedRatio(n_singleton_SNPs, n_singleton_indels);
        indel_novelty_rate = Utils.formattedNoveltyRate(nKnownIndels, n_indels);
        ratio_of_1_to_2_bp_indels = Utils.formattedRatio(countByLength[1], countByLength[2]);
        ratio_of_1_to_3_bp_indels = Utils.formattedRatio(countByLength[1], countByLength[3]);
        ratio_of_2_to_3_bp_indels = Utils.formattedRatio(countByLength[2], countByLength[3]);
        ratio_of_1_and_2_to_3_bp_indels = Utils.formattedRatio(countByLength[1] + countByLength[2], countByLength[3]);
        frameshift_rate_for_coding_indels = Utils.formattedPercent(n_coding_indels_frameshifting, n_coding_indels_in_frame + n_coding_indels_frameshifting);

        SNP_het_to_hom_ratio = Utils.formattedRatio(nSNPHets, nSNPHoms);
        indel_het_to_hom_ratio = Utils.formattedRatio(nIndelHets, nIndelHoms);

        insertion_to_deletion_ratio = Utils.formattedRatio(nInsertions, n_indels - nInsertions);
        insertion_to_deletion_ratio_for_1bp_indels = Utils.formattedRatio(n1bpInsertions, n1bpDeletions);
        insertion_to_deletion_ratio_for_large_indels = Utils.formattedRatio(n_large_insertions, n_large_deletions);

    }
}

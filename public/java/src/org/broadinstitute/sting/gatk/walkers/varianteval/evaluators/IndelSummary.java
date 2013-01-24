/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
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
import org.broadinstitute.variant.variantcontext.Allele;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;

@Analysis(description = "Evaluation summary for indels")
public class IndelSummary extends VariantEvaluator implements StandardEval {
    final protected static Logger logger = Logger.getLogger(IndelSummary.class);

    //
    // counts of snps and indels
    //
    @DataPoint(description = "Number of SNPs", format = "%d")
    public int n_SNPs = 0;

    @DataPoint(description = "Number of singleton SNPs", format = "%d")
    public int n_singleton_SNPs = 0;

    @DataPoint(description = "Number of indels", format = "%d")
    public int n_indels = 0;

    @DataPoint(description = "Number of singleton indels", format = "%d")
    public int n_singleton_indels = 0;

    //
    // gold standard
    //
    @DataPoint(description = "Number of Indels overlapping gold standard sites", format = "%d")
    public int n_indels_matching_gold_standard = 0;

    @DataPoint(description = "Percent of indels overlapping gold standard sites")
    public String gold_standard_matching_rate;

    //
    // multi-allelics
    //
    // Number of Indels Sites (counts one for any number of alleles at site)
    public int nIndelSites = 0;

    @DataPoint(description = "Number of sites with where the number of alleles is greater than 2")
    public int n_multiallelic_indel_sites = 0;

    @DataPoint(description = "Percent of indel sites that are multi-allelic")
    public String percent_of_sites_with_more_than_2_alleles;

    //
    // snp : indel ratios
    //
    @DataPoint(description = "SNP to indel ratio")
    public String SNP_to_indel_ratio;

    @DataPoint(description = "Singleton SNP to indel ratio")
    public String SNP_to_indel_ratio_for_singletons;

    //
    // novelty
    //
    @DataPoint(description = "Number of novel indels", format = "%d")
    public int n_novel_indels = 0;

    @DataPoint(description = "Indel novelty rate")
    public String indel_novelty_rate;

    //
    // insertions to deletions
    //
    @DataPoint(description = "Number of insertion indels")
    public int n_insertions = 0;

    @DataPoint(description = "Number of deletion indels")
    public int n_deletions = 0;

    @DataPoint(description = "Insertion to deletion ratio")
    public String insertion_to_deletion_ratio;

    @DataPoint(description = "Number of large (>10 bp) deletions")
    public int n_large_deletions = 0;

    @DataPoint(description = "Number of large (>10 bp) insertions")
    public int n_large_insertions = 0;

    @DataPoint(description = "Ratio of large (>10 bp) insertions to deletions")
    public String insertion_to_deletion_ratio_for_large_indels;

    //
    // Frameshifts
    //
    @DataPoint(description = "Number of indels in protein-coding regions labeled as frameshift")
    public int n_coding_indels_frameshifting = 0;

    @DataPoint(description = "Number of indels in protein-coding regions not labeled as frameshift")
    public int n_coding_indels_in_frame = 0;

    @DataPoint(description = "Frameshift percent")
    public String frameshift_rate_for_coding_indels;

    //
    // Het : hom ratios
    //
    @DataPoint(description = "Het to hom ratio for SNPs")
    public String SNP_het_to_hom_ratio;

    @DataPoint(description = "Het to hom ratio for indels")
    public String indel_het_to_hom_ratio;
    
    int nSNPHets = 0, nSNPHoms = 0, nIndelHets = 0, nIndelHoms = 0;

    int[] insertionCountByLength = new int[]{0, 0, 0, 0}; // note that the first element isn't used
    int[] deletionCountByLength = new int[]{0, 0, 0, 0}; // note that the first element isn't used

    // - Since 1 & 2 bp insertions and 1 & 2 bp deletions are equally likely to cause a
    // downstream frameshift, if we make the simplifying assumptions that 3 bp ins
    // and 3bp del (adding/subtracting 1 AA in general) are roughly comparably
    // selected against, we should see a consistent 1+2 : 3 bp ratio for insertions
    // as for deletions, and certainly would expect consistency between in/dels that
    // multiple methods find and in/dels that are unique to one method  (since deletions
    // are more common and the artifacts differ, it is probably worth looking at the totals,
    // overlaps and ratios for insertions and deletions separately in the methods
    // comparison and in this case don't even need to make the simplifying in = del functional assumption

    @DataPoint(description = "ratio of 1 and 2 bp insertions to 3 bp insertions")
    public String ratio_of_1_and_2_to_3_bp_insertions;

    @DataPoint(description = "ratio of 1 and 2 bp deletions to 3 bp deletions")
    public String ratio_of_1_and_2_to_3_bp_deletions;

    public final static int LARGE_INDEL_SIZE_THRESHOLD = 10;

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
                final VariantContext gold = getWalker().goldStandard == null ? null : tracker.getFirstValue(getWalker().goldStandard);

                nIndelSites++;
                if ( ! eval.isBiallelic() ) n_multiallelic_indel_sites++;

                // collect information about het / hom ratio
                for ( final Genotype g : eval.getGenotypes() ) {
                    if ( g.isHet() ) nIndelHets++;
                    if ( g.isHomVar() ) nIndelHoms++;
                }

                for ( Allele alt : eval.getAlternateAlleles() ) {
                    n_indels++; // +1 for each alt allele
                    if ( variantWasSingleton(eval) ) n_singleton_indels++;
                    if ( comp == null ) n_novel_indels++; // TODO -- make this test allele specific?
                    if ( gold != null ) n_indels_matching_gold_standard++;

                    // ins : del ratios
                    final int alleleSize = alt.length() - eval.getReference().length();
                    if ( alleleSize == 0 ) throw new ReviewedStingException("Allele size not expected to be zero for indel: alt = " + alt + " ref = " + eval.getReference());
                    if ( alleleSize > 0 ) n_insertions++;
                    if ( alleleSize < 0 ) n_deletions++;

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
                    final int[] countByLength = alleleSize < 0 ? deletionCountByLength : insertionCountByLength;
                    final int absSize = Math.abs(alleleSize);
                    if ( absSize < countByLength.length ) countByLength[absSize]++;

                }

                break;
            default:
                // TODO - MIXED, SYMBOLIC, and MNP records are skipped over
                //throw new UserException.BadInput("Unexpected variant context type: " + eval);
                break;
        }

        return;
    }

    public void finalizeEvaluation() {
        percent_of_sites_with_more_than_2_alleles = Utils.formattedPercent(n_multiallelic_indel_sites, nIndelSites);
        SNP_to_indel_ratio = Utils.formattedRatio(n_SNPs, n_indels);
        SNP_to_indel_ratio_for_singletons = Utils.formattedRatio(n_singleton_SNPs, n_singleton_indels);

        gold_standard_matching_rate = Utils.formattedPercent(n_indels_matching_gold_standard, n_indels);
        indel_novelty_rate = Utils.formattedNoveltyRate(n_indels - n_novel_indels, n_indels);
        frameshift_rate_for_coding_indels = Utils.formattedPercent(n_coding_indels_frameshifting, n_coding_indels_in_frame + n_coding_indels_frameshifting);

        ratio_of_1_and_2_to_3_bp_deletions = Utils.formattedRatio(deletionCountByLength[1] + deletionCountByLength[2], deletionCountByLength[3]);
        ratio_of_1_and_2_to_3_bp_insertions = Utils.formattedRatio(insertionCountByLength[1] + insertionCountByLength[2], insertionCountByLength[3]);

        SNP_het_to_hom_ratio = Utils.formattedRatio(nSNPHets, nSNPHoms);
        indel_het_to_hom_ratio = Utils.formattedRatio(nIndelHets, nIndelHoms);

        insertion_to_deletion_ratio = Utils.formattedRatio(n_insertions, n_deletions);
        insertion_to_deletion_ratio_for_large_indels = Utils.formattedRatio(n_large_insertions, n_large_deletions);

    }
}

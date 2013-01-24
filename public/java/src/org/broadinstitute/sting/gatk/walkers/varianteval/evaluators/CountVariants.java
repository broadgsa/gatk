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

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.variant.variantcontext.Genotype;
import org.broadinstitute.variant.variantcontext.VariantContext;

@Analysis(description = "Counts different classes of variants in the sample")
public class CountVariants extends VariantEvaluator implements StandardEval {
    // the following fields are in output order:

    // basic counts on various rates found
    @DataPoint(description = "Number of processed loci", format = "%d")
    public long nProcessedLoci = 0;
    @DataPoint(description = "Number of called loci", format = "%d")
    public long nCalledLoci = 0;
    @DataPoint(description = "Number of reference loci", format = "%d")
    public long nRefLoci = 0;
    @DataPoint(description = "Number of variant loci", format = "%d")
    public long nVariantLoci = 0;

    // the following two calculations get set in the finalizeEvaluation
    @DataPoint(description = "Variants per loci rate", format = "%.8f")
    public double variantRate = 0;
    @DataPoint(description = "Number of variants per base", format = "%.8f")
    public double variantRatePerBp = 0;

    @DataPoint(description = "Number of snp loci", format = "%d")
    public long nSNPs = 0;
    @DataPoint(description = "Number of mnp loci", format = "%d")
    public long nMNPs = 0;
    @DataPoint(description = "Number of insertions", format = "%d")
    public long nInsertions = 0;
    @DataPoint(description = "Number of deletions", format = "%d")
    public long nDeletions = 0;
    @DataPoint(description = "Number of complex indels", format = "%d")
    public long nComplex = 0;
    @DataPoint(description = "Number of symbolic events", format = "%d")
    public long nSymbolic = 0;

    @DataPoint(description = "Number of mixed loci (loci that can't be classified as a SNP, Indel or MNP)", format = "%d")
    public long nMixed = 0;

    @DataPoint(description = "Number of no calls loci", format = "%d")
    public long nNoCalls = 0;
    @DataPoint(description = "Number of het loci", format = "%d")
    public long nHets = 0;
    @DataPoint(description = "Number of hom ref loci", format = "%d")
    public long nHomRef = 0;
    @DataPoint(description = "Number of hom var loci", format = "%d")
    public long nHomVar = 0;
    @DataPoint(description = "Number of singletons", format = "%d")
    public long nSingletons = 0;
    @DataPoint(description = "Number of derived homozygotes", format = "%d")
    public long nHomDerived = 0;

    // calculations that get set in the finalizeEvaluation method
    @DataPoint(description = "heterozygosity per locus rate", format = "%.2e")
    public double heterozygosity = 0;
    @DataPoint(description = "heterozygosity per base pair", format = "%.2f")
    public double heterozygosityPerBp = 0;
    @DataPoint(description = "heterozygosity to homozygosity ratio", format = "%.2f")
    public double hetHomRatio = 0;
    @DataPoint(description = "indel rate (insertion count + deletion count)", format = "%.2e")
    public double indelRate = 0;
    @DataPoint(description = "indel rate per base pair", format = "%.2f")
    public double indelRatePerBp = 0;
    @DataPoint(description = "insertion  to deletion ratio", format = "%.2f")
    public double insertionDeletionRatio = 0;
    
    private double perLocusRate(long n) {
        return rate(n, nProcessedLoci);
    }

    private long perLocusRInverseRate(long n) {
        return inverseRate(n, nProcessedLoci);
    }


    public int getComparisonOrder() {
        return 1;   // we only need to see each eval track
    }

    public void update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        nCalledLoci++;

        // Note from Eric:
        // This is really not correct.  What we really want here is a polymorphic vs. monomorphic count (i.e. on the Genotypes).
        // So in order to maintain consistency with the previous implementation (and the intention of the original author), I've
        // added in a proxy check for monomorphic status here.
        // Protect against case when vc only as no-calls too - can happen if we strafity by sample and sample as a single no-call.
       if ( getWalker().ignoreAC0Sites() && vc1.isMonomorphicInSamples() ) {
            nRefLoci++;
        } else {
             switch (vc1.getType()) {
                case NO_VARIATION:
                    // shouldn't get here
                    break;
                case SNP:
                    nVariantLoci++;
                    nSNPs++;
                    if (variantWasSingleton(vc1)) nSingletons++;
                    break;
                case MNP:
                    nVariantLoci++;
                    nMNPs++;
                    if (variantWasSingleton(vc1)) nSingletons++;
                    break;
                case INDEL:
                    nVariantLoci++;
                    if (vc1.isSimpleInsertion())
                        nInsertions++;
                    else if (vc1.isSimpleDeletion())
                        nDeletions++;
                    else
                        nComplex++;
                    break;
                case MIXED:
                    nVariantLoci++;
                    nMixed++;
                    break;
                case SYMBOLIC:
                    nSymbolic++;
                    break;
                default:
                    throw new ReviewedStingException("Unexpected VariantContext type " + vc1.getType());
            }
        }

        // these operations are ordered to ensure that we don't get the base string of the ref unless we need it
        final String aaStr = vc1.hasAttribute("ANCESTRALALLELE") ? vc1.getAttributeAsString("ANCESTRALALLELE", null).toUpperCase() : null;
        final String refStr = aaStr != null ? vc1.getReference().getBaseString().toUpperCase() : null;

        // ref  aa  alt  class
        // A    C   A    der homozygote
        // A    C   C    anc homozygote

        // A    A   A    ref homozygote
        // A    A   C
        // A    C   A
        // A    C   C

        for (final Genotype g : vc1.getGenotypes()) {
            final String altStr = vc1.getAlternateAlleles().size() > 0 ? vc1.getAlternateAllele(0).getBaseString().toUpperCase() : null;

            switch (g.getType()) {
                case NO_CALL:
                    nNoCalls++;
                    break;
                case HOM_REF:
                    nHomRef++;

                    if ( aaStr != null && altStr != null && !refStr.equalsIgnoreCase(aaStr) ) {
                        nHomDerived++;
                    }

                    break;
                case HET:
                    nHets++;
                    break;
                case HOM_VAR:
                    nHomVar++;

                    if ( aaStr != null && altStr != null && !altStr.equalsIgnoreCase(aaStr) ) {
                        nHomDerived++;
                    }

                    break;
                case MIXED:
                    break;
                default:
                    throw new ReviewedStingException("BUG: Unexpected genotype type: " + g);
            }
        }
    }

    public void finalizeEvaluation() {
        nProcessedLoci = getWalker().getnProcessedLoci();
        variantRate = perLocusRate(nVariantLoci);
        variantRatePerBp = perLocusRInverseRate(nVariantLoci);
        heterozygosity = perLocusRate(nHets);
        heterozygosityPerBp = perLocusRInverseRate(nHets);
        hetHomRatio = ratio(nHets, nHomVar);
        indelRate = perLocusRate(nDeletions + nInsertions + nComplex);
        indelRatePerBp = perLocusRInverseRate(nDeletions + nInsertions + nComplex);
        insertionDeletionRatio = ratio(nInsertions, nDeletions);
    }
}
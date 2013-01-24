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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import com.google.java.contract.Ensures;
import com.google.java.contract.Invariant;
import com.google.java.contract.Requires;
import org.broad.tribble.util.ParsingUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.vcf.VCFHeader;

import java.util.*;

/**
 * A class for tabulating and evaluating a callset-by-callset genotype concordance table
 * */
public class ConcordanceMetrics {

    private Map<String,GenotypeConcordanceTable> perSampleGenotypeConcordance;
    private GenotypeConcordanceTable overallGenotypeConcordance;
    private SiteConcordanceTable overallSiteConcordance;

    public ConcordanceMetrics(VCFHeader evaluate, VCFHeader truth) {
        HashSet<String> overlappingSamples = new HashSet<String>(evaluate.getGenotypeSamples());
        overlappingSamples.retainAll(truth.getGenotypeSamples());
        perSampleGenotypeConcordance = new HashMap<String, GenotypeConcordanceTable>(overlappingSamples.size());
        for ( String sample : overlappingSamples ) {
            perSampleGenotypeConcordance.put(sample,new GenotypeConcordanceTable());
        }
        overallGenotypeConcordance = new GenotypeConcordanceTable();
        overallSiteConcordance = new SiteConcordanceTable();
    }

    public GenotypeConcordanceTable getOverallGenotypeConcordance() {
        return overallGenotypeConcordance;
    }

    public SiteConcordanceTable getOverallSiteConcordance() {
        return overallSiteConcordance;
    }

    public GenotypeConcordanceTable getGenotypeConcordance(String sample) {
        GenotypeConcordanceTable table = perSampleGenotypeConcordance.get(sample);
        if ( table == null )
            throw new ReviewedStingException("Attempted to request the concordance table for sample "+sample+" on which it was not calculated");
        return table;
    }

    public Map<String,GenotypeConcordanceTable> getPerSampleGenotypeConcordance() {
        return Collections.unmodifiableMap(perSampleGenotypeConcordance);
    }

    public Map<String,Double> getPerSampleNRD() {
        Map<String,Double> nrd = new HashMap<String,Double>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            nrd.put(sampleTable.getKey(),calculateNRD(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(nrd);
    }

    public Double getOverallNRD() {
        return calculateNRD(overallGenotypeConcordance);
    }

    public Map<String,Double> getPerSampleNRS() {
        Map<String,Double> nrs = new HashMap<String,Double>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            nrs.put(sampleTable.getKey(),calculateNRS(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(nrs);
    }

    public Double getOverallNRS() {
        return calculateNRS(overallGenotypeConcordance);
    }

    @Requires({"eval != null","truth != null"})
    public void update(VariantContext eval, VariantContext truth) {
        overallSiteConcordance.update(eval,truth);
        Set<String> alleleTruth = new HashSet<String>(8);
        alleleTruth.add(truth.getReference().getBaseString());
        for ( Allele a : truth.getAlternateAlleles() ) {
            alleleTruth.add(a.getBaseString());
        }
        for ( String sample : perSampleGenotypeConcordance.keySet() ) {
            Genotype evalGenotype = eval.getGenotype(sample);
            Genotype truthGenotype = truth.getGenotype(sample);
            perSampleGenotypeConcordance.get(sample).update(evalGenotype,truthGenotype,alleleTruth);
            overallGenotypeConcordance.update(evalGenotype,truthGenotype,alleleTruth);
        }
    }

    private static double calculateNRD(GenotypeConcordanceTable table) {
        return calculateNRD(table.getTable());
    }

    private static double calculateNRD(int[][] concordanceCounts) {
        int correct = 0;
        int total = 0;
        correct += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HET.ordinal()];
        correct += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += correct;
        total += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HET.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HOM_REF.ordinal()];
        total += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HOM_REF.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HET.ordinal()];
        // NRD is by definition incorrec/total = 1.0-correct/total
        // note: if there are no observations (so the ratio is NaN), set this to 100%
        return total == 0 ? 1.0 : 1.0 - ( (double) correct)/( (double) total);
    }

    private static double calculateNRS(GenotypeConcordanceTable table) {
        return calculateNRS(table.getTable());
    }

    private static double calculateNRS(int[][] concordanceCounts) {
        long confirmedVariant = 0;
        long unconfirmedVariant = 0;
        for ( GenotypeType truthState : Arrays.asList(GenotypeType.HET,GenotypeType.HOM_VAR) ) {
            for ( GenotypeType evalState : GenotypeType.values() ) {
                if ( evalState == GenotypeType.MIXED )
                    continue;
                if ( evalState.equals(GenotypeType.HET) || evalState.equals(GenotypeType.HOM_VAR) )
                    confirmedVariant += concordanceCounts[evalState.ordinal()][truthState.ordinal()];
                else
                    unconfirmedVariant += concordanceCounts[evalState.ordinal()][truthState.ordinal()];
            }
        }

        long total = confirmedVariant + unconfirmedVariant;
        // note: if there are no observations (so the ratio is NaN) set this to 0%
        return total == 0l ? 0.0 : ( (double) confirmedVariant ) / ( (double) ( total ) );
    }


    class GenotypeConcordanceTable {

        private int[][] genotypeCounts;
        private int nMismatchingAlt;

        public GenotypeConcordanceTable() {
            genotypeCounts = new int[GenotypeType.values().length][GenotypeType.values().length];
            nMismatchingAlt = 0;
        }

        @Requires({"eval!=null","truth != null","truthAlleles != null"})
        public void update(Genotype eval, Genotype truth, Set<String> truthAlleles) {
            // this is slow but correct
            boolean matchingAlt = true;
            if ( eval.isCalled() && truth.isCalled() ) {
                // by default, no-calls "match" between alleles, so if
                // one or both sites are no-call or unavailable, the alt alleles match
                // otherwise, check explicitly: if the eval has an allele that's not ref, no-call, or present in truth
                // the alt allele is mismatching - regardless of whether the genotype is correct.
                for ( Allele evalAllele : eval.getAlleles() ) {
                    matchingAlt &= truthAlleles.contains(evalAllele.getBaseString());
                }
            }

            if ( matchingAlt ) {
                genotypeCounts[eval.getType().ordinal()][truth.getType().ordinal()]++;
            } else {
                nMismatchingAlt++;
            }
        }

        public int[][] getTable() {
            return genotypeCounts;
        }

        public int getnMismatchingAlt() {
            return nMismatchingAlt;
        }

        public int getnEvalGenotypes(GenotypeType type) {
            int nGeno = 0;
            for ( GenotypeType comptype : GenotypeType.values() )
                nGeno += genotypeCounts[type.ordinal()][comptype.ordinal()];
            return nGeno;
        }

        public int getnCalledEvalGenotypes() {
            int nGeno = 0;
            for ( GenotypeType evalType : Arrays.asList(GenotypeType.HOM_REF,GenotypeType.HOM_VAR,GenotypeType.HET) ) {
                nGeno += getnEvalGenotypes(evalType);
            }

            return nGeno + nMismatchingAlt;
        }

        public int getnCompGenotypes(GenotypeType type) {
            int nGeno = 0;
            for ( GenotypeType evaltype : GenotypeType.values() )
                nGeno += genotypeCounts[evaltype.ordinal()][type.ordinal()];
            return nGeno;
        }

        public int getnCalledCompGenotypes() {
            int nGeno = 0;
            for ( GenotypeType compType : Arrays.asList(GenotypeType.HOM_REF,GenotypeType.HOM_VAR,GenotypeType.HET) ) {
                nGeno += getnCompGenotypes(compType);
            }
            return nGeno;
        }

        public int get(GenotypeType evalType, GenotypeType compType) {
            return genotypeCounts[evalType.ordinal()][compType.ordinal()];
        }
    }

    class SiteConcordanceTable {

        private int[] siteConcordance;

        public SiteConcordanceTable() {
            siteConcordance = new int[SiteConcordanceType.values().length];
        }

        public void update(VariantContext evalVC, VariantContext truthVC) {
            SiteConcordanceType matchType = getMatchType(evalVC,truthVC);
            siteConcordance[matchType.ordinal()]++;
        }

        @Requires({"evalVC != null","truthVC != null"})
        private SiteConcordanceType getMatchType(VariantContext evalVC, VariantContext truthVC) {
            return SiteConcordanceType.getConcordanceType(evalVC,truthVC);
        }

        public int[] getSiteConcordance() {
            return siteConcordance;
        }

        public int get(SiteConcordanceType type) {
            return getSiteConcordance()[type.ordinal()];
        }
    }

    enum SiteConcordanceType {
        ALLELES_MATCH,
        EVAL_SUPERSET_TRUTH,
        EVAL_SUBSET_TRUTH,
        ALLELES_DO_NOT_MATCH,
        EVAL_ONLY,
        TRUTH_ONLY;

        public static SiteConcordanceType getConcordanceType(VariantContext eval, VariantContext truth) {
            if ( eval.isMonomorphicInSamples() )
                return TRUTH_ONLY;
            if ( truth.isMonomorphicInSamples() )
                return EVAL_ONLY;

            boolean evalSubsetTruth = VariantContextUtils.allelesAreSubset(eval,truth);
            boolean truthSubsetEval = VariantContextUtils.allelesAreSubset(truth,eval);

            if ( evalSubsetTruth && truthSubsetEval )
                return ALLELES_MATCH;

            if ( evalSubsetTruth )
                return EVAL_SUBSET_TRUTH;

            if ( truthSubsetEval )
                return EVAL_SUPERSET_TRUTH;

            return ALLELES_DO_NOT_MATCH;
        }
    }
}

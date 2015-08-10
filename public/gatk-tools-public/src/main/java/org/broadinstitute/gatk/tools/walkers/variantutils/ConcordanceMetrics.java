/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.tools.walkers.variantutils;

import com.google.java.contract.Requires;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFHeader;

import java.io.PrintStream;
import java.util.*;

/**
 * A class for tabulating and evaluating a callset-by-callset genotype concordance table
 * */
public class ConcordanceMetrics {

    final private Map<String,GenotypeConcordanceTable> perSampleGenotypeConcordance;
    final private GenotypeConcordanceTable overallGenotypeConcordance;
    final private SiteConcordanceTable overallSiteConcordance;
    final PrintStream sitesFile;

    public ConcordanceMetrics(VCFHeader evaluate, VCFHeader truth, PrintStream inputSitesFile) {
        HashSet<String> overlappingSamples = new HashSet<>(evaluate.getGenotypeSamples());
        overlappingSamples.retainAll(truth.getGenotypeSamples());
        perSampleGenotypeConcordance = new HashMap<>(overlappingSamples.size());
        for ( String sample : overlappingSamples ) {
            perSampleGenotypeConcordance.put(sample,new GenotypeConcordanceTable());
        }
        overallGenotypeConcordance = new GenotypeConcordanceTable();
        overallSiteConcordance = new SiteConcordanceTable();
        sitesFile = inputSitesFile;
        if (sitesFile != null) printSitesFileHeader();
    }

    private void printSitesFileHeader() {
        sitesFile.println("Locus\tSample\tTruth Genotype\tEval Genotype");
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
            throw new ReviewedGATKException("Attempted to request the concordance table for sample "+sample+" on which it was not calculated");
        return table;
    }

    public Map<String,GenotypeConcordanceTable> getPerSampleGenotypeConcordance() {
        return Collections.unmodifiableMap(perSampleGenotypeConcordance);
    }

    public Map<String,Double> getPerSampleNRD() {
        Map<String,Double> nrd = new HashMap<>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            nrd.put(sampleTable.getKey(),calculateNRD(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(nrd);
    }

    public Map<String,Double> getPerSampleOGC() {
        Map<String,Double> ogc = new HashMap<>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            ogc.put(sampleTable.getKey(),calculateOGC(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(ogc);
    }

    public Double getOverallNRD() {
        return calculateNRD(overallGenotypeConcordance);
    }

    public Double getOverallOGC() {
        return calculateOGC(overallGenotypeConcordance);
    }

    public Map<String,Double> getPerSampleNRS() {
        Map<String,Double> nrs = new HashMap<>(perSampleGenotypeConcordance.size());
        for ( Map.Entry<String,GenotypeConcordanceTable> sampleTable : perSampleGenotypeConcordance.entrySet() ) {
            nrs.put(sampleTable.getKey(),calculateNRS(sampleTable.getValue()));
        }

        return Collections.unmodifiableMap(nrs);
    }

    public Double getOverallNRS() {
        return calculateNRS(overallGenotypeConcordance);
    }

    @Requires({"eval != null","truth != null"})
    public void update(final VariantContext eval, final VariantContext truth) {
        overallSiteConcordance.update(eval,truth);
        final Set<Allele> truthAlleles = new HashSet<>(truth.getAlleles());
        for ( final String sample : perSampleGenotypeConcordance.keySet() ) {
            final Genotype evalGenotype = eval.getGenotype(sample);
            final Genotype truthGenotype = truth.getGenotype(sample);
            // ensure genotypes are either no-call ("."), missing (empty alleles), or diploid
            if ( ( ! evalGenotype.isNoCall() && evalGenotype.getPloidy() != 2 && evalGenotype.getPloidy() > 0) ||
                 ( ! truthGenotype.isNoCall() && truthGenotype.getPloidy() != 2 && truthGenotype.getPloidy() > 0) ) {
                throw new UserException(String.format("Concordance Metrics is currently only implemented for DIPLOID genotypes, found eval ploidy: %d, comp ploidy: %d",evalGenotype.getPloidy(),truthGenotype.getPloidy()));
            }
            final boolean allelesMatch = doAllelesMatch(evalGenotype, truthGenotype, truth.getReference(), truthAlleles);
            perSampleGenotypeConcordance.get(sample).update(allelesMatch, evalGenotype, truthGenotype);
            final boolean doPrint = overallGenotypeConcordance.update(allelesMatch, evalGenotype, truthGenotype);
            if(sitesFile != null && doPrint)
                sitesFile.println(eval.getChr() + ":" + eval.getStart() + "\t" + sample + "\t" + truthGenotype.getType() + "\t" + evalGenotype.getType());
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

    private static double calculateOGC(int[][] concordanceCounts) {
        int correct = 0;
        int total = 0;
        correct += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HOM_REF.ordinal()];
        correct += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HET.ordinal()];
        correct += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += correct;
        total += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HET.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_REF.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HOM_REF.ordinal()];
        total += concordanceCounts[GenotypeType.HET.ordinal()][GenotypeType.HOM_VAR.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HOM_REF.ordinal()];
        total += concordanceCounts[GenotypeType.HOM_VAR.ordinal()][GenotypeType.HET.ordinal()];
        // OGC is by definition correct/total
        // note: if there are no observations (so the ratio is NaN), set this to 100%
        return total == 0 ? 1.0 : ( (double) correct)/( (double) total);
    }

    private static double calculateNRS(GenotypeConcordanceTable table) {
        return calculateNRS(table.getTable());
    }

    private static double calculateOGC(GenotypeConcordanceTable table) {
        return calculateOGC(table.getTable());
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

    private boolean doAllelesMatch(final Genotype eval, final Genotype truth,
                                   final Allele truthRef, final Set<Allele> truthSiteAlleles) {
        // When determining if alleles match, there are a number of cases to consider.  In order:
        //  1) If either genotype is uncalled or unavailable, the alleles MATCH
        //  2) If the truth genotype is hom ref, then:
        //     a) If the truth variant is mononallelic (no alternate alleles), the alleles MATCH
        //     b) Otherwise, the alleles match IFF the alleles in the eval genotype are a subset
        //        of the alleles in the truth VARIANT
        //  3) Otherwise, the alleles match IFF the alleles in the eval genotype are a subset
        //     of the alleles in (the truth GENOTYPE + the truth REF allele)
        boolean matching = true;
        if (eval.isCalled() && truth.isCalled()) { // Case 1
            if (truth.isHomRef()) { // Case 2
                matching = truthSiteAlleles.size() == 1 || truthSiteAlleles.containsAll(eval.getAlleles());
            } else { // Case 3
                final Set<Allele> truthAlleles = new HashSet<>(truth.getAlleles());
                truthAlleles.add(truthRef);
                matching = truthAlleles.containsAll(eval.getAlleles());
            }
        }
        return matching;
    }


    class GenotypeConcordanceTable {

        private int[][] genotypeCounts;
        private int nMismatchingAlt;

        public GenotypeConcordanceTable() {
            genotypeCounts = new int[GenotypeType.values().length][GenotypeType.values().length];
            nMismatchingAlt = 0;
        }

        @Requires({"eval!=null","truth != null","truthAlleles != null"})
        public Boolean update(final boolean matchingAlt, final Genotype eval, final Genotype truth) {
            if ( matchingAlt ) {
                final int evalGT = eval.getType().ordinal();
                final int truthGT = truth.getType().ordinal();
                genotypeCounts[evalGT][truthGT]++;
                if(evalGT != truthGT)  //report variants where genotypes don't match
                    return true;
            } else {
                nMismatchingAlt++;
                return false;
                //return true; //alternatively, report variants where alt alleles don't match
            }
            return false;
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

            boolean evalSubsetTruth = GATKVariantContextUtils.allelesAreSubset(eval, truth);
            boolean truthSubsetEval = GATKVariantContextUtils.allelesAreSubset(truth, eval);

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

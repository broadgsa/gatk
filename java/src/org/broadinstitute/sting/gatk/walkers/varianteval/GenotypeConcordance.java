package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.playground.utils.report.tags.Analysis;
import org.broadinstitute.sting.playground.utils.report.tags.DataPoint;
import org.broadinstitute.sting.playground.utils.report.utils.TableType;
import org.broadinstitute.sting.utils.StingException;
import org.apache.log4j.Logger;

import java.util.*;

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

@Analysis(name = "Genotype Concordance", description = "Determine the genotype concordance between the genotypes in difference tracks")
public class GenotypeConcordance extends VariantEvaluator implements StandardEval {
    private static final boolean PRINT_INTERESTING_SITES = true;

    protected final static Logger logger = Logger.getLogger(GenotypeConcordance.class);

    // a mapping from allele count to stats
    @DataPoint(description = "the frequency statistics for each allele")
    FrequencyStats alleleFreqStats = new FrequencyStats();

    // a mapping from sample to stats
    @DataPoint(name="samples", description = "the concordance statistics for each sample")
    SampleStats sampleStats = null;

    // a mapping from sample to stats summary
    @DataPoint(name="summary", description = "the concordance statistics summary for each sample")
    SampleSummaryStats sampleSummaryStats = null;

    // two histograms of variant quality scores, for true positive and false positive calls
    @DataPoint(description = "the variant quality score histograms for true positive and false positive calls")
    QualityScoreHistograms qualityScoreHistograms = null;

    private static final int MAX_MISSED_VALIDATION_DATA = 100;

    private VariantEvalWalker.EvaluationContext group = null;

    static class FrequencyStats implements TableType {
        class Stats {
            public Stats(int found, int missed) { nFound = found; nMissed = missed; }
            public long nFound = 0;
            public long nMissed = 0;
        }
        public HashMap<Integer, Stats> alleleCountStats = new HashMap<Integer, Stats>();

        public Object[] getRowKeys() {
            String rows[] = new String[alleleCountStats.size()];
            int index = 0;
            for (int i : alleleCountStats.keySet()) rows[index++] = "Allele Count " + i;
            return rows;
        }

        public Object[] getColumnKeys() {
            return new String[]{"number_found", "number_missing"};
        }

        public String getName() {
            return "FrequencyStats";
        }

        public String getCell(int x, int y) {
            if (x >= alleleCountStats.size()) throw new IllegalStateException(x + " is greater than the max index of " + (alleleCountStats.size()-1));
            if (y == 0) return String.valueOf(alleleCountStats.get(alleleCountStats.keySet().toArray(new Integer[alleleCountStats.size()])[x]).nFound);
            else return String.valueOf(alleleCountStats.get(alleleCountStats.keySet().toArray(new Integer[alleleCountStats.size()])[x]).nMissed);
        }

        public void incrementFoundCount(int alleleFreq) {
            if (!alleleCountStats.containsKey(alleleFreq))
                alleleCountStats.put(alleleFreq,new Stats(1,0));
            else
                alleleCountStats.get(alleleFreq).nFound++;
        }

        public void incrementMissedCount(int alleleFreq) {
            if (!alleleCountStats.containsKey(alleleFreq))
                alleleCountStats.put(alleleFreq,new Stats(0,1));
            else
                alleleCountStats.get(alleleFreq).nMissed++;
        }
    }

    static class QualityScoreHistograms implements TableType {
        final static int NUM_BINS = 20;
        final HashMap<Integer,Integer> truePositiveQualityScoreMap = new HashMap<Integer,Integer>(); // A HashMap holds all the quality scores until we are able to bin them appropriately
        final HashMap<Integer,Integer> falsePositiveQualityScoreMap = new HashMap<Integer,Integer>();
        final int truePositiveHist[] = new int[NUM_BINS]; // the final histograms that get reported out
        final int falsePositiveHist[] = new int[NUM_BINS];
        final String[] rowKeys = new String[]{"true_positive_hist", "false_positive_hist"};

        public Object[] getRowKeys() {
            return rowKeys;
        }

        public Object[] getColumnKeys() {
            final String columnKeys[] = new String[NUM_BINS];
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                columnKeys[iii] = "histBin" + iii;
            }
            return columnKeys;
        }

        public String getName() {
            return "QualityScoreHistogram";
        }

        public String getCell(int x, int y) {
            if( x == 0 ) {
                return String.valueOf(truePositiveHist[y]);
            } else if ( x == 1 ) {
                return String.valueOf(falsePositiveHist[y]);
            } else {
                throw new StingException( "Unknown row in " + getName() + ", row = " + x );
            }
        }

        public String toString() {
            String returnString = "";
            // output both histogram arrays
            returnString += "TP: ";
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                returnString += truePositiveHist[iii] + " ";
            }
            returnString += "\nFP: ";
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                returnString += falsePositiveHist[iii] + " ";
            }
            return returnString;
        }

        public void incrValue( final double qual, final boolean isTruePositiveCall ) {
            HashMap<Integer,Integer> qualScoreMap;
            if( isTruePositiveCall ) {
                qualScoreMap = truePositiveQualityScoreMap;
            } else {
                qualScoreMap = falsePositiveQualityScoreMap;
            }
            final Integer qualKey = Math.round((float) qual);
            if( qualScoreMap.containsKey(qualKey) ) {
                qualScoreMap.put(qualKey, qualScoreMap.get(qualKey) + 1);
            } else {
                qualScoreMap.put(qualKey, 1);
            }
        }

        public void organizeHistogramTables() {
            for( int iii = 0; iii < NUM_BINS; iii++ ) {
                truePositiveHist[iii] = 0;
                falsePositiveHist[iii] = 0;
            }

            int maxQual = 0;

            // Calculate the maximum quality score for both TP and FP calls in order to normalize and histogram
            for( final Integer qual : truePositiveQualityScoreMap.keySet()) {
                if( qual > maxQual ) {
                    maxQual = qual;
                }
            }
            for( final Integer qual : falsePositiveQualityScoreMap.keySet()) {
                if( qual > maxQual ) {
                    maxQual = qual;
                }
            }

            final double binSize = ((double)maxQual) / ((double) (NUM_BINS-1)); //BUGBUG: should be normalized max to min, not max to 0

            for( final Integer qual : truePositiveQualityScoreMap.keySet()) {
                final int index = (int)Math.floor( ((double)qual) / binSize );
                if(index >= 0) { //BUGBUG: problem when maxQual is zero?
                    truePositiveHist[ index ] += truePositiveQualityScoreMap.get(qual);
                }
            }
            for( final Integer qual : falsePositiveQualityScoreMap.keySet()) {
                final int index = (int)Math.floor( ((double)qual) / binSize );
                if(index >= 0) {
                    falsePositiveHist[ index ] += falsePositiveQualityScoreMap.get(qual);
                }
            }
        }
    }

    // keep a list of the validation data we saw before the first eval data
    private HashSet<VariantContext> missedValidationData = new HashSet<VariantContext>();


    public GenotypeConcordance(VariantEvalWalker parent) {
        super(parent);
    }

    public String getName() {
        return "genotypeConcordance";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public boolean enabled() {
        return true;
    }

    public String toString() {
        return getName() + ": <table>";
    }

    private boolean warnedAboutValidationData = false;

    public String update2(VariantContext eval, VariantContext validation, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context, VariantEvalWalker.EvaluationContext group) {
        this.group = group;

        String interesting = null;

        // sanity check that we at least have either eval or validation data
        if (eval == null && !isValidVC(validation)) {
            return interesting;
        }

        if( qualityScoreHistograms == null ) {
            qualityScoreHistograms = new QualityScoreHistograms();
        }

        if (sampleStats == null) {
            if (eval != null) {
                // initialize the concordance table
                sampleStats = new SampleStats(eval,Genotype.Type.values().length);
                sampleSummaryStats = new SampleSummaryStats(eval);
                for (final VariantContext vc : missedValidationData) {
                    determineStats(null, vc);
                }
                missedValidationData = null;
            } else {
                // todo -- Eric, this results in a memory problem when eval is WEx data but you are using CG calls genome-wide
                // todo -- perhaps you need should extend the evaluators with an initialize
                // todo -- method that gets the header (or samples) for the first eval sites?
                if (missedValidationData.size() > MAX_MISSED_VALIDATION_DATA) {
                    if (!warnedAboutValidationData) {
                        logger.warn("Too many genotype sites missed before eval site appeared; ignoring");
                        warnedAboutValidationData = true;
                    }
                } else {
                    missedValidationData.add(validation);
                }
                return interesting;
            }
        }

        interesting = determineStats(eval, validation);

        return interesting; // we don't capture any interesting sites
    }

    private String determineStats(final VariantContext eval, final VariantContext validation) {

        String interesting = null;

        final boolean validationIsValidVC = isValidVC(validation);

        // determine concordance for eval data
        if (eval != null) {
           for (final String sample : eval.getSampleNames()) {
                final Genotype.Type called = eval.getGenotype(sample).getType();
                final Genotype.Type truth;

                if (!validationIsValidVC || !validation.hasGenotype(sample)) {
                    truth = Genotype.Type.NO_CALL;
                } else {
                    truth = validation.getGenotype(sample).getType();
                    // interesting = "ConcordanceStatus=FP";
                }

                sampleStats.incrValue(sample, truth, called);
            }
        }
        // otherwise, mark no-calls for all samples
        else {
            final Genotype.Type called = Genotype.Type.NO_CALL;

            for (final String sample : validation.getSampleNames()) {
                final Genotype.Type truth = validation.getGenotype(sample).getType();
                sampleStats.incrValue(sample, truth, called);
                if ( (truth == Genotype.Type.HOM_VAR || truth == Genotype.Type.HET) && called == Genotype.Type.NO_CALL ) {
                    if ( PRINT_INTERESTING_SITES && super.getVEWalker().printInterestingSites() )
                        System.out.printf("%s: HM3 FN => %s%n", group, validation);
                }
            }
        }

        // determine allele count concordance ()
        if (validationIsValidVC && validation.isPolymorphic()) {
            int trueAlleleCount = 0;
            for (final Allele a : validation.getAlternateAlleles()) {
                trueAlleleCount += validation.getChromosomeCount(a);
            }
            if (eval != null) {
                alleleFreqStats.incrementFoundCount(trueAlleleCount);
            } else {
                alleleFreqStats.incrementMissedCount(trueAlleleCount);
            }
        }

        // TP & FP quality score histograms
        if( eval != null && eval.isPolymorphic() && validationIsValidVC ) {
            if( eval.getSampleNames().size() == 1 ) { // single sample calls
                for( final String sample : eval.getSampleNames() ) { // only one sample
                    if( validation.hasGenotype(sample) ) {
                        final Genotype truth = validation.getGenotype(sample);
                        qualityScoreHistograms.incrValue( eval.getPhredScaledQual(), !truth.isHomRef() );
                    }
                }
            } else { // multi sample calls
                qualityScoreHistograms.incrValue( eval.getPhredScaledQual(), validation.isPolymorphic() );
            }

        }

        return interesting;
    }

    private static boolean isValidVC(final VariantContext vc) {
        return (vc != null && !vc.isFiltered());
    }

    public void finalizeEvaluation() {
        if( qualityScoreHistograms != null ) {
            qualityScoreHistograms.organizeHistogramTables();
        }

        if( sampleSummaryStats != null && sampleStats != null ) {
            sampleSummaryStats.generateSampleSummaryStats( sampleStats );
        }
    }
}

/**
 * a table of sample names to genotype concordance figures
 */
class SampleStats implements TableType {
    private final int nGenotypeTypes;

    // sample to concordance stats object
    public final HashMap<String, long[][]> concordanceStats = new HashMap<String, long[][]>();

    /**
     *
     * @return one row per sample
     */
    public Object[] getRowKeys() {
        return concordanceStats.keySet().toArray(new String[concordanceStats.size()]);
    }

    /**
     * increment the specified value
     * @param sample the sample name
     * @param truth the truth type
     * @param called the called type
     */
    public void incrValue(String sample, Genotype.Type truth, Genotype.Type called) {
        if ( concordanceStats.containsKey(sample) )
            concordanceStats.get(sample)[truth.ordinal()][called.ordinal()]++;
        else if ( called != Genotype.Type.NO_CALL )
            throw new StingException("Sample " + sample + " has not been seen in a previous eval; this analysis module assumes that all samples are present in each variant context");
    }

    /**
     * get the column keys
     * @return a list of objects, in this case strings, that are the column names
     */
    public Object[] getColumnKeys() {
        return new String[]{"total_true_ref","%_ref/ref","n_ref/no-call",
                            "n_ref/ref","n_ref/het","n_ref/hom",
                            "total_true_het","%_het/het","n_het/no-call",
                            "n_het/ref","n_het/het","n_het/hom",
                            "total_true_hom","%_hom/hom","n_hom/no-call",
                            "n_hom/ref","n_hom/het","n_hom/hom"};
    }

    public SampleStats(VariantContext vc, int nGenotypeTypes) {
        this.nGenotypeTypes = nGenotypeTypes;
        for (String sample : vc.getSampleNames())
            concordanceStats.put(sample, new long[nGenotypeTypes][nGenotypeTypes]);
    }

    public Object getCell(int x, int y) {
        // we have three rows of 6 right now for output (rows: ref, het, hom)
        Genotype.Type type = Genotype.Type.values()[(y/6)+1]; // get the row type
        // save some repeat work, get the total every time
        long total = 0;
        Object[] rowKeys = getRowKeys();
        for (int called = 0; called < nGenotypeTypes; called++) {
            total += concordanceStats.get(rowKeys[x])[type.ordinal()][called];
        }

        // now get the cell they're interested in
        switch (y % 6) {
            case (0): // get the total_true for this type
                return total;
            case (1):
                return total == 0 ? 0.0 : (100.0 * (double) concordanceStats.get(rowKeys[x])[type.ordinal()][type.ordinal()] / (double) total);
            default:
                return concordanceStats.get(rowKeys[x])[type.ordinal()][(y % 6) - 2];
        }
    }

    public String getName() {
        return "Sample Statistics";
    }
}

/**
 * a table of sample names to genotype concordance summary statistics
 */
class SampleSummaryStats implements TableType {
    private final static String ALL_SAMPLES_KEY = "allSamples";
    private final static String[] COLUMN_KEYS = new String[]{
            "percent_comp_ref_called_var",
            "percent_comp_het_called_het",
            "percent_comp_het_called_var",
            "percent_comp_hom_called_hom",
            "percent_comp_hom_called_var",
            "percent_non-reference_sensitivity",
            "percent_overall_genotype_concordance",
            "percent_non-reference_discrepancy_rate"};

    // sample to concordance stats object
    private final HashMap<String, double[]> concordanceSummary = new HashMap<String, double[]>();

    /**
     *
     * @return one row per sample
     */
    public Object[] getRowKeys() {
        return concordanceSummary.keySet().toArray(new String[concordanceSummary.size()]);
    }

    /**
     * get the column keys
     * @return a list of objects, in this case strings, that are the column names
     */
    public Object[] getColumnKeys() {
        return COLUMN_KEYS;
    }

    public SampleSummaryStats(final VariantContext vc) {
        concordanceSummary.put(ALL_SAMPLES_KEY, new double[COLUMN_KEYS.length]);
        for( final String sample : vc.getSampleNames() ) {
            concordanceSummary.put(sample, new double[COLUMN_KEYS.length]);
        }
    }

    public Object getCell(int x, int y) {
        final Object[] rowKeys = getRowKeys();
        return String.format("%.2f",concordanceSummary.get(rowKeys[x])[y]);
    }

    /**
     * Helper routine that sums up all columns / rows found in stats specified by all pairs in d1 x d2
     *
     * @param stats
     * @param d1
     * @param d2
     * @return
     */
    private long sumStatsAllPairs( final long[][] stats, EnumSet<Genotype.Type> d1, EnumSet<Genotype.Type> d2 ) {
        long sum = 0L;

        for(final Genotype.Type e1 : d1 ) {
            for(final Genotype.Type e2 : d2 ) {
                sum += stats[e1.ordinal()][e2.ordinal()];
            }
        }

        return sum;
    }

    private long sumStatsDiag( final long[][] stats, EnumSet<Genotype.Type> d1) {
        long sum = 0L;

        for(final Genotype.Type e1 : d1 ) {
            sum += stats[e1.ordinal()][e1.ordinal()];
        }

        return sum;
    }

    private double ratio(long numer, long denom) {
        return denom != 0L ? 100.0 * ( ((double)numer) / ((double)denom) ) : 0.0;
    }

    final long[] allSamplesNumerators = new long[COLUMN_KEYS.length];
    final long[] allSamplesDenominators = new long[COLUMN_KEYS.length];

    private void updateSummaries(int i, double[] summary, long numer, long denom ) {
        allSamplesNumerators[i] += numer;
        allSamplesDenominators[i] += denom;
        summary[i] = ratio(numer, denom);
    }


    /**
     * Calculate the five summary stats per sample
     * @param sampleStats The Map which holds concordance values per sample
     */
    public void generateSampleSummaryStats( final SampleStats sampleStats ) {
        EnumSet<Genotype.Type> allVariantGenotypes = EnumSet.of(Genotype.Type.HOM_VAR, Genotype.Type.HET);
        EnumSet<Genotype.Type> allCalledGenotypes = EnumSet.of(Genotype.Type.HOM_VAR, Genotype.Type.HET, Genotype.Type.HOM_REF);
        EnumSet<Genotype.Type> allGenotypes = EnumSet.allOf(Genotype.Type.class);

        for( final String sample : concordanceSummary.keySet() ) {
            if ( sample.equals(ALL_SAMPLES_KEY) ) continue;

            final long[][] stats = sampleStats.concordanceStats.get(sample);
            final double[] summary = concordanceSummary.get(sample);
            if( stats == null ) { throw new StingException( "SampleStats and SampleSummaryStats contain different samples! sample = " + sample ); }

            long numer, denom;

            // Summary 0: % ref called as var
            numer = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_REF), allVariantGenotypes);
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_REF), allGenotypes);
            updateSummaries(0, summary,  numer, denom);

            // Summary 1: % het called as het
            numer = stats[Genotype.Type.HET.ordinal()][Genotype.Type.HET.ordinal()];
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HET), allGenotypes);
            updateSummaries(1, summary,  numer, denom);

            // Summary 2: % het called as var
            numer = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HET), allVariantGenotypes);
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HET), allGenotypes);
            updateSummaries(2, summary,  numer, denom);

            // Summary 3: % homVar called as homVar
            numer = stats[Genotype.Type.HOM_VAR.ordinal()][Genotype.Type.HOM_VAR.ordinal()];
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_VAR), allGenotypes);
            updateSummaries(3, summary,  numer, denom);

            // Summary 4: % homVars called as var
            numer = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_VAR), allVariantGenotypes);
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_VAR), allGenotypes);
            updateSummaries(4, summary,  numer, denom);

            // Summary 5: % non-ref called as non-ref
            // MAD: this is known as the non-reference sensitivity (# non-ref according to comp found in eval / # non-ref in comp)
            numer = sumStatsAllPairs(stats, allVariantGenotypes, allVariantGenotypes);
            denom = sumStatsAllPairs(stats, allVariantGenotypes, allGenotypes);
            updateSummaries(5, summary,  numer, denom);

            // Summary 6: overall genotype concordance of sites called in eval track
            // MAD: this is the tradition genotype concordance
            numer = sumStatsDiag(stats, allCalledGenotypes);
            denom = sumStatsAllPairs(stats, allCalledGenotypes, allCalledGenotypes);
            updateSummaries(6, summary,  numer, denom);

            // Summary 7: overall genotype concordance of sites called non-ref in eval track
            long homrefConcords = stats[Genotype.Type.HOM_REF.ordinal()][Genotype.Type.HOM_REF.ordinal()];
            long diag = sumStatsDiag(stats, allVariantGenotypes);
            long allNoHomRef = sumStatsAllPairs(stats, allCalledGenotypes, allCalledGenotypes) - homrefConcords;
            numer = allNoHomRef - diag;
            denom = allNoHomRef;
            updateSummaries(7, summary,  numer, denom);
        }

        // update the final summary stats
        final double[] allSamplesSummary = concordanceSummary.get(ALL_SAMPLES_KEY);
        for ( int i = 0; i < allSamplesSummary.length; i++) {
            allSamplesSummary[i] = ratio(allSamplesNumerators[i], allSamplesDenominators[i]);
        }

    }

    public String getName() {
        return "Sample Summary Statistics";
    }
}

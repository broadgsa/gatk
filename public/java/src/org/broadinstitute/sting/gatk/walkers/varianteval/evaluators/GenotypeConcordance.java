package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
public class GenotypeConcordance extends VariantEvaluator {
    private static final boolean PRINT_INTERESTING_SITES = true;

    protected final static Logger logger = Logger.getLogger(GenotypeConcordance.class);

    // a mapping from sample to stats
    @DataPoint(description = "the detailed concordance statistics for each sample")
    SampleStats detailedStats = null;

    // a mapping from sample to stats summary
    @DataPoint(description = "the simplified concordance statistics for each sample")
    SampleSummaryStats simplifiedStats = null;

    private static final int MAX_MISSED_VALIDATION_DATA = 100;

    private boolean discordantInteresting = false;

    static class FrequencyStats implements TableType {
        class Stats {
            public Stats(int found, int missed) { nFound = found; nMissed = missed; }
            public long nFound = 0;
            public long nMissed = 0;
        }
        public HashMap<Integer, Stats> foundMissedByAC = new HashMap<Integer, Stats>();

        public Object[] getRowKeys() {
            String rows[] = new String[foundMissedByAC.size()];
            int index = 0;
            for (int i : foundMissedByAC.keySet()) rows[index++] = "AlleleCount_" + i;
            return rows;
        }

        public Object[] getColumnKeys() {
            return new String[]{"number_found", "number_missing"};
        }

        public String getName() {
            return "FrequencyStats";
        }

        public String getCell(int x, int y) {
            if (x >= foundMissedByAC.size()) throw new IllegalStateException(x + " is greater than the max index of " + (foundMissedByAC.size()-1));
            if (y == 0) return String.valueOf(foundMissedByAC.get(foundMissedByAC.keySet().toArray(new Integer[foundMissedByAC.size()])[x]).nFound);
            else return String.valueOf(foundMissedByAC.get(foundMissedByAC.keySet().toArray(new Integer[foundMissedByAC.size()])[x]).nMissed);
        }

        public void incrementFoundCount(int alleleFreq) {
            if (!foundMissedByAC.containsKey(alleleFreq))
                foundMissedByAC.put(alleleFreq,new Stats(1,0));
            else
                foundMissedByAC.get(alleleFreq).nFound++;
        }

        public void incrementMissedCount(int alleleFreq) {
            if (!foundMissedByAC.containsKey(alleleFreq))
                foundMissedByAC.put(alleleFreq,new Stats(0,1));
            else
                foundMissedByAC.get(alleleFreq).nMissed++;
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
                throw new ReviewedStingException( "Unknown row in " + getName() + ", row = " + x );
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


    //public GenotypeConcordance(VariantEvalWalker parent) {
    //    super(parent);
	//	discordantInteresting = parent.DISCORDANT_INTERESTING;
    //}

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

    public String update2(VariantContext eval, VariantContext validation, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        String interesting = null;

        // sanity check that we at least have either eval or validation data
        if ( (validation != null && !validation.hasGenotypes()) || eval == null && !isValidVC(validation)) {
            return interesting;
        }

        if (detailedStats == null) {
            if (eval != null) {
                // initialize the concordance table
                detailedStats = new SampleStats(eval,Genotype.Type.values().length);
                simplifiedStats = new SampleSummaryStats(eval);
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
                        //logger.warn("Too many genotype sites missed before eval site appeared; ignoring");
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
        final String evalAC = ( vcHasGoodAC(eval) ) ? String.format("evalAC%d",getAC(eval)) : null ;
        final String validationAC = ( vcHasGoodAC(validation) ) ? String.format("compAC%d",getAC(validation)) : null;

        // determine concordance for eval data
        if (eval != null) {
           for (final String sample : eval.getGenotypes().keySet()) {
                final Genotype.Type called = eval.getGenotype(sample).getType();
                final Genotype.Type truth;

                if (!validationIsValidVC || !validation.hasGenotype(sample)) {
                    truth = Genotype.Type.NO_CALL;
                } else {
                    truth = validation.getGenotype(sample).getType();
                    // interesting = "ConcordanceStatus=FP";
                    if (discordantInteresting && truth.ordinal() != called.ordinal())
                    {
                        interesting = "ConcordanceStatus=" + truth + "/" + called;
                    }
                }

                detailedStats.incrValue(sample, truth, called);
            }
        }
        // otherwise, mark no-calls for all samples
        else {
            final Genotype.Type called = Genotype.Type.NO_CALL;

            for (final String sample : validation.getGenotypes().keySet()) {
                final Genotype.Type truth = validation.getGenotype(sample).getType();
                detailedStats.incrValue(sample, truth, called);

                // print out interesting sites
                /*
                if ( PRINT_INTERESTING_SITES && super.getVEWalker().gcLog != null ) {
                    if ( (truth == Genotype.Type.HOM_VAR || truth == Genotype.Type.HET) && called == Genotype.Type.NO_CALL ) {
                        super.getVEWalker().gcLog.printf("%s FN %s%n", group, validation);
                    }
                    if ( (called == Genotype.Type.HOM_VAR || called == Genotype.Type.HET) && truth == Genotype.Type.HOM_REF ) {
                        super.getVEWalker().gcLog.printf("%s FP %s%n", group, validation);
                    }
                }
                */
            }
        }

        return interesting;
    }

    private static boolean isValidVC(final VariantContext vc) {
        return (vc != null && !vc.isFiltered());
    }

    public void finalizeEvaluation() {
        if( simplifiedStats != null && detailedStats != null ) {
            simplifiedStats.generateSampleSummaryStats(detailedStats);
        }
    }

    private boolean vcHasGoodAC(VariantContext vc) {
        return ( vc != null && vc.getAlternateAlleles().size() == 1 && vc.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) );

    }

    private int getAC(VariantContext vc) {
        if ( List.class.isAssignableFrom(vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY).getClass()) ) {
            return ((List<Integer>) vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY)).get(0);
        } else if ( Integer.class.isAssignableFrom(vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY).getClass())) {
            return (Integer) vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY);
        } else if ( String.class.isAssignableFrom(vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY).getClass()) ) {
            // two ways of parsing
            String ac = (String) vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY);
            if ( ac.startsWith("[") ) {
                return Integer.parseInt(ac.replaceAll("\\[","").replaceAll("\\]",""));
            } else {
                try {
                    return Integer.parseInt(ac);
                } catch ( NumberFormatException e ) {
                    throw new UserException(String.format("The format of the AC field is improperly formatted: AC=%s",ac));
                }
            }
        } else {
            throw new UserException(String.format("The format of the AC field does not appear to be of integer-list or String format, class was %s",vc.getAttribute(VCFConstants.ALLELE_COUNT_KEY).getClass()));
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
            throw new UserException.CommandLineException("Sample " + sample + " has not been seen in a previous eval; this analysis module assumes that all samples are present in each variant context");
    }

    /**
     * get the column keys
     * @return a list of objects, in this case strings, that are the column names
     */
    public Object[] getColumnKeys() {
//        return new String[]{"total_true_ref","%_ref/ref","n_ref/no-call",
//                            "n_ref/ref","n_ref/het","n_ref/hom",
//                            "total_true_het","%_het/het","n_het/no-call",
//                            "n_het/ref","n_het/het","n_het/hom",
//                            "total_true_hom","%_hom/hom","n_hom/no-call",
//                            "n_hom/ref","n_hom/het","n_hom/hom"};
        return new String[]{"total_true_ref","pct_ref_vs_ref","n_ref_vs_no_call",
                "n_ref_vs_ref","n_ref_vs_het","n_ref_vs_hom",
                "total_true_het","pct_het_vs_het","n_het_vs_no_call",
                "n_het_vs_ref","n_het_vs_het","n_het_vs_hom",
                "total_true_hom","pct_hom_vs_hom","n_hom_vs_no_call",
                "n_hom_vs_ref","n_hom_vs_het","n_hom_vs_hom"};
    }


    public SampleStats(VariantContext vc, int nGenotypeTypes) {
        this.nGenotypeTypes = nGenotypeTypes;
        for (String sample : vc.getGenotypes().keySet())
            concordanceStats.put(sample, new long[nGenotypeTypes][nGenotypeTypes]);
    }

    public SampleStats(int genotypeTypes) {
        nGenotypeTypes = genotypeTypes;
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
 * Sample stats, but for AC
 */
class ACStats extends SampleStats {
    private String[] rowKeys;

    public ACStats(VariantContext evalvc, VariantContext compvc, int nGenotypeTypes) {
        super(nGenotypeTypes);
        rowKeys = new String[1+2*evalvc.getGenotypes().size()+1+2*compvc.getGenotypes().size()];
        for ( int i = 0; i <= 2*evalvc.getGenotypes().size(); i++ ) { // todo -- assuming ploidy 2 here...
            concordanceStats.put(String.format("evalAC%d",i),new long[nGenotypeTypes][nGenotypeTypes]);
            rowKeys[i] = String.format("evalAC%d",i);

        }

        for ( int i = 0; i <= 2*compvc.getGenotypes().size(); i++ ) {
            concordanceStats.put(String.format("compAC%d",i), new long[nGenotypeTypes][nGenotypeTypes]);
            rowKeys[1+2*evalvc.getGenotypes().size()+i] = String.format("compAC%d",i);
        }
    }

    public String getName() {
        return "Allele Count Statistics";
    }

    public Object[] getRowKeys() {
        if ( rowKeys == null ) {
            throw new StingException("RowKeys is null!");
        }
        return rowKeys;
    }
}

/**
 * a table of sample names to genotype concordance summary statistics
 */
class SampleSummaryStats implements TableType {
    protected final static String ALL_SAMPLES_KEY = "allSamples";
    protected final static String[] COLUMN_KEYS = new String[]{
            "percent_comp_ref_called_ref",
            "percent_comp_het_called_het",
            "percent_comp_hom_called_hom",
            "percent_non_reference_sensitivity",
            "percent_overall_genotype_concordance",
            "percent_non_reference_discrepancy_rate"};

    // sample to concordance stats object
    protected final HashMap<String, double[]> concordanceSummary = new HashMap<String, double[]>();

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
        for( final String sample : vc.getGenotypes().keySet() ) {
            concordanceSummary.put(sample, new double[COLUMN_KEYS.length]);
        }
    }

    public SampleSummaryStats() {

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
            if( stats == null ) { throw new ReviewedStingException( "SampleStats and SampleSummaryStats contain different samples! sample = " + sample ); }

            long numer, denom;

            // Summary 0: % ref called as ref
            numer = stats[Genotype.Type.HOM_REF.ordinal()][Genotype.Type.HOM_REF.ordinal()];
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_REF), allGenotypes);
            updateSummaries(0, summary,  numer, denom);

            // Summary 1: % het called as het
            numer = stats[Genotype.Type.HET.ordinal()][Genotype.Type.HET.ordinal()];
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HET), allGenotypes);
            updateSummaries(1, summary,  numer, denom);

            // Summary 2: % homVar called as homVar
            numer = stats[Genotype.Type.HOM_VAR.ordinal()][Genotype.Type.HOM_VAR.ordinal()];
            denom = sumStatsAllPairs(stats, EnumSet.of(Genotype.Type.HOM_VAR), allGenotypes);
            updateSummaries(2, summary,  numer, denom);

            // Summary 3: % non-ref called as non-ref
            // MAD: this is known as the non-reference sensitivity (# non-ref according to comp found in eval / # non-ref in comp)
            numer = sumStatsAllPairs(stats, allVariantGenotypes, allVariantGenotypes);
            denom = sumStatsAllPairs(stats, allVariantGenotypes, allGenotypes);
            updateSummaries(3, summary,  numer, denom);

            // Summary 4: overall genotype concordance of sites called in eval track
            // MAD: this is the tradition genotype concordance
            numer = sumStatsDiag(stats, allCalledGenotypes);
            denom = sumStatsAllPairs(stats, allCalledGenotypes, allCalledGenotypes);
            updateSummaries(4, summary,  numer, denom);

            // Summary 5: overall genotype concordance of sites called non-ref in eval track
            long homrefConcords = stats[Genotype.Type.HOM_REF.ordinal()][Genotype.Type.HOM_REF.ordinal()];
            long diag = sumStatsDiag(stats, allVariantGenotypes);
            long allNoHomRef = sumStatsAllPairs(stats, allCalledGenotypes, allCalledGenotypes) - homrefConcords;
            numer = allNoHomRef - diag;
            denom = allNoHomRef;
            updateSummaries(5, summary,  numer, denom);
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

/**
 * SampleSummaryStats .. but for allele counts
 */
class ACSummaryStats extends SampleSummaryStats {
    private String[] rowKeys;

    public ACSummaryStats (final VariantContext evalvc, final VariantContext compvc) {
        concordanceSummary.put(ALL_SAMPLES_KEY, new double[COLUMN_KEYS.length]);
        rowKeys = new String[3+2*evalvc.getGenotypes().size() + 2*compvc.getGenotypes().size()];
        rowKeys[0] = ALL_SAMPLES_KEY;
        for( int i = 0; i <= 2*evalvc.getGenotypes().size() ; i ++ ) {
            concordanceSummary.put(String.format("evalAC%d",i), new double[COLUMN_KEYS.length]);
            rowKeys[i+1] = String.format("evalAC%d",i);
        }
        for( int i = 0; i <= 2*compvc.getGenotypes().size() ; i ++ ) {
            concordanceSummary.put(String.format("compAC%d",i), new double[COLUMN_KEYS.length]);
            rowKeys[2+2*evalvc.getGenotypes().size()+i] = String.format("compAC%d",i);
        }

    }

    public String getName() {
        return "Allele Count Summary Statistics";
    }

    public Object[] getRowKeys() {
        if ( rowKeys == null) {
            throw new StingException("rowKeys is null!!");
        }
        return rowKeys;
    }
}

class CompACNames implements Comparator{

    final Logger myLogger;
    private boolean info = true;

    public CompACNames(Logger l) {
        myLogger = l;
    }

    public boolean equals(Object o) {
        return ( o.getClass() == CompACNames.class );
    }

    public int compare(Object o1, Object o2) {
        if ( info ) {
            myLogger.info("Sorting AC names");
            info = false;
        }
        //System.out.printf("Objects %s %s get ranks %d %d%n",o1.toString(),o2.toString(),getRank(o1),getRank(o2));
        return getRank(o1) - getRank(o2);
    }

    public int getRank(Object o) {
        if ( o.getClass() != String.class ) {
            return Integer.MIN_VALUE/4;
        } else {
            String s = (String) o;
            if ( s.startsWith("eval") ) {
                return Integer.MIN_VALUE/4 + 1 + parseAC(s);
            } else if ( s.startsWith("comp") ) {
                return 1+ parseAC(s);
            } else {
                return Integer.MIN_VALUE/4;
            }
        }
    }

    public int parseAC(String s) {
        String[] g = s.split("AC");
        return Integer.parseInt(g[1]);
    }
}


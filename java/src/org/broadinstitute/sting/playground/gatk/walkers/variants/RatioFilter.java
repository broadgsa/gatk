package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.gatk.refdata.TabularROD;
import org.broadinstitute.sting.utils.*;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import edu.mit.broad.picard.util.MathUtil;

class GenotypeFeatureData extends ArrayList<Double> {
    public enum Tail { LeftTailed, RightTailed, TwoTailed }

    private String genotype = null;
    private Integer[] permutation = null;
    private Logger logger = null;

    private Tail tail = null;

    double lowThreshold = -1;
    double highThreshold = -1;

    public GenotypeFeatureData(final String genotype, Logger logger, Tail tail ) {
        super();
        this.genotype = genotype;
        this.logger = logger;
        this.tail = tail;
    }

    public Tail getTail() {
        return tail;
    }

    public double getLowThreshold() {
        return lowThreshold;
    }

    public void setLowThreshold(double lowThreshold) {
        this.lowThreshold = lowThreshold;
    }

    public double getHighThreshold() {
        return highThreshold;
    }

    public void setHighThreshold(double highThreshold) {
        this.highThreshold = highThreshold;
    }

    public void finalizeData() {
        permutation = Utils.SortPermutation(this);
    }

    public void determineThresholds(final double pvalueLimit) {
        if ( pvalueLimit <= 0.0 || pvalueLimit >= 1.0 )
            throw new RuntimeException(String.format("Invalid pValue limit = %f", pvalueLimit));

        switch ( tail ) {
            case LeftTailed:
                lowThreshold = highThreshold = determineOneTailThreshold(pvalueLimit);
                break;
            case RightTailed:
                lowThreshold = highThreshold = determineOneTailThreshold(1 - pvalueLimit);
                break;
            case TwoTailed:
                lowThreshold = determineOneTailThreshold(pvalueLimit/2);
                highThreshold = determineOneTailThreshold(1-pvalueLimit/2);
                break;
        }

        //logger.info(String.format("  %d of %d elements are >= threshold = %f (%f percent)",
        //        nElementsAboveThreshold(), nElements(), getHighThreshold(), nElementsAboveThreshold() / ( 1.0 * nElements())));
        //logger.info(String.format("  %d of %d elements are <= threshold = %f (%f percent)",
        //        nElementsBelowThreshold(), nElements(), getLowThreshold(), nElementsBelowThreshold() / ( 1.0 * nElements())));
    }

    private double determineOneTailThreshold(final double pvalueLimit) {
        double trueSortedIndexLimit = pvalueLimit * nElements();
        int sortedIndexLimit = (int)(pvalueLimit > 0.5 ? Math.ceil(trueSortedIndexLimit) : Math.floor(trueSortedIndexLimit));
        double threshold = get(itemIndex(sortedIndexLimit));

        logger.debug(String.format("determineTailThreshold(%s, %f) => %f => %d => %f", genotype, pvalueLimit, trueSortedIndexLimit, sortedIndexLimit, threshold));

        return threshold;
    }

    public int nElementsBelowThreshold() { return nElementsPassingThreshold(true, false); }
    public int nElementsAboveThreshold() { return nElementsPassingThreshold(false, true); }

    public int nElementsPassingThreshold(boolean keepBelow, boolean keepAbove) {
        int count = 0;

        for ( double stat : this ) {
            if ( passesThreshold(stat, keepBelow, keepAbove) ) count++;
        }

        return count;
    }

    public boolean passesThreshold(double value, boolean keepBelow, boolean keepAbove) {
        //System.out.printf("passThreshold %s, %b, %b = %b%n", value, keepBelow, keepAbove, value < threshold && keepBelow || value >= threshold && keepAbove);
        return value <= getHighThreshold() && keepBelow || value >= getLowThreshold() && keepAbove;
    }

    public boolean passesThreshold(double value) {
        switch ( tail ) {
            case LeftTailed:
                return value >= getLowThreshold();
            case RightTailed:
                return value <= getHighThreshold();
            case TwoTailed:
                return value >= getLowThreshold() && value < getHighThreshold();
            default:
                return true;
        }
    }

    public double pValue(double value) {
        return 1.0;
    }

    public int itemIndex(int sortedIndex) {
        return permutation[sortedIndex];
    }

    public int nElements() { return size(); }

    public String toString() {
        return String.format("[GenotypeFeatureData: genotype=%s, tail=%s, lowThreshold=%f, highThreshold=%f]", genotype, tail, lowThreshold, highThreshold);
    }
}


class ObservationTable extends HashMap<String, GenotypeFeatureData> {
    String statField = null;
    File observationsFile = null;
    Logger logger = null;
    GenotypeFeatureData.Tail tail;

    public ObservationTable(final File observationsFile, final String statField, Logger logger, GenotypeFeatureData.Tail tail ) throws NoSuchFieldException, FileNotFoundException, IOException {
        this.observationsFile = observationsFile;
        this.statField = statField;
        this.logger = logger;
        this.tail = tail;
        readAllData(observationsFile, statField);
    }

    public Set<String> genotypes() {
        return keySet();
    }

    public void readAllData(final File observationsFile, final String statField) throws NoSuchFieldException, FileNotFoundException, IOException {
        ArrayList<String> header = TabularROD.readHeader(observationsFile);

        //logger.info(String.format("Starting to read table %s", observationsFile));

        for ( String line : new xReadLines(observationsFile) ) {
            TabularROD d = new TabularROD("ignoreMe", header);
            String[] parts = line.split("\\s+");
            if ( d.parseLine(header, parts) ) {
                if (! d.containsKey(statField)) {
                    throw new NoSuchFieldException(String.format("Could not find field %s in line %s", statField, line));
                }
                double stat = Double.valueOf(d.get(statField));

                final String genotype = d.get("genotype");
                GenotypeFeatureData gfd = getData(genotype);
                gfd.add(stat);
            }
        }

        for ( GenotypeFeatureData gfd : this.values() ) {
            gfd.finalizeData();
        }

        //logger.info(String.format("Read table %s, found %d genotypes/data pairs in %s field",
        //        observationsFile, size(), statField));
    }

    public GenotypeFeatureData getData(final String genotype) {
        if ( ! containsKey(genotype) ) {
            GenotypeFeatureData gfd = new GenotypeFeatureData(genotype, logger, tail);
            put(genotype, gfd);
        }

        return get(genotype);
    }

    public String toString() {
        return String.format("[ObservationTable: file=%s, field=%s, nGenotypes=%d]", observationsFile, statField, size());
    }
}

public abstract class RatioFilter implements VariantExclusionCriterion {
    protected double pvalueLimit = -1;
    protected File observationsFile = null;
    protected String statField = null; // "AlleleRatio";
    protected Logger logger = null; // Logger.getLogger(RatioFilter.class);
    protected ObservationTable dataTable = null;
    protected String name = null;
    protected GenotypeFeatureData.Tail tail = null;


    /**
     * A short-term hack to stop the systme from rejecting poorly covered sites, under the assumption
     * that the ratio test is poorly determined without at least MIN_COUNTS_TO_APPLY_TEST.  To be
     * replaced by the more robust sampling approach outlined below
     */
    final private static int MIN_COUNTS_TO_APPLY_TEST = 20;

    final static boolean integrateOverSamplingProbabilities = true;

    public RatioFilter(final String name, final String statField, Class myClass, GenotypeFeatureData.Tail tail ) {
        this.name = name;
        this.statField = statField;
        this.tail = tail;
        logger = Logger.getLogger(myClass);
    }

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            String[] argPieces = arguments.split(",");
            try {
                pvalueLimit = Double.valueOf(argPieces[0]);
                observationsFile = new File(argPieces[1]);

                dataTable = new ObservationTable(observationsFile, statField, logger, tail);

                logger.info(String.format("Initialized data for ratio filter %s: %s", name, dataTable));
                for ( String genotype : dataTable.genotypes() ) {
                    GenotypeFeatureData gfd = dataTable.get(genotype);
                    gfd.determineThresholds(pvalueLimit);
                    logger.info(gfd.toString());
                }
            } catch ( NumberFormatException e ) {
                throw new RuntimeException(String.format("Couldn't parse pValue limit from %s", argPieces[0]), e);
            } catch ( FileNotFoundException e ) {
                throw new RuntimeException("Couldn't open ObservationTable " + observationsFile, e);
            } catch ( NoSuchFieldException e ) {
                throw new RuntimeException("Couldn't parse ObservationTable " + observationsFile, e);
            } catch ( IOException e ) {
                throw new RuntimeException("Couldn't parse ObservationTable " + observationsFile,e );
            }
        }
    }

    protected abstract boolean applyToVariant(rodVariants variant);
    protected abstract Pair<Integer, Integer> scoreVariant(char ref, ReadBackedPileup pileup, rodVariants variant);

    public boolean exclude(char ref, LocusContext context, rodVariants variant) {
        boolean exclude = false;

        //
        // todo -- need to calculate a significance value for the chance that the
        // todo -- observed first and second counts are reliably not passing the threshold
        // todo -- basically, we need to sample with replacement from all first / second pairs
        // todo -- and only conclude that the variant fails the test if the probability of failing
        // todo -- across all samples is greater than some P value, like 0.05
        //
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        if ( applyToVariant(variant) ) {
            Pair<Integer, Integer> counts = scoreVariant(ref, pileup, variant);
            GenotypeFeatureData gfd = dataTable.get(variant.getBestGenotype());

            if (integrateOverSamplingProbabilities)
                exclude = integralExclude(gfd, counts);
            else
                exclude = pointEstimateExclude(gfd, counts);

            //
            // for printing only
            //
            int n = counts.first + counts.second;
            double value = counts.first / (1.0 * counts.first + counts.second);
            logger.info(String.format("%s: counts1=%d (%.2f), counts2=%d (%.2f), n=%d, value=%f, %s, exclude=%b, bases=%s",
                    name, counts.first, counts.first / (0.01 * n), counts.second, counts.second / (0.01 * n), n, 
                    value, gfd, exclude, pileup.getBases()));
        }

        return exclude;
    }

    private final static double SEARCH_INCREMENT = 0.01;
    private final static double integralPValueThreshold = 0.05;

    private boolean integralExclude(GenotypeFeatureData gfd, Pair<Integer, Integer> counts ) {
        double sumExclude = 0.0, sumP = 0.0;
        int n = counts.first + counts.second;
        for ( double r = 0.0; r <= 1.0; r += SEARCH_INCREMENT ) {
            double p = MathUtils.binomialProbability(counts.first, n, r);
            sumP += p;
            boolean exclude = ! gfd.passesThreshold(r);
            sumExclude += p * (exclude ? 1.0 : 0.0);

            //System.out.printf("integral: k=%d, n=%d, r=%f, p=%f, sumP = %f, exclude=%b | sum=%f, percentExcluded=%f%n",
            //        counts.first, n, r, p, sumP, exclude, sumExclude, sumExclude / sumP);
        }

        double percentExcluded = sumExclude / sumP;
        return 1 - percentExcluded <= integralPValueThreshold ;
    }

    private boolean pointEstimateExclude(GenotypeFeatureData gfd, Pair<Integer, Integer> counts ) {
        int n = counts.first + counts.second;
        if ( n < MIN_COUNTS_TO_APPLY_TEST ) return false;
        double value = counts.first / (1.0 * counts.first + counts.second);
        //double value = counts.first / (1.0 * Math.max(counts.second, 1.0));
        return ! gfd.passesThreshold(value);
    }
}
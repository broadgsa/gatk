package edu.mit.broad.picard.quality;

import edu.mit.broad.picard.util.Histogram;

import java.util.TreeMap;
import java.util.Map;
import java.util.SortedMap;

/**
 * <p>Holds all the information necessary to perform quality score calibration for a single
 * end/read for a lane or run of sequencing.  General usage is to construct an instance
 * an call {@link #addObservation(int, int, boolean)} repeatedly and when all input data
 * is consumed call {@link #computeCalibratedQualities()}.</p>
 *
 * <p>Once this is done then {@link #getCalibratedQualities()} can be called to get a matrix
 * of quality score calibrations by cycle and input quality.  However it is preferred to call
 * {@link #getCalibratedQuality(int, int)} which will attempt to infer the correct value in the
 * case that the input quality was not observed in the training data.</p>
 *
 * @author Tim Fennell
 */
public class QualityScoreMatrix {
    // Maps by cycle, histograms by quality
    private SortedMap<Integer, Histogram<Integer>> observations = new TreeMap<Integer, Histogram<Integer>>();
    private SortedMap<Integer, Histogram<Integer>> errors       = new TreeMap<Integer, Histogram<Integer>>();

    private int[][] calibratedQualities = null;

    /**
     * Adds an observation to the matrix.
     * @param cycle the cycle in the read (1-based)
     * @param quality the uncalibrated quality
     * @param error true if the base did not match the reference, false otherwise
     */
    public void addObservation(int cycle, int quality, boolean error) {
        Histogram<Integer> obs = this.observations.get(cycle);
        if (obs == null) {
            obs = new Histogram<Integer>();
            this.observations.put(cycle, obs);
        }
        obs.increment(quality);

        if (error) {
            Histogram<Integer> errs = this.errors.get(cycle);
            if (errs == null) {
                errs = new Histogram<Integer>();
                this.errors.put(cycle, errs);
            }
            errs.increment(quality);
        }
    }

    /**
     * Takes the input observations so far and builds a matrix of input cycle and
     * uncalibrated quality to calibrated quality value.
     */
    public void computeCalibratedQualities() {
        this.calibratedQualities = new int[this.observations.lastKey() + 1][];

        for (int cycle=1; cycle<this.calibratedQualities.length; ++cycle) {
            Histogram<Integer> obs = this.observations.get(cycle);
            Histogram<Integer> err = this.errors.get(cycle);

            this.calibratedQualities[cycle] = new int[obs.lastKey() + 1];

            for (Integer qual : obs.keySet()) {
                double o = obs.get(qual).getValue();
                Histogram<Integer>.Bin errBin = err.get(qual);
                double e = (errBin == null) ? 1 : errBin.getValue();

                this.calibratedQualities[cycle][qual] = computePhredScore(e, o);
            }
        }
    }

    /**
     * Returns the set of calibrated quality scores from the training data. The array is
     * indexed first by the cycle (1-based, index 0 is empty) and then by input quality
     * (again, the actualy quality, not shifted).
     *
     * @return an array of calibrated qualities for the read
     */
    public int[][] getCalibratedQualities() {
        return calibratedQualities;
    }

    /**
     * Accesses the calibrated quality for the given input cycle and quality. If the quality
     * is outside the range given in the training data then the upper or lower bound of
     * the calibrated qualities is used instead.
     *
     * @param cycle the input cycle (1-based)
     * @param quality the uncalibrated quality
     * @return the calibrated quality for the cycle and uncalibrated quality
     */
    public final int getCalibratedQuality(int cycle, int quality) {
        final int[] quals = this.calibratedQualities[cycle];

        // TODO: proper iterpolation where we don't have the right quality
        try {
            int retval = quals[quality];

            // If we didn't calibrate this quality value, search up and down for non-zero
            for (int i=quality; i>0 && retval == 0; --i) {
                if (quals[i] != 0)  retval = quals[i];
            }

            for (int i=quality; i<quals.length && retval == 0; ++i) {
                if (quals[i] != 0)  retval = quals[i];
            }

            return retval;
        }
        catch (IndexOutOfBoundsException ioobe) {
            // If we try to fetch a quality out of the calibrted range use either
            // 1 or max quality based on which side we were out of range on
            if (quality < 1) return 1;
            else return quals[quals.length - 1];
        }
    }

    /** Returns true if no observations were made, otherwise false. */
    public boolean isEmpty() {
        return this.observations.isEmpty();
    }

    /** Just does the simple phred scaling given the errors and observations. */
    private int computePhredScore(double errors, double observations) {
        return (int) Math.round(-10d * Math.log10(errors / observations));
    }


}

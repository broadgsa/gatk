package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;

/**
 * BasecallingStats is a utility class to aggregate and emit basecalling
 * stats (total bases seen and consistency between basecalling methods).
 *
 * @author Kiran Garimella
 */
public class BasecallingStats {
    private int basesConsistent = 0;
    private int basesTotal = 0;
    private int readsTotal = 0;

    /**
     * Constructor that does nothing.
     */
    public BasecallingStats() {}

    /**
     * Return the number of bases called identically by two different methods.
     *
     * @return the number of consistent bases.
     */
    public int getBasesConsistent() {
        return basesConsistent;
    }

    /**
     * Return the total number of bases seen.
     *
     * @return the total number of bases seen.
     */
    public int getBasesTotal() {
        return basesTotal;
    }

    /**
     * Return the total number of reads seen.
     *
     * @return the total number of reads seen.
     */
    public int getReadsTotal() {
        return readsTotal;
    }

    /**
     * Return the percent of bases called consistently by two different methods.
     *
     * @return the percent of bases called consistently
     */
    public double getPercentConsistent() {
        return 100.0*((double) getBasesConsistent())/((double) getBasesTotal());
    }

    /**
     * Updates the number of bases seen, the number of reads seen, and the number of consistent bases.
     *
     * @param rr   the raw Illumina read
     * @param fpr  the FourProb read
     */
    public void update(RawRead rr, FourProbRead fpr) {
        for (int cycle = 0; cycle < fpr.size(); cycle++) {
            int rawBaseIndex = BaseUtils.simpleBaseToBaseIndex((char) rr.getSequence()[cycle]);
            int fpBaseIndex = fpr.get(cycle).indexAtRank(0);

            if (rawBaseIndex >= 0 && fpBaseIndex >= 0) {
                basesTotal++;

                if (rawBaseIndex == fpBaseIndex) {
                    basesConsistent++;
                }
            }
        }

        readsTotal++;
    }

    /**
     * Returns basecalling stats info in a nicely formatted string
     *
     * @return nicely formatted string containing basecalling stats
     */
    public String toString() {
        return String.format("%% bases consistent: %d/%d (%2.2f%%)", getBasesConsistent(), getBasesTotal(), getPercentConsistent());
    }

    /**
     * Periodically print a line containing basecalling stats
     *
     * @param interval  the periodicity of the messages given in number of bases observed
     */
    public void notifyOnInterval(int interval) {
        if (getBasesTotal() > 0 && getBasesTotal() % interval == 0) {
            System.out.printf("%s\r", toString());
        }
    }

    /**
     * Immediately print a line containing basecalling stats
     */
    public void notifyNow() {
        System.out.printf("%s\n", toString());
    }
}

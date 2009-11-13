package org.broadinstitute.sting.utils.genotype;

/**
 * @author ebanks
 * Interface AlleleFrequencyBacked
 *
 * this interface indicates that the genotype is
 * backed up by allele frequency information.
 */
public interface AlleleFrequencyBacked {

    /**
     *
     * @return returns the best allele frequency for this genotype
     */
    public double getAlleleFrequency();

    /**
     *
     * @param   frequency    the allele frequency for this genotype
     */
    public void setAlleleFrequency(double frequency);

    /**
     *
     * @return returns the allele frequency for this genotype
     */
    public AlleleFrequencyRange getAlleleFrequencyRange();

    /**
     *
     * @param   range    the allele frequency range for this genotype
     */
    public void setAlleleFrequencyRange(AlleleFrequencyRange range);


    /**
     * A class representing a range of allele frequencies that make up a given
     * fraction of the total probability over all frequencies.
     */
    public class AlleleFrequencyRange {

        private double mLow, mHigh, mFraction;

        public AlleleFrequencyRange(double low, double high, double fraction) {
            mLow = low;
            mHigh = high;
            mFraction = fraction;
        }

        public double getLowEnd() { return mLow; }
        public double getHighEnd() { return mHigh; }
        public double getPercentageOfProbability() { return mFraction; }

        public String toString() {
            return String.format("%.2f-%.2f,%.0f%%", mLow, mHigh, 100*mFraction);
        }
    }
}
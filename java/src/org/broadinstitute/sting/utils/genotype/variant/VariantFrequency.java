package org.broadinstitute.sting.utils.genotype.variant;


/**
 * 
 * @author aaron 
 * 
 * Class VariantFrequency
 *
 * a class that represents the variant frequency, and could serve as a base for any other
 * variant frequency information (i.e. is it pop gen, from chip, etc).
 */
public class VariantFrequency {

    private double mFrequency;

    /**
     * create a variant frequency
     * @param frequency
     */
    public VariantFrequency(double frequency) {
        this.mFrequency = frequency;
    }

    public double getFrequency() {
        return mFrequency;
    }
}

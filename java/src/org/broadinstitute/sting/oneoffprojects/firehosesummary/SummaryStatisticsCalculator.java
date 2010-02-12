package org.broadinstitute.sting.oneoffprojects.firehosesummary;

/**
 * This object calculates the first three sample moments of a data set on-the-fly
 * @Author chartl
 * @Date Feb 12, 2010
 */
public class SummaryStatisticsCalculator {
    // todo -- can median / quantiles be estimated on-the-fly?
    // todo -- can we divine a metric for modality?
    private double mean; // mean of the samples fed to it
    private double var; // variance of the samples fed to it
    private double skew; // skew of the samples fed to it
    private int sampleSize; // number of samples given to calculator
    private String name; // name to associate calculator with, if any

    public SummaryStatisticsCalculator() {
        mean = 0;
        var = 0;
        skew = 0;
        sampleSize = 0;
        name = "SUMMARY";
    }

    public SummaryStatisticsCalculator(String name) {
        mean = 0;
        var = 0;
        skew = 0;
        sampleSize = 0;
        this.name = name;
    }

    public void update(double datum) {
        mean = (sampleSize*mean + datum)/(sampleSize + 1);
        var = (sampleSize*var + Math.pow( mean - datum , 2 ) )/( sampleSize + 1);
        // note: the variance is not re-calculated over all data points each time the mean changes
        // this means the convergence will be slower than x^-(1/2); but it will still converge with the mean
        // whole exome is ~3.3 million datapoints, which will be plenty for this to converge to within a very
        // small proportion of the true sample variance
        skew  =  ( skew*Math.pow(var,2/3)*sampleSize+ Math.pow( datum - mean, 3) ) / ( ( sampleSize + 1 )* Math.pow(var,2/3) );
        // again, new mean/variance estimates are not propagated over all previous data points, so convergence is a bit slower
        // but that's the price one pays for on-the-fly calculation
        sampleSize++;
    }
    
    public void update(int datum) {
        update( (double) datum);
    }

    public double getMean() {
        return mean;
    }

    public double getVar() {
        return var;
    }

    public double getSkew() {
        return skew;
    }

    public String getName() {
        return name;
    }

    public int getSampleSize() {
        return sampleSize;
    }
}

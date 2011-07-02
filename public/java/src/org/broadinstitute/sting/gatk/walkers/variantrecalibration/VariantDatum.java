package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDatum implements Comparable<VariantDatum> {

    public double[] annotations;
    public boolean[] isNull;
    public boolean isKnown;
    public double lod;
    public boolean atTruthSite;
    public boolean atTrainingSite;
    public boolean isTransition;
    public boolean isSNP;
    public boolean failingSTDThreshold;
    public double originalQual;
    public double prior;
    public int consensusCount;
    public int usedForTraining;
    public String contig;
    public int start;
    public int stop;
    public int worstAnnotation;
    public MultivariateGaussian assignment; // used in K-means implementation 

    public int compareTo( final VariantDatum other ) {
        return Double.compare(this.lod, other.lod);
    }
}

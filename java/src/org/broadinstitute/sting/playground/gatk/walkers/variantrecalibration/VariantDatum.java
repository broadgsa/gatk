package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantDatum implements Comparable<VariantDatum> {

    public double[] annotations;
    public boolean isKnown;
    public double lod;
    public double pVarGivenModel;
    public boolean atTruthSite;
    public boolean isTransition;
    public MultivariateGaussian assignment; // used in K-means implementation 

    public int compareTo( final VariantDatum other ) {
        return Double.compare(this.lod, other.lod);
    }
}

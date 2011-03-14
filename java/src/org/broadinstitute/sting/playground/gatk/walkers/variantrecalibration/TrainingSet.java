package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/12/11
 */

public class TrainingSet implements Comparable<TrainingSet> {

    public String name;
    public boolean isKnown;
    public boolean isTraining;
    public boolean isTruth;
    public double prior;

    public TrainingSet( final String name, final boolean isKnown, final boolean isTraining, final boolean isTruth, final double prior ) {
        this.name = name;
        this.isKnown = isKnown;
        this.isTraining = isTraining;
        this.isTruth = isTruth;
        this.prior = prior;
    }

    public int compareTo( final TrainingSet other ) {
        return Double.compare(this.prior, other.prior);
    }
}

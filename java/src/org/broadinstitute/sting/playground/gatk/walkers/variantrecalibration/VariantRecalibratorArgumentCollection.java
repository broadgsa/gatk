package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.broadinstitute.sting.commandline.Argument;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 4, 2011
 */

public class VariantRecalibratorArgumentCollection {

    @Argument(fullName="maxGaussians", shortName="mG", doc="The maximum number of Gaussians to try during variational Bayes algorithm", required=false)
    public int MAX_GAUSSIANS = 32;
    @Argument(fullName="maxIterations", shortName="mI", doc="The maximum number of VBEM iterations to be performed in variational Bayes algorithm. Procedure will normally end when convergence is detected.", required=false)
    public int MAX_ITERATIONS = 100;
    @Argument(fullName="numKMeans", shortName="nKM", doc="The number of k-means iterations to perform in order to initialize the means of the Gaussians in the Gaussian mixture model.", required=false)
    public int NUM_KMEANS_ITERATIONS = 10;
    @Argument(fullName="stdThreshold", shortName="std", doc="If a variant has annotations more than -std standard deviations away from mean then don't use it for building the Gaussian mixture model.", required=false)
    public double STD_THRESHOLD = 4.5;
    @Argument(fullName="qualThreshold", shortName="qual", doc="If a known variant has raw QUAL value less than -qual then don't use it for building the Gaussian mixture model.", required=false)
    public double QUAL_THRESHOLD = 80.0;
    @Argument(fullName="shrinkage", shortName="shrinkage", doc="The shrinkage parameter in variational Bayes algorithm.", required=false)
    public double SHRINKAGE = 0.0001;
    @Argument(fullName="dirichlet", shortName="dirichlet", doc="The dirichlet parameter in variational Bayes algorithm.", required=false)
    public double DIRICHLET_PARAMETER = 0.0001;
    @Argument(fullName="percentBadVariants", shortName="percentBad", doc="What percentage of the worst scoring variants to use when building the Gaussian mixture model of bad variants. 0.07 means bottom 7 percent.", required=false)
    public double PERCENT_BAD_VARIANTS = 0.07;
}

package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * This class emcompasses all the basic information about a genotype.  It is immutable.
 *
 * @author ebanks
 */
public class CalledGenotype extends MutableGenotype {

    public static final String LIKELIHOODS_ATTRIBUTE_KEY = "Likelihoods";
    public static final String POSTERIORS_ATTRIBUTE_KEY = "Posteriors";
    public static final String READBACKEDPILEUP_ATTRIBUTE_KEY = "ReadBackedPileup";

    public CalledGenotype(String sampleName, List<Allele> alleles, double negLog10PError, Set<String> filters, Map<String, ?> attributes, boolean genotypesArePhased) {
        super(sampleName, alleles, negLog10PError, filters, attributes, genotypesArePhased);
    }

    public CalledGenotype(String sampleName, List<Allele> alleles, double negLog10PError) {
        super(sampleName, alleles, negLog10PError);
    }

    public CalledGenotype(String sampleName, List<Allele> alleles) {
        super(sampleName, alleles, NO_NEG_LOG_10PERROR, null, null, false);
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // routines to modify useful attribute fields
    //
    // ---------------------------------------------------------------------------------------------------------

    public void setLikelihoods(double likelihoods[])  { putAttribute(LIKELIHOODS_ATTRIBUTE_KEY, likelihoods); }
    public void setPosteriors(double posteriors[])  { putAttribute(POSTERIORS_ATTRIBUTE_KEY, posteriors); }
    public void setReadBackedPileup(ReadBackedPileup pileup)  { putAttribute(READBACKEDPILEUP_ATTRIBUTE_KEY, pileup); }

    public double[] getLikelihoods()  { return (double[])getAttribute(LIKELIHOODS_ATTRIBUTE_KEY); }
    public double[] getPosteriors()  { return (double[])getAttribute(POSTERIORS_ATTRIBUTE_KEY); }
    public ReadBackedPileup getReadBackedPileup()  { return (ReadBackedPileup)getAttribute(READBACKEDPILEUP_ATTRIBUTE_KEY); }

}
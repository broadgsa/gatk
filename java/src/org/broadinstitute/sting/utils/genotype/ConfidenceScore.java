package org.broadinstitute.sting.utils.genotype;


/**
 * @author aaron
 *         <p/>
 *         Class ConfidenceScore
 *         <p/>
 *         this class represents the confidence in a genotype, and the method we used to obtain it
 */
public class ConfidenceScore implements Comparable<ConfidenceScore> {

    private static final Double EPSILON = 1.0e-15;

    public enum SCORE_METHOD {
        BEST_NEXT, BEST_REF, OTHER;
    }

    private Double mScore;
    private SCORE_METHOD mMethod;

    public ConfidenceScore(double score, SCORE_METHOD method) {
        this.mScore = score;
        this.mMethod = method;
    }

    /**
     * generate a confidence score, given the two likelihoods, and the method used
     *
     * @param likelihoodOne the first likelihood
     * @param likelihoodTwo the second likelihood
     * @param method        the method used to determine the likelihood
     */
    public ConfidenceScore(double likelihoodOne, double likelihoodTwo, SCORE_METHOD method) {
        this.mScore = likelihoodOne / likelihoodTwo;
        this.mMethod = method;
    }

    /**
     * compare this ConfidenceScore to another, throwing an exception if they're not the same scoring method
     * @param o the other confidence score if 
     * @return 0 if equal
     */
    @Override
    public int compareTo(ConfidenceScore o) {
        if (o.mMethod != this.mMethod) {
            throw new UnsupportedOperationException("Attemped to compare Confidence scores with different methods");
        }
        double diff = mScore - o.mScore;
        if (Math.abs(diff) < (EPSILON * Math.abs(mScore)))
            return 0;
        else if (diff < 0)
            return -1;
        else
            return 1;
    }

    public Double getScore() {
        return mScore;
    }

    public SCORE_METHOD getMethod() {
        return mMethod;
    }
}

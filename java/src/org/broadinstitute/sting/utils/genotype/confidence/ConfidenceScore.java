package org.broadinstitute.sting.utils.genotype.confidence;

/**
 * @author aaron
 *         <p/>
 *         Class ConfidenceScore
 *         <p/>
 *         this class represents the confidence in a genotype, and the method we used to obtain it
 */
public abstract class ConfidenceScore implements Comparable<ConfidenceScore> {

    // the general error of a floating point value
    private static final Double EPSILON = 1.0e-15;

    public enum SCORE_METHOD {
        LOD_SCORE, UNKNOWN;
    }

    private Double mScore;

    public ConfidenceScore(double score) {
        this.mScore = score;
    }

    /**
     * compare this ConfidenceScore to another, throwing an exception if they're not the same scoring method
     * @param o the other confidence score if 
     * @return 0 if equal
     */
    @Override
    public int compareTo(ConfidenceScore o) {
        if (o.getConfidenceType() != this.getConfidenceType()) {
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

    /**
     * get the score
     * @return a double representing the genotype score
     */
    public Double getScore() {
        return mScore;
    }

    /**
     * return the confidence method we're employing, UNKNOWN is an option
     * @return the method of confidence we represent
     */
    public abstract SCORE_METHOD getConfidenceType();

    /**
     * get the confidence score, normalized to the range of [0-1]
     * @return a confidence score
     */
    public abstract double normalizedConfidenceScore();
}

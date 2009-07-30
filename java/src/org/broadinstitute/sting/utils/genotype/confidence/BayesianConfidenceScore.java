package org.broadinstitute.sting.utils.genotype.confidence;


/**
 * 
 * @author aaron 
 * 
 * Class LikelihoodConfidenceScore
 *
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class BayesianConfidenceScore extends ConfidenceScore {
    public BayesianConfidenceScore(double score) {
        super(score);
    }

    /**
     * return the confidence method we're employing, UNKNOWN is an option
     *
     * @return the method of confidence we represent
     */
    @Override
    public SCORE_METHOD getConfidenceType() {
        return SCORE_METHOD.LOD_SCORE;
    }

    /**
     * get the confidence score, normalized to the range of [0-1]
     *
     * @return a confidence score
     */
    @Override
    public double normalizedConfidenceScore() {
        return Math.pow(10,this.getScore());
    }
}

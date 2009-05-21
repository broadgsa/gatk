package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;

import java.io.File;
import java.util.ArrayList;

/**
 * BasecallingReadModel represents the statistical models for
 * all bases in all cycles.  It allows for easy training via
 * the addTrainingPoint() method, and for the computation of
 * the 4x4 likelihood matrix or the 1x4 probability vector
 * (with contextual components marginalized out of the
 * likelihood matrix).
 *
 * @author Kiran Garimella
 */
public class BasecallingReadModel {
    private BasecallingBaseModel[] basemodels = null;
    private boolean correctForContext = true;

    /**
     * Constructs a BasecallingReadModel with space for a given read length.
     * 
     * @param readLength  the length of the reads to which this model will apply.
     */
    public BasecallingReadModel(int readLength) {
        initialize(readLength);
    }

    /**
     * Constructs a BasecallingReadModel and trains it using the specified training data.
     *
     * @param trainingData  a set of RawReads from which the model will be trained.
     */
    public BasecallingReadModel(ArrayList<RawRead> trainingData) {
        initialize(trainingData.get(0).getReadLength());

        train(trainingData);
    }

    /**
     * Initialize the model and set default parameters for each cycle appropriately.
     *
     * @param readLength  the length of the reads to which this model will apply.
     */
    public void initialize(int readLength) {
        basemodels = new BasecallingBaseModel[readLength];

        for (int cycle = 0; cycle < readLength; cycle++) {
            basemodels[cycle] = new BasecallingBaseModel(cycle != 0 && correctForContext);
        }
    }

    /**
     * Train the model using the specified training data.
     *
     * @param trainingData  a set of RawReads from which the model will be trained.
     */
    public void train(ArrayList<RawRead> trainingData) {
        for ( RawRead read : trainingData ) {
            addMeanPoints(read);
        }

        for ( RawRead read : trainingData ) {
            addCovariancePoints(read);
        }
    }

    /**
     * Add a training point for the mean intensity values per base and per cycle.
     *
     * @param cycle          the cycle number (0-based)
     * @param probMatrix     the probability matrix for the base
     * @param fourintensity  the four raw intensities for the base
     */
    public void addMeanPoint(int cycle, double[][] probMatrix, double[] fourintensity) {
        basemodels[cycle].addMeanPoint(probMatrix, fourintensity);
    }

    /**
     * Add a training point for the mean intensity values per base in all cycles.
     *
     * @param read  the raw read
     */
    public void addMeanPoints(RawRead read) {
        byte[] seqs = read.getSequence();
        byte[] quals = read.getQuals();
        short[][] ints = read.getIntensities();

        for (int cycle = 0; cycle < seqs.length; cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : seqs[cycle - 1]);
            char baseCur = (char) seqs[cycle];
            double probCur = QualityUtils.qualToProb(quals[cycle]);

            double[][] probMatrix = getBaseProbabilityMatrix(cycle, basePrev, baseCur, probCur);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                fourIntensity[channel] = (double) ints[cycle][channel];
            }

            basemodels[cycle].addMeanPoint(probMatrix, fourIntensity);
        }
    }

    /**
     * Add a training point for the intensity covariance matrix per base and per cycle.
     *
     * @param cycle          the cycle number (0-based)
     * @param probMatrix     the probability matrix for the base
     * @param fourintensity  the four raw intensities for the base
     */
    public void addCovariancePoint(int cycle, double[][] probMatrix, double[] fourintensity) {
        basemodels[cycle].addCovariancePoint(probMatrix, fourintensity);
    }

    /**
     * Add a training point for the intensity covariance matrix per base in all cycles.
     *
     * @param read  the raw read
     */
    public void addCovariancePoints(RawRead read) {
        byte[] seqs = read.getSequence();
        byte[] quals = read.getQuals();
        short[][] ints = read.getIntensities();

        for (int cycle = 0; cycle < seqs.length; cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : seqs[cycle - 1]);
            char baseCur = (char) seqs[cycle];
            double probCur = QualityUtils.qualToProb(quals[cycle]);

            double[][] probMatrix = getBaseProbabilityMatrix(cycle, basePrev, baseCur, probCur);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                fourIntensity[channel] = (double) ints[cycle][channel];
            }

            basemodels[cycle].addCovariancePoint(probMatrix, fourIntensity);
        }
    }

    /**
     * Compute the likelihoods that a given set of intensities yields each possible base.
     *
     * @param cycle          the cycle number (0-based)
     * @param fourintensity  the four raw intensities for the base
     * @return               the matrix of likelihoods
     */
    public double[][] computeLikelihoods(int cycle, double[] fourintensity) {
        return basemodels[cycle].computeLikelihoods(cycle, fourintensity);
    }

    /**
     * Compute the probabilities that a given set of intensities yields each possible base.
     *
     * @param cycle          the cycle number (0-based)
     * @param basePrev       the previous base
     * @param qualPrev       the previous base's quality score
     * @param fourintensity  the four raw intensities for the base
     * @return the probability distribution over the four base possibilities
     */
    public FourProb computeProbabilities(int cycle, char basePrev, byte qualPrev, double[] fourintensity) {
        double[][] likes = computeLikelihoods(cycle, fourintensity);

        double total = 0;

        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                double prior = 1.0;
                if (correctForContext) {
                    double prob = QualityUtils.qualToProb(qualPrev);
                    if (basePrevIndex == BaseUtils.simpleBaseToBaseIndex(basePrev)) {
                        prior = prob;
                    } else {
                        prior = (1.0 - prob)/((double) (4*likes.length - 1));
                    }
                }
                likes[basePrevIndex][baseCurIndex] = prior*likes[basePrevIndex][baseCurIndex];
                total += likes[basePrevIndex][baseCurIndex];
            }
        }

        for (int basePrevIndex = 0; basePrevIndex < likes.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < 4; baseCurIndex++) {
                likes[basePrevIndex][baseCurIndex] /= total;
            }
        }

        return new FourProb(likes);
    }

    /**
     * Call the bases in the given RawRead.
     *
     * @param read  the RawRead
     * @return  the basecalled read
     */
    public FourProbRead call(RawRead read) {
        FourProbRead fpr = new FourProbRead(read.getReadLength());
        
        for (int cycle = 0; cycle < read.getReadLength(); cycle++) {
            char basePrev = (char) ((cycle == 0) ? '.' : read.getSequence()[cycle - 1]);
            byte qualPrev = ((cycle == 0) ? 0 : read.getQuals()[cycle - 1]);

            double[] fourIntensity = new double[4];
            for (int channel = 0; channel < 4; channel++) {
                fourIntensity[channel] = (double) read.getIntensities()[cycle][channel];
            }

            fpr.add(cycle, computeProbabilities(cycle, basePrev, qualPrev, fourIntensity));

        }

        return fpr;
    }

    /**
     * Return the probability matrix given the previous cycle's base, the current cycle's base, and the current base's probability.
     *
     * @param cycle     the cycle number (0-based)
     * @param basePrev  the previous base
     * @param baseCur   the current base
     * @param probCur   the probability of the current base
     * @return the probability matrix of the base
     */
    public double[][] getBaseProbabilityMatrix(int cycle, char basePrev, char baseCur, double probCur) {
        double[][] dist = new double[(correctForContext && cycle > 0) ? 4 : 1][4];

        int actualBasePrevIndex = (correctForContext && cycle > 0) ? BaseUtils.simpleBaseToBaseIndex(basePrev) : 0;
        int actualBaseCurIndex = BaseUtils.simpleBaseToBaseIndex(baseCur);

        double residualTheories = (double) (dist.length*dist[0].length - 1);

        for (int basePrevIndex = 0; basePrevIndex < dist.length; basePrevIndex++) {
            for (int baseCurIndex = 0; baseCurIndex < dist[basePrevIndex].length; baseCurIndex++) {
                dist[basePrevIndex][baseCurIndex] = (basePrevIndex == actualBasePrevIndex && baseCurIndex == actualBaseCurIndex) ? probCur : ((1.0 - probCur)/residualTheories);
            }
        }

        return dist;
    }

    /**
     * Write model parameters to disk.
     *
     * @param dir  the directory in which model parameters should be stored.
     */
    public void write(File dir) {
        for (int cycle = 0; cycle < basemodels.length; cycle++) {
            File outparam = new File(dir.getPath() + "/param." + cycle + ".r");
            basemodels[cycle].write(outparam);
        }
    }
}

package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.containers.BoundedScoringSet;
import org.broadinstitute.sting.utils.StingException;

import java.util.ArrayList;
import java.util.Arrays;

public class SecondaryBaseAnnotator {
    private static final int TRAINING_LIMIT = 10000;
    private boolean trained;
    private final BoundedScoringSet<RawRead> trainingAggregator;
    public BasecallingReadModel model;

    public SecondaryBaseAnnotator() {
        trained = false;
        trainingAggregator = new BoundedScoringSet<RawRead>(TRAINING_LIMIT);
    }

    public void addTrainingRead(RawRead rawRead) { trainingAggregator.add(rawRead); }

    public boolean haveEnoughTrainingReads() { return false; }

    public void doneTraining() {
        ArrayList<RawRead> trainingData = new ArrayList<RawRead>(trainingAggregator.size());

        trainingData.addAll(Arrays.asList(trainingAggregator.toArray(new RawRead[0])));

        model = new BasecallingReadModel(trainingData);

        trained = true;
    }

    public FourProbRead getFourProbRead(RawRead rawRead) {
        return model.call(rawRead);
    }

    public byte[] getSqTagValue(RawRead rawRead) {
        if (!trained) {
            throw new StingException("Model must be trained via addTrainingRead() before getSqTagValue() can be called");
        }

        FourProbRead fpr = model.call(rawRead);

        return fpr.getSQTag(rawRead);
    }

    private void train() {}

    private byte[] getSQTag(FourProbRead fourProbRead, RawRead rawRead) { return null; }

    private static boolean isGoodTrainingRead(RawRead rawRead) { return false; }

    private static double getAverageQualityScore(RawRead rawRead) { return 0.0; }
}

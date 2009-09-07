package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.HashMap;

public interface ConcordanceType {

    public void initialize(HashMap<String,String> args);
    public void computeConcordance(RefMetaDataTracker tracker, ReferenceContext ref);
    public void cleanup();

}
package org.broadinstitute.sting.playground.gatk.walkers.diagnostics.newvarianteval;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;

import java.util.List;

public interface VariantEvaluation {
    public void update(RodVCF variant, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context);

    //public Integer map(

    public String getName();

    public List<String> getResult();
}

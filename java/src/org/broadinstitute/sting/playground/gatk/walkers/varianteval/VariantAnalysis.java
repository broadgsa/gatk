package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;

public interface VariantAnalysis {
    public String getName();
    public PrintStream getPrintStream();
    public List<String> getParams();
    public void initialize(VariantEvalWalker master, PrintStream out);
    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context);
    public List<String> done();
}

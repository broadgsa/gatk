package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;

public interface VariantAnalysis {
    public String getName();
    public PrintStream getSummaryPrintStream();
    public PrintStream getCallPrintStream();
    public List<String> getParams();
    public void initialize(VariantEvalWalker master, PrintStream out, String filename);
    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context);
    public void finalize(long nSites);
    public List<String> done();
}

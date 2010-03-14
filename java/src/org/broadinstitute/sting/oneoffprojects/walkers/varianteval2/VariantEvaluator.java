package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;

import java.util.List;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
abstract class VariantEvaluator {
//    protected boolean accumulateInterestingSites = false, printInterestingSites = false;
//    protected String interestingSitePrefix = null;
    protected boolean processedASite = false;
//    protected List<VariantContext> interestingSites = new ArrayList<VariantContext>();

    public abstract String getName();

    // do we want to keep track of things that are interesting
//    public void accumulateInterestingSites(boolean enable)              { accumulateInterestingSites = enable; }
//    public void printInterestingSites(String prefix)                    { printInterestingSites = true; interestingSitePrefix = prefix; }
//    public boolean isAccumulatingInterestingSites()                     { return accumulateInterestingSites; }
//    public List<VariantContext> getInterestingSites()                   { return interestingSites; }

//    protected void addInterestingSite(String why, VariantContext vc) {
//        if ( accumulateInterestingSites )
//            interestingSites.add(vc);
//        if ( printInterestingSites )
//            System.out.printf("%40s %s%n", interestingSitePrefix, why);
//    }

    public abstract boolean enabled();
    //public boolean processedAnySites()                                  { return processedASite; }
    //protected void markSiteAsProcessed()                                { processedASite = true; }

    // Should return the number of VariantContexts expected as inputs to update.  Can be 1 or 2
    public abstract int getComparisonOrder();

    // called at all sites, regardless of eval context itself; useful for counting processed bases
    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { }

    public String update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return null;
    }

    public String update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return null;
    }

    public void finalize() {}

    public abstract String toString();

    // making it a table
    public abstract List<String> getTableHeader();
    public abstract List<List<String>> getTableRows();

    //
    // useful common utility routines
    //
    protected double rate(long n, long d) {
        return n / (1.0 * Math.max(d, 1));
    }

    protected long inverseRate(long n, long d) {
        return n == 0 ? 0 : d / Math.max(n, 1);
    }

    protected double ratio(long num, long denom) {
        return ((double)num) / (Math.max(denom, 1));
    }

}

package org.broadinstitute.sting.oneoffprojects.walkers.varianteval2;

import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.List;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jan 29, 2010
 * Time: 3:38:02 PM
 * To change this template use File | Settings | File Templates.
 */
abstract class VariantEvaluator {
    protected boolean accumulateInterestingSites = false, printInterestingSites = false;
    protected String interestingSitePrefix = null;
    protected boolean processedASite = false;
    protected List<VariantContext> interestingSites = new ArrayList<VariantContext>();

    public abstract String getName();

    // do we want to keep track of things that are interesting
    public void accumulateInterestingSites(boolean enable)              { accumulateInterestingSites = enable; }
    public void printInterestingSites(String prefix)                    { printInterestingSites = true; interestingSitePrefix = prefix; }
    public boolean isAccumulatingInterestingSites()                     { return accumulateInterestingSites; }
    public List<VariantContext> getInterestingSites()                   { return interestingSites; }
    protected void addInterestingSite(String why, VariantContext vc) {
        if ( accumulateInterestingSites )
            interestingSites.add(vc);
        if ( printInterestingSites )
            System.out.printf("%40s %s%n", interestingSitePrefix, why);
    }

    public boolean processedAnySites()                                  { return processedASite; }
    protected void markSiteAsProcessed()                                { processedASite = true; }

    /** Should return the number of VariantContexts expected as inputs to update.  Can be 1 or 2 */
    public abstract int getComparisonOrder();

    /** called at all sites, regardless of eval context itself; useful for counting processed bases */
    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { }
    public void update1(VariantContext vc1, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { }
    public void update2(VariantContext vc1, VariantContext vc2, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) { }

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

package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

public class VariantDBCoverage extends BasicVariantAnalysis {
    private String dbName;
    private int nDBObs = 0;
    private int nEvalObs = 0;
    private int nOverlapping = 0;

    public VariantDBCoverage(final String name) {
        super("db_coverage");
        dbName = name;
    }

    public void inc(boolean inDB, boolean inEval) {
        if (inDB) nDBObs++;
        if (inEval) nEvalObs++;
        if (inDB && inEval) nOverlapping++;
    }

    public int nDBSites()           { return nDBObs; }
    public int nEvalSites()         { return nEvalObs; }
    public int nOverlappingSites()  { return nOverlapping; }
    public int nNovelSites()        { return Math.abs(nEvalSites() - nOverlappingSites()); }

    /**
     * What fraction of the evaluated site variants were also found in the db?
     *
     * @return
     */
    public double fractionEvalSitesCoveredByDB() {
        return nOverlappingSites() / (1.0 * nEvalSites());
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        // There are four cases here:
        AllelicVariant dbsnp = (AllelicVariant)tracker.lookup(dbName, null);
        inc(dbsnp != null, eval != null);
        return dbsnp == null && eval != null ? "Novel " + eval : null;
    }

    /**
     * What fraction of the DB sites were discovered in the evalution calls?
     *
     * @return
     */
    public double fractionDBSitesDiscoveredInEval() {
        return nOverlappingSites() / (1.0 * nDBSites());
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("%d\t%d\t%d\t%.2f\t%.2f", nDBSites(), nEvalSites(), nOverlappingSites(), fractionEvalSitesCoveredByDB(), fractionDBSitesDiscoveredInEval()));
        s.add(String.format("name                 %s", dbName));
        s.add(String.format("n_db_sites           %d", nDBSites()));
        s.add(String.format("n_eval_sites         %d", nEvalSites()));
        s.add(String.format("n_overlapping_sites  %d", nOverlappingSites()));
        s.add(String.format("n_novel_sites        %d", nNovelSites()));
        s.add(String.format("per_eval_sites_in_db %.2f", 100*fractionEvalSitesCoveredByDB()));
        s.add(String.format("per_db_sites_in_eval %.2f", 100*fractionDBSitesDiscoveredInEval()));
        return s;
    }
}
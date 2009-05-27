package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

public class VariantDBCoverage {
    private String dbName;
    private int nDBObs = 0;
    private int nEvalObs = 0;
    private int nOverlapping = 0;

    public VariantDBCoverage(final String name) {
        dbName = name;
    }

    public void inc(boolean inDB, boolean inEval) {
        if (inDB) nDBObs++;
        if (inEval) nEvalObs++;
        if (inDB && inEval) nOverlapping++;
    }

    public int nDBSites() {
        return nDBObs;
    }

    public int nEvalSites() {
        return nEvalObs;
    }

    public int nOverlappingSites() {
        return nOverlapping;
    }

    /**
     * What fraction of the evaluated site variants were also found in the db?
     *
     * @return
     */
    public double fractionEvalSitesCoveredByDB() {
        return nOverlappingSites() / (1.0 * nEvalSites());
    }

    /**
     * What fraction of the DB sites were discovered in the evalution calls?
     *
     * @return
     */
    public double fractionDBSitesDiscoveredInEval() {
        return nOverlappingSites() / (1.0 * nDBSites());
    }

    public String toSingleLineString(final String prefix) {
        return String.format("%s %d\t%d\t%d\t%.2f\t%.2f", prefix, nDBSites(), nEvalSites(), nOverlappingSites(), fractionEvalSitesCoveredByDB(), fractionDBSitesDiscoveredInEval());
    }

    public String toMultiLineString(final String prefix) {
        StringBuilder s = new StringBuilder();
        s.append(String.format("%s name                 %s%n", prefix, dbName));
        s.append(String.format("%s n_db_sites           %d%n", prefix, nDBSites()));
        s.append(String.format("%s n_eval_sites         %d%n", prefix, nEvalSites()));
        s.append(String.format("%s n_overlapping_sites  %d%n", prefix, nOverlappingSites()));
        s.append(String.format("%s per_eval_sites_in_db %.2f%n", prefix, 100*fractionEvalSitesCoveredByDB()));
        s.append(String.format("%s per_db_sites_in_eval %.2f%n", prefix, 100*fractionDBSitesDiscoveredInEval()));
        return s.toString();
    }
}
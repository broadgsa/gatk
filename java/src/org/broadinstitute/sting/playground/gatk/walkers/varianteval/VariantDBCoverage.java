package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.util.List;
import java.util.ArrayList;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */
public class VariantDBCoverage extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    private String dbName;
    private int nDBObs = 0;
    private int nEvalObs = 0;
    private int nOverlapping = 0;
    private int nConcordant = 0;

    public VariantDBCoverage(final String name) {
        super("db_coverage");
        dbName = name;
    }

    public void inc(AllelicVariant dbSNP, AllelicVariant eval) {
        boolean inDB = dbSNP != null;
        boolean inEval = eval != null;

        if (inDB) nDBObs++;
        if (inEval) nEvalObs++;
        if (inDB && inEval) {
            nOverlapping++;

            if ( ! discordantP(dbSNP, eval) )
                nConcordant++;
        }
    }

    public int nDBSites()           { return nDBObs; }
    public int nEvalSites()         { return nEvalObs; }
    public int nOverlappingSites()  { return nOverlapping; }
    public int nConcordant()        { return nConcordant; }
    public int nNovelSites()        { return Math.abs(nEvalSites() - nOverlappingSites()); }

    public boolean discordantP(AllelicVariant dbSNP, AllelicVariant eval) {
        if (dbSNP != null && dbSNP.isSNP() && eval != null ) {
            boolean concordant = dbSNP.getAltSnpFWD() == eval.getAltSnpFWD() || dbSNP.getRefSnpFWD() == eval.getAltSnpFWD();

            //System.out.printf("dbSNP=%s | %c, eval=%s | %c, concordant=%b %s %s%n",
            //        dbSNP.getGenotype().get(0), dbSNP.getAltSnpFWD(),
            //        eval.getGenotype().get(0), eval.getAltSnpFWD(),
            //        concordant,
            //        dbSNP, eval);

            return ! concordant;
        } else {
            return false;
        }
    }


    /**
     * What fraction of the evaluated site variants were also found in the db?
     *
     * @return
     */
    public double fractionEvalSitesCoveredByDB() {
        return nOverlappingSites() / (1.0 * nEvalSites());
    }

    public double concordanceRate() {
        return nConcordant() / (1.0 * nOverlappingSites());
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        // There are four cases here:
        AllelicVariant dbsnp = (AllelicVariant)tracker.lookup(dbName, null);
        boolean isSNP = dbsnp != null && dbsnp.isSNP();
        inc(isSNP ? dbsnp : null, eval);
        return ! isSNP && eval != null ? "Novel      " + eval : (discordantP(dbsnp, eval) ? (String.format("Discordant DBSNP=%s %s", dbsnp.getGenotype().get(0), eval)) : null);
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
        s.add(String.format("n_concordant         %d", nConcordant()));
        s.add(String.format("n_novel_sites        %d", nNovelSites()));
        s.add(String.format("per_eval_sites_in_db %.2f", 100*fractionEvalSitesCoveredByDB()));
        s.add(String.format("concordance_rate     %.2f", 100*concordanceRate()));
        s.add(String.format("per_db_sites_in_eval %.2f", 100*fractionDBSitesDiscoveredInEval()));
        return s;
    }
}
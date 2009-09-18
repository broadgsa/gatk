package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.ArrayList;
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
public class VariantDBCoverage extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    private String dbName;
    private int nDBSNPs = 0;
    private int nDBIndels = 0;
    private int nEvalObs = 0;
    private int nOverlapping = 0;
    private int nConcordant = 0;
    private int nSNPsCalledAtIndels = 0;

    public VariantDBCoverage(final String name) {
        super("db_coverage");
        dbName = name;
    }

    public void inc(Variation dbSNP, Variation eval) {
        boolean inDB = dbSNP != null;
        boolean inEval = eval != null;

        if (inDB) {
            if (dbSNP.isSNP()) nDBSNPs++;
            if (dbSNP.isIndel()) nDBIndels++;
        }

        if (inEval) nEvalObs++;
        if (inDB && inEval) {
            if (dbSNP.isSNP()) { // changes the calculation a bit
                nOverlapping++;

                if (!discordantP(dbSNP, eval))
                    nConcordant++;
            }

            if (dbSNP.isIndel() && eval.isSNP())
                nSNPsCalledAtIndels++;
        }
    }

    public int callCount = 0;

    public int nDBSNPs() {
        return nDBSNPs;
    }

    public int nDBIndels() {
        return nDBIndels;
    }

    public int nEvalSites() {
        return nEvalObs;
    }

    public int nOverlappingSites() {
        return nOverlapping;
    }

    public int nConcordant() {
        return nConcordant;
    }

    public int nNovelSites() {
        return Math.abs(nEvalSites() - nOverlappingSites());
    }

    public int nSNPsAtIndels() {
        return nSNPsCalledAtIndels;
    }
    public boolean discordantP(Variation dbSNP, Variation eval) {
        if (eval != null) {
            char alt = (eval.isSNP()) ? eval.getAlternativeBaseForSNP() : eval.getReference().charAt(0);
            if (dbSNP != null && dbSNP.isSNP())
                return !dbSNP.getAlternateBase().contains(String.valueOf(alt));
        }
        return false;
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

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {

        // There are four cases here:
        Variation dbsnp = (Variation) tracker.lookup(dbName, null);
        boolean isSNP = dbsnp != null && dbsnp.isSNP();
        inc(dbsnp, eval);

        if (dbsnp != null && eval != null) {

            if (dbsnp.isSNP() && eval.isSNP() && discordantP(dbsnp, eval)) {
                return String.format("Discordant [DBSNP %s] [EVAL %s]", dbsnp, eval);
            } else if (dbsnp.isIndel() && eval.isSNP()) {
                return String.format("SNP-at-indel DBSNP=%s %s", dbsnp.getAlternateBase(), eval);
            } else {
                return null;
            }
        } else {
            return null;
        }
    }

    /**
     * What fraction of the DB sites were discovered in the evalution calls?
     *
     * @return
     */
    public double fractionDBSitesDiscoveredInEval() {
        return nOverlappingSites() / (1.0 * nDBSNPs());
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        //s.add(String.format("%d\t%d\t%d\t%.2f\t%.2f", nDBSNPs(), nEvalSites(), nOverlappingSites(), fractionEvalSitesCoveredByDB(), fractionDBSitesDiscoveredInEval()));
        s.add(String.format("name                     %s", dbName));

        s.add(String.format("n_db_sites               %d", nDBSNPs() + nDBIndels()));
        s.add(String.format("n_db_snps                %d", nDBSNPs()));
        s.add(String.format("n_db_indels              %d", nDBIndels()));
        s.add(String.format("n_eval_sites             %d", nEvalSites()));
        s.add(String.format("n_overlapping_sites      %d", nOverlappingSites()));
        s.add(String.format("n_concordant             %d", nConcordant()));
        s.add(String.format("n_novel_sites            %d", nNovelSites()));

        s.add(String.format("percent_eval_sites_in_db %.2f", 100 * fractionEvalSitesCoveredByDB()));
        s.add(String.format("concordance_rate         %.2f", 100 * concordanceRate()));

        s.add(String.format("percent_db_sites_in_eval %.2f", 100 * fractionDBSitesDiscoveredInEval()));

        s.add(String.format("n_snp_calls_at_indels    %d", nSNPsAtIndels()));
        s.add(String.format("percent_calls_at_indels  %.2f", nSNPsAtIndels() / (0.01 * nEvalSites())));

        return s;
    }
}
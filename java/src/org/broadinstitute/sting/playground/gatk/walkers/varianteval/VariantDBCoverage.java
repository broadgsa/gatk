package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.Utils;
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
    private int nEvalObs = 0;
    private int nSNPsAtdbSNPs = 0;
    private int nConcordant = 0;

    public VariantDBCoverage(final String name) {
        super("db_coverage");
        dbName = name;
    }

    public int nDBSNPs()        { return nDBSNPs; }
    public int nEvalSites()     { return nEvalObs; }
    public int nSNPsAtdbSNPs()  { return nSNPsAtdbSNPs; }
    public int nConcordant()    { return nConcordant; }
    public int nNovelSites()    { return Math.abs(nEvalSites() - nSNPsAtdbSNPs()); }

    public boolean discordantP(Variation dbSNP, Variation eval) {
        if (eval != null) {
            char alt = (eval.isSNP()) ? eval.getAlternativeBaseForSNP() : Utils.stringToChar(eval.getReference());
            if (dbSNP != null && dbSNP.isSNP())
                return !dbSNP.getAlleleList().contains(String.valueOf(alt));
        }
        return false;
    }


    /**
     * What fraction of the evaluated site variants were also found in the db?
     *
     * @return
     */
    public double dbSNPRate() {
        return nSNPsAtdbSNPs() / (1.0 * nEvalSites());
    }

    public double concordanceRate() {
        return nConcordant() / (1.0 * nSNPsAtdbSNPs());
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        rodDbSNP dbSNP = rodDbSNP.getFirstRealSNP(tracker.getTrackData( dbName, null ));
        String result = null;

        if (dbSNP != null) nDBSNPs++;               // count the number of real dbSNP events
        if ( eval != null && eval.isSNP() ) {       // ignore indels right now
            nEvalObs++;                             // count the number of eval snps we've seen

            if (dbSNP != null) {                    // both eval and dbSNP have real snps
                nSNPsAtdbSNPs++;

                if (!discordantP(dbSNP, eval))      // count whether we're concordant or not with the dbSNP value
                    nConcordant++;
                else
                    result = String.format("Discordant [DBSNP %s] [EVAL %s]", dbSNP, eval);
            }
        }

        if ( dbSNP != null && dbSNP.isSNP() ) {
            BrokenRODSimulator.attach("dbSNP");
            rodDbSNP dbsnp = (rodDbSNP) BrokenRODSimulator.simulate_lookup("dbSNP", context.getLocation(), tracker);
            if ( ! dbSNP.getRS_ID().equals(dbsnp.getRS_ID()) && dbsnp.isSNP() ) {
                System.out.printf("Discordant site! %n%s%n vs.%n%s%n", dbSNP, dbsnp);
            }
        }

        return result;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("name                     %s", dbName));

        s.add(String.format("n_db_snps                %d", nDBSNPs()));
        s.add(String.format("n_eval_sites             %d", nEvalSites()));
        s.add(String.format("n_overlapping_sites      %d", nSNPsAtdbSNPs()));
        s.add(String.format("n_concordant             %d", nConcordant()));
        s.add(String.format("n_novel_sites            %d", nNovelSites()));

        s.add(String.format("dbsnp_rate               %.2f       # percent eval snps at dbsnp snps", 100 * dbSNPRate()));
        s.add(String.format("concordance_rate         %.2f", 100 * concordanceRate()));

        return s;
    }
}
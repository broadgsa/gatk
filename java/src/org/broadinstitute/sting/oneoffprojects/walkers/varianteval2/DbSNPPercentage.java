package org.broadinstitute.sting.oneoffprojects.variantcontext.varianteval2;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.oneoffprojects.variantcontext.VariantContext;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Genotype;
import org.broadinstitute.sting.oneoffprojects.variantcontext.Allele;

import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class DbSNPPercentage extends VariantEvaluator {
    private long nDBSNPs       = 0;
    private long nEvalSNPs     = 0;
    private long nSNPsAtdbSNPs = 0;
    private long nConcordant   = 0;

    public DbSNPPercentage(VariantEval2Walker parent) {
        // don't do anything
    }

    public String getName() {
        return "dbsnp_percentage";
    }

    public int getComparisonOrder() {
        return 2;   // we need to see each eval track and each comp track
    }

    public long nDBSNPs()        { return nDBSNPs; }
    public long nEvalSNPs()      { return nEvalSNPs; }
    public long nSNPsAtdbSNPs()  { return nSNPsAtdbSNPs; }
    public long nConcordant()    { return nConcordant; }
    public long nNovelSites()    { return Math.abs(nEvalSNPs() - nSNPsAtdbSNPs()); }

    public String toString() {
        return getName() + ": " + summaryLine();
    }

    private String summaryLine() {
        return String.format("%d %d %d %d %d %.2f %.2f",
                nDBSNPs(), nEvalSNPs(), nSNPsAtdbSNPs(), nConcordant(), nNovelSites(), 100 * dbSNPRate(), 100 * concordanceRate());
    }

    private static List<String> HEADER =
            Arrays.asList("n_dbsnps", "n_eval_snps", "n_overlapping_snps", "n_concordant", "n_novel_snps", "dbsnp_rate", "concordance_rate");

    // making it a table
    public List<String> getTableHeader() {
        return HEADER;
    }

    public List<List<String>> getTableRows() {
        return Arrays.asList(Arrays.asList(summaryLine().split(" ")));
    }

    /**
     * Returns true if every allele in eval is also in dbsnp
     *
     * @param eval
     * @param dbsnp
     * @return
     */
    public boolean discordantP(VariantContext eval, VariantContext dbsnp ) {
        for ( Allele a : eval.getAlleles() ) {
            if ( ! dbsnp.hasAllele(a, true) )
                return true;
        }

        return false;
    }


    /**
     * What fraction of the evaluated site variants were also found in the db?
     *
     * @return
     */
    public double dbSNPRate()           { return rate(nSNPsAtdbSNPs(), nEvalSNPs()); }
    public double concordanceRate()     { return rate(nConcordant(), nSNPsAtdbSNPs()); }

    public void update2(VariantContext eval, VariantContext dbsnp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        boolean dbSNPIsGood = dbsnp != null && dbsnp.isSNP() && dbsnp.isNotFiltered();
        boolean evalIsGood = eval != null && eval.isSNP();

        if ( dbSNPIsGood ) nDBSNPs++;             // count the number of dbSNP events
        if ( evalIsGood )  nEvalSNPs++;           // count the number of dbSNP events

        if ( dbSNPIsGood && evalIsGood )  {
            nSNPsAtdbSNPs++;

            if ( ! discordantP(eval, dbsnp) )     // count whether we're concordant or not with the dbSNP value
                nConcordant++;
        }
    }
}

package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.genotype.Variation;

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
public class VariantMatcher extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    String dbName;

    public VariantMatcher(final String name) {
        super("variant_matches");
        dbName = name;
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        String r = null;
        Variation db = (Variation)tracker.lookup(dbName, null);
        
        if ( eval != null || db != null ) {
            String matchFlag = "    ";
            if ( eval != null && db != null ) matchFlag = "*** ";
            if ( eval == null && db != null ) matchFlag = ">>> ";
            if ( eval != null && db == null ) matchFlag = "<<< ";

            r = String.format("%s  %s: %s <=> EVAL: %s", dbName, matchFlag, db, eval);
        }
        return r;
    }
}
package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
import java.util.List;
import java.util.Arrays;
import java.util.ArrayList;

public class VariantMatcher extends BasicVariantAnalysis {
    String dbName;

    public VariantMatcher(final String name) {
        super("variant_matches");
        dbName = name;
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        String r = null;
        AllelicVariant db = (AllelicVariant)tracker.lookup(dbName, null);
        
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
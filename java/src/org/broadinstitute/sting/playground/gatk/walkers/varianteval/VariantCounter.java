package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;

import java.io.PrintStream;
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
public class VariantCounter extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    long nBasesCovered = 0;
    int nSNPs = 0;

    public VariantCounter() {
        super("variant_counts");
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        nSNPs += eval == null ? 0 : 1;
        return null;
    }

    /**
     * No need to finalize the data in general
     * @param nSites
     */
    public void finalize(long nSites) {
        nBasesCovered = nSites;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("n bases covered: %d", nBasesCovered));
        s.add(String.format("sites: %d", nSNPs));
        s.add(String.format("variant rate: %.5f confident variants per base", nSNPs / (1.0 * Math.max(nBasesCovered, 1))));
        s.add(String.format("variant rate: 1 / %d confident variants per base", nBasesCovered / Math.max(nSNPs, 1)));
        return s;
    }
}
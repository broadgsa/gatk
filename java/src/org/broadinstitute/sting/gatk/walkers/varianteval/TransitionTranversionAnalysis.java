package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.utils.BaseUtils;
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
public class TransitionTranversionAnalysis extends BasicVariantAnalysis implements GenotypeAnalysis, PopulationAnalysis {
    long nTransitions = 0, nTransversions = 0;
    long novelTransitions = 0, novelTransversions = 0;

    public TransitionTranversionAnalysis() {
        super("transitions_transversions");
    }

    public String update(Variation eval, RefMetaDataTracker tracker, char ref, AlignmentContext context) {
        if (eval != null && eval.isSNP()) {
            boolean isNovel = false;
            if( eval instanceof RodVCF && ((RodVCF)eval).isNovel() ) { isNovel = true; }
            final char altBase = eval.getAlternativeBaseForSNP();
            final BaseUtils.BaseSubstitutionType subType = BaseUtils.SNPSubstitutionType((byte)ref, (byte)altBase);
            if (subType == BaseUtils.BaseSubstitutionType.TRANSITION) {
                nTransitions++;
                if( isNovel ) { novelTransitions++; }
            } else {
                nTransversions++;
                if( isNovel ) { novelTransversions++; }
            }
        }
        return null;
    }

    public List<String> done() {
        List<String> s = new ArrayList<String>();
        s.add(String.format("transitions     %d", nTransitions));
        s.add(String.format("transversions   %d", nTransversions));
        if( novelTransitions != 0 || novelTransversions != 0 ) {
            s.add(String.format("overall ratio   %.2f", nTransitions / (1.0 * nTransversions)));
            s.add(String.format("known var ratio %.2f", (nTransitions-novelTransitions) / (1.0 * (nTransversions-novelTransversions))));
            s.add(String.format("novel var ratio %.2f", novelTransitions / (1.0 * novelTransversions)));
        } else {
            s.add(String.format("ratio           %.2f", nTransitions / (1.0 * nTransversions)));
        }
        return s;
    }
}

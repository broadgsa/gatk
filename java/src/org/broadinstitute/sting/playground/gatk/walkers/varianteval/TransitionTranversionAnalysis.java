package org.broadinstitute.sting.playground.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.BaseUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Jun 4, 2009
 * Time: 4:38:00 PM
 * To change this template use File | Settings | File Templates.
 */
public class TransitionTranversionAnalysis extends BasicVariantAnalysis {
    int N_TRANSITION_TRANVERSION_BINS = 100;
    Histogram<Integer> transitions;
    Histogram<Integer> transversions;

    public TransitionTranversionAnalysis() {
        super("transitions_transversions");
        transitions   = new Histogram<Integer>(N_TRANSITION_TRANVERSION_BINS, 0.0, 1.0, 0);
        transversions = new Histogram<Integer>(N_TRANSITION_TRANVERSION_BINS, 0.0, 1.0, 0);
    }

    public String update(AllelicVariant eval, RefMetaDataTracker tracker, char ref, LocusContext context) {
        if ( eval != null && eval.isSNP() ) {
            char refBase = eval.getRefSnpFWD();
            char altBase = eval.getAltSnpFWD();
            //System.out.printf("%c %c%n", refBase, altBase);
            //int i = transition_transversion_bin(dbsnp.getHeterozygosity());
            //System.out.printf("MAF = %f => %d%n", dbsnp.getMAF(), i);
            //EnumMap<BaseUtils.BaseSubstitutionType, Integer> bin = transition_transversion_counts[i];

            BaseUtils.BaseSubstitutionType subType = BaseUtils.SNPSubstitutionType(refBase, altBase);
            Histogram<Integer> h = subType == BaseUtils.BaseSubstitutionType.TRANSITION ? transitions : transversions;
            double het = eval.getHeterozygosity();
            h.setBin(het, h.getBin(het) + 1);
            //int sit = bin.get(BaseUtils.BaseSubstitutionType.TRANSITION);
            //int ver = bin.get(BaseUtils.BaseSubstitutionType.TRANSVERSION);
            //System.out.printf("%s %.2f %s%n", subType, dbsnp.getHeterozygosity(), h.x2bin(logHet), dbsnp.toString());
        }

        return null;
    }

    public List<String> done() {
        int nTransitions = 0;
        int nTransversions = 0;

        List<String> s = new ArrayList<String>();
        for ( int i = 0; i < N_TRANSITION_TRANVERSION_BINS; i++ ) {
            //double avHet = Math.pow(10, transitions.bin2x(i));
            double avHet = transitions.bin2x(i);
            if ( avHet > 0.5 ) break;

            int sit = transitions.getBin(i);
            int ver = transversions.getBin(i);
            nTransitions += sit;
            nTransversions += ver;
            double ratio = (float)sit/ver;
            //s.append(String.format("%s %s %.2f %d %d %f%n", commentLine, prefix, avHet, sit, ver, ratio));
        }

        s.add(String.format("transitions    %d", nTransitions));
        s.add(String.format("transversions  %d", nTransversions));
        s.add(String.format("ratio          %.2f", nTransitions / (1.0 * nTransversions)));
        return s;
    }
}

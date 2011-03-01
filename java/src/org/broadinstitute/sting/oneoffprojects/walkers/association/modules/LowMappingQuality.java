package org.broadinstitute.sting.oneoffprojects.walkers.association.modules;

import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.AssociationContext;
import org.broadinstitute.sting.oneoffprojects.walkers.association.AssociationContextAtom;
import org.broadinstitute.sting.oneoffprojects.walkers.association.interfaces.TStatistic;
import org.broadinstitute.sting.oneoffprojects.walkers.association.interfaces.ZStatistic;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 1:49:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class LowMappingQuality extends AssociationContext<LowMappingQualityAtom> implements ZStatistic, TStatistic {
    private double df = -1;

    public LowMappingQuality(Class<LowMappingQualityAtom> atomClass ) {
        super(atomClass);
    }

    public LowMappingQuality() {
        this(LowMappingQualityAtom.class);
    }

    public int getWindowSize() { return 20; }
    public int slideByValue() { return 5; }

    public double getZStatistic() {
        int case_low= 0;
        int case_tot = 0;
        int ctrl_low = 0;
        int ctrl_tot = 0;
        for ( LowMappingQualityAtom atom : window ) {
            Map<String, Pair<Integer,Integer>> ccc = atom.getCaseControlCounts();
            case_low += ccc.get("case").first;
            case_tot += ccc.get("case").second;
            ctrl_low += ccc.get("control").first;
            ctrl_tot += ccc.get("control").second;
        }
        double p_case = ( (double) case_low)/case_tot;
        double p_ctrl = ((double) ctrl_low)/ctrl_tot;
        double p_pool = ( (double) case_low + (double) ctrl_low )/(case_tot+ctrl_tot);
        double SE = Math.sqrt( p_pool*(1-p_pool)*( 1/((double)ctrl_tot) + 1/((double)case_tot) ));
        return (p_case - p_ctrl)/SE;
    }

    public double getTStatistic() {
        // todo -- running average is unnecessary and inefficient here, but was quick to write
        MathUtils.RunningAverage cases =  new MathUtils.RunningAverage();
        MathUtils.RunningAverage controls = new MathUtils.RunningAverage();
        Set<Sample> caseSamples = new HashSet<Sample>();
        Set<Sample> controlSamples = new HashSet<Sample>();
        for ( LowMappingQualityAtom atom : window ) {
            Map<String,List<Integer>> caseControlQuals = atom.getCaseControlMappingQuals();
            cases.addAll((Collection) caseControlQuals.get("case"));
            controls.addAll((Collection) caseControlQuals.get("control"));
            for ( Sample s : atom.getSamples() ) {
                if ( s.getProperty("cohort").equals("case") ) {
                    caseSamples.add(s);
                }
                if ( s.getProperty("cohort").equals("control")) {
                    controlSamples.add(s);
                }
            }
        }
        double a = cases.var()/caseSamples.size();
        double b = controls.var()/controlSamples.size();

        df = ( a + b )/( a*a/(caseSamples.size()-1) + b*b/(controlSamples.size()-1));

        return ( cases.mean() - controls.mean() )/Math.sqrt( (cases.var()/caseSamples.size()) + (controls.var()/controlSamples.size()));
    }

    // todo -- this is super hacky
    public double getDegreesOfFreedom() {
        if ( df < 0 ) {
            getTStatistic();
        }
        double holder = df;
        df = -1;
        return holder;
    }
}

package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import cern.jet.random.Normal;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.Dichotomizable;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.StatisticalTest;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 12:58:38 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ProportionTest extends CaseControl<Pair<Number,Number>> implements Dichotomizable, StatisticalTest {
    public static final Normal standardNormal = new Normal(0.0,1.0,null);

    public abstract Pair<Number,Number> map(ReadBackedPileup rbp );

    public Pair<Number,Number> add(Pair<Number,Number> left, Pair<Number,Number> right) {
        return new Pair<Number,Number>(left.first.doubleValue()+right.first.doubleValue(),
                                       left.second.doubleValue()+right.second.doubleValue());
    }

    public Pair<Collection<Number>,Collection<Number>> getDichotomizedData() {
        Collection<Number> caseProps = new ArrayList<Number>();
        Collection<Number> controlProps = new ArrayList<Number>();

        for (Map<Sample,Pair<Number,Number>> sampleNumberMap : window ) {
            for ( Map.Entry<Sample,Pair<Number,Number>> samplePairEntry : sampleNumberMap.entrySet() ) {
                Pair<Number,Number> val = samplePairEntry.getValue();
                if ( samplePairEntry.getKey().getProperty("cohort").equals("case")) {
                    caseProps.add(val.first.doubleValue()/val.second.doubleValue());
                } else if ( samplePairEntry.getKey().getProperty("cohort").equals("control") ) {
                    controlProps.add(val.first.doubleValue()/val.second.doubleValue());
                }
            }
        }

        return new Pair<Collection<Number>,Collection<Number>>(caseProps,controlProps);
    }

    public Pair<Double,Double> getStatisticAndPValue() {
        Map<CaseControl.Cohort,Pair<Number,Number>> caseControlCounts = getCaseControl();
        if ( caseControlCounts == null || caseControlCounts.get(CaseControl.Cohort.CASE) == null || caseControlCounts.get(CaseControl.Cohort.CONTROL) == null ) {
            return new Pair<Double,Double>(Double.NaN,Double.NaN);
        }
        double pCase = caseControlCounts.get(CaseControl.Cohort.CASE).first.doubleValue()/caseControlCounts.get(CaseControl.Cohort.CASE).second.doubleValue();
        double pControl = caseControlCounts.get(CaseControl.Cohort.CONTROL).first.doubleValue()/caseControlCounts.get(CaseControl.Cohort.CONTROL).second.doubleValue();
        double nCase = caseControlCounts.get(CaseControl.Cohort.CASE).second.doubleValue();
        double nControl = caseControlCounts.get(CaseControl.Cohort.CONTROL).second.doubleValue();

        double p2 = (caseControlCounts.get(CaseControl.Cohort.CASE).first.doubleValue()+caseControlCounts.get(CaseControl.Cohort.CONTROL).first.doubleValue())/
                (caseControlCounts.get(CaseControl.Cohort.CASE).second.doubleValue()+caseControlCounts.get(CaseControl.Cohort.CONTROL).second.doubleValue());
        double se = Math.sqrt(p2*(1-p2)*(1/nCase + 1/nControl));

        double z = (pCase-pControl)/se;
        double p = z < 0 ? 2*standardNormal.cdf(z) : 2*(1-standardNormal.cdf(z));

        return new Pair<Double,Double>(z,p);
    }

    public String getStatisticName() { return "Z"; }

}
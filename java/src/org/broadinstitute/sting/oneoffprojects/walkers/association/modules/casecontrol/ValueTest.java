package org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol;

import cern.jet.random.StudentT;
import org.broadinstitute.sting.gatk.datasources.sample.Sample;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.Dichotomizable;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.StatisticalTest;
import org.broadinstitute.sting.utils.MannWhitneyU;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 3/2/11
 * Time: 1:53 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class ValueTest extends CaseControl<Collection<Number>> implements StatisticalTest, Dichotomizable {

    public abstract Collection<Number> map(ReadBackedPileup rbp );

    public boolean useTStatistic() { return false; } // default to using the U statistic

    public Collection<Number> add(Collection<Number> left, Collection<Number> right) {
        if ( left instanceof ArrayList ) {
            ((ArrayList) left).addAll(right);
            return left;
        } else if ( left instanceof Set) {
            ((Set) left).addAll(right);
            return left;
        } else {
            List<Number> newList = new ArrayList<Number>(left.size()+right.size());
            newList.addAll(left);
            newList.addAll(right);
            return newList;
        }
    }

    public Pair<Collection<Number>,Collection<Number>> getDichotomizedData() {
        Collection<Number> caseMeans = new ArrayList<Number>();
        Collection<Number> controlMeans = new ArrayList<Number>();

        for ( Map<Sample,Collection<Number>> sampleMap : window ) {
            for ( Map.Entry<Sample,Collection<Number>> sampleEntry : sampleMap.entrySet() ) {
                if ( sampleEntry.getKey().getProperty("cohort").equals("case") ) {
                    caseMeans.add(MathUtils.average(sampleEntry.getValue()));
                } else if ( sampleEntry.getKey().getProperty("cohort").equals("control") ) {
                    controlMeans.add(MathUtils.average(sampleEntry.getValue()));
                }
            }
        }

        return new Pair<Collection<Number>,Collection<Number>>(caseMeans,controlMeans);
    }

    public Pair<Double,Double> getUStatisticAndPValue() {
        MannWhitneyU mwu = new MannWhitneyU();

        for ( Map.Entry<Cohort,Collection<Number>> cohortEntry : getCaseControl().entrySet() ) {
            if ( cohortEntry.getKey().equals(Cohort.CASE) ) {
                for ( Number n : cohortEntry.getValue() ) {
                    mwu.add(n,MannWhitneyU.USet.SET1);
                }
            } else if ( cohortEntry.getKey().equals(Cohort.CONTROL)) {
                for ( Number n : cohortEntry.getValue() ) {
                    mwu.add(n,MannWhitneyU.USet.SET2);
                }
            }
        }

        return mwu.runTwoSidedTest();
    }

    public Pair<Double,Double> getStatisticAndPValue() { return useTStatistic() ? getTStatisticAndPValue() : getUStatisticAndPValue(); }

    public String getStatisticName() { return useTStatistic() ? "T" : "V"; }

    public Pair<Double,Double> getTStatisticAndPValue() {
        Map<CaseControl.Cohort,Collection<Number>> caseControlVectors = getCaseControl();
        if ( caseControlVectors == null || caseControlVectors.get(CaseControl.Cohort.CASE) == null || caseControlVectors.get(CaseControl.Cohort.CONTROL) == null ) {
            return new Pair<Double,Double>(Double.NaN,Double.NaN);
        }
        double meanCase = MathUtils.average(caseControlVectors.get(CaseControl.Cohort.CASE));
        double varCase = MathUtils.variance(caseControlVectors.get(CaseControl.Cohort.CASE),meanCase);
        double nCase =  caseControlVectors.get(CaseControl.Cohort.CASE).size();
        double meanControl = MathUtils.average(caseControlVectors.get(CaseControl.Cohort.CONTROL));
        double varControl =  MathUtils.variance(caseControlVectors.get(CaseControl.Cohort.CONTROL),meanControl);
        double nControl = caseControlVectors.get(CaseControl.Cohort.CONTROL).size();

        double df_num = Math.pow(varCase/nCase + varControl/nControl,2);
        double df_denom = Math.pow(varCase/nCase,2)/(nCase-1) + Math.pow(varControl/nControl,2)/(nControl-1);
        double t = (meanCase-meanControl)/Math.sqrt(varCase/nCase+varControl/nControl);

        StudentT studentT = new StudentT(df_num/df_denom,null);
        double p = t < 0 ? 2*studentT.cdf(t) : 2*(1-studentT.cdf(t));

        return new Pair<Double,Double>(t,p);
    }

}

package org.broadinstitute.sting.oneoffprojects.walkers.association;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import cern.jet.random.Normal;
import cern.jet.random.StudentT;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.casecontrol.*;
import org.broadinstitute.sting.utils.MannWhitneyU;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.WilcoxonRankSum;
import org.broadinstitute.sting.utils.collections.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 11:48:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class AssociationTestRunner {
    static Normal standardNormal = new Normal(0.0,1.0,null);

    public static List<String> runTests(AssociationContext context) {
        List<String> results = new ArrayList<String>();
        if ( context instanceof TStatistic) {
            results.add(runStudentT((TStatistic) context));
        }

        if ( context instanceof ZStatistic) {
            results.add(runZ((ZStatistic) context));
        }

        if ( context instanceof UStatistic ) {
            results.add(runU((UStatistic) context));
        }

        return results;
    }

    public static String runStudentT(TStatistic context) {
        Map<CaseControl.Cohort,Collection<Number>> caseControlVectors = context.getCaseControl();
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
        return String.format("T: %.2f\tP: %.2e",t,p);
    }

    public static String runZ(ZStatistic context) {
        Map<CaseControl.Cohort,Pair<Number,Number>> caseControlCounts = context.getCaseControl();
        double pCase = caseControlCounts.get(CaseControl.Cohort.CASE).first.doubleValue()/caseControlCounts.get(CaseControl.Cohort.CASE).second.doubleValue();
        double pControl = caseControlCounts.get(CaseControl.Cohort.CONTROL).first.doubleValue()/caseControlCounts.get(CaseControl.Cohort.CONTROL).second.doubleValue();
        double nCase = caseControlCounts.get(CaseControl.Cohort.CASE).second.doubleValue();
        double nControl = caseControlCounts.get(CaseControl.Cohort.CONTROL).second.doubleValue();

        double p2 = (caseControlCounts.get(CaseControl.Cohort.CASE).first.doubleValue()+caseControlCounts.get(CaseControl.Cohort.CONTROL).first.doubleValue())/
                     (caseControlCounts.get(CaseControl.Cohort.CASE).second.doubleValue()+caseControlCounts.get(CaseControl.Cohort.CONTROL).first.doubleValue());
        double se = Math.sqrt(p2*(1-p2)*(1/nCase + 1/nControl));

        double z = (pCase-pControl)/se;
        double p = z < 0 ? 2*standardNormal.cdf(z) : 2*(1-standardNormal.cdf(z));
        return String.format("Z: %.2f\tP: %.2e",z,p);
    }

    public static String runU(UStatistic context) {
        Map<CaseControl.Cohort,Collection<Number>> caseControlVectors = context.getCaseControl();
        MannWhitneyU mwu = new MannWhitneyU();
        for ( Number n : caseControlVectors.get(CaseControl.Cohort.CASE) ) {
            mwu.add(n, MannWhitneyU.USet.SET1);
        }
        for ( Number n : caseControlVectors.get(CaseControl.Cohort.CONTROL) ) {
            mwu.add(n,MannWhitneyU.USet.SET2);
        }
        Pair<Integer,Double> results = mwu.runTwoSidedTest();
        return String.format("U: %d\tP: %.2e",results.first,results.second);
    }

    public static String runFisherExact(AssociationContext context) {
        return "Test not yet implemented";
    }
}

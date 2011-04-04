package org.broadinstitute.sting.oneoffprojects.walkers.association;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import cern.jet.random.Normal;
import cern.jet.random.StudentT;
import org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol.*;
import org.broadinstitute.sting.utils.MannWhitneyU;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.collections.Pair;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Feb 24, 2011
 * Time: 11:48:26 AM
 * To change this template use File | Settings | File Templates.
 */
public class AssociationTestRunner {
    final static int MAX_Q_VALUE = Integer.MAX_VALUE;
    // todo -- this was written when ACs could implement interfaces, now that they extend, there's no multiple inheritance
    static Normal standardNormal = new Normal(0.0,1.0,null);

    public static int pToQ(double p) {
        return Math.min((int) Math.floor(QualityUtils.phredScaleErrorRate(p)),MAX_Q_VALUE);
    }

    public static Pair<Double,Pair<Double,Integer>> getTestValues(AssociationContext context) {
        if ( context instanceof TStatistic) {
            Pair<Double,Double> t = testStudentT((TStatistic) context);
            return new Pair<Double,Pair<Double,Integer>> (t.first,new Pair<Double,Integer>(t.second,pToQ(t.second)));
        }

        if ( context instanceof ZStatistic) {
            Pair<Double,Double> z = testZ((ZStatistic) context);
            return new Pair<Double,Pair<Double,Integer>> (z.first,new Pair<Double,Integer>(z.second,pToQ(z.second)));
        }

        if ( context instanceof UStatistic ) {
            Pair<Double,Double> u = mannWhitneyUTest((UStatistic) context);
            return new Pair<Double,Pair<Double,Integer>> ( u.first, new Pair<Double,Integer>(u.second,pToQ(u.second)));
        }

        return null;
    }

    public static String runTests(AssociationContext context) {
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

        StringBuffer buf = new StringBuffer();
        if ( results.size() > 0 ) {
            buf.append(results.remove(0));
            for ( String s : results ) {
                buf.append('\t');
                buf.append(s);
            }
        }

        return buf.toString();
    }


    public static String runStudentT(TStatistic context) {
        Pair<Double,Double> stats = testStudentT(context);
        double t = stats.first;
        double p = stats.second;
        return String.format("T: %.2f\tP: %.2e\tQ: %d",t,p,pToQ(p));
    }

    public static Pair<Double,Double> testStudentT(TStatistic context) {
        Map<CaseControl.Cohort,Collection<Number>> caseControlVectors = context.getCaseControl();
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

    public static String runZ(ZStatistic context) {
        Pair<Double,Double> stats = testZ(context);
        double z = stats.first;
        double p = stats.second;
        return String.format("Z: %.2f\tP: %.2e\tQ: %d",z,p,pToQ(p));
    }

    public static Pair<Double,Double> testZ(ZStatistic context) {
        Map<CaseControl.Cohort,Pair<Number,Number>> caseControlCounts = context.getCaseControl();
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

    public static String runU(UStatistic context) {
        // note: u statistic (U) is relatively useless for recalibrating outside of the context of m and n
        // thus we report V = (U - (m*n+1)/2)/(n*m*(n+m+1)/12)
        Pair<Double,Double> results = mannWhitneyUTest(context);
        return String.format("V: %.2f\tP: %.2e\tQ: %d",results.first,results.second,pToQ(results.second));
    }

    public static Pair<Double,Double> mannWhitneyUTest(UStatistic context) {
        Map<CaseControl.Cohort,Collection<Number>> caseControlVectors = context.getCaseControl();
        if ( caseControlVectors == null || caseControlVectors.get(CaseControl.Cohort.CASE) == null || caseControlVectors.get(CaseControl.Cohort.CONTROL) == null ) {
            return new Pair<Double,Double>(Double.NaN,Double.NaN);
        }
        MannWhitneyU mwu = new MannWhitneyU();
        for ( Number n : caseControlVectors.get(CaseControl.Cohort.CASE) ) {
            mwu.add(n, MannWhitneyU.USet.SET1);
        }
        for ( Number n : caseControlVectors.get(CaseControl.Cohort.CONTROL) ) {
            mwu.add(n,MannWhitneyU.USet.SET2);
        }
        return mwu.runTwoSidedTest();
    }

    public static String runFisherExact(AssociationContext context) {
        return "Test not yet implemented";
    }
}

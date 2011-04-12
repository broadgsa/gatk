package org.broadinstitute.sting.oneoffprojects.walkers.association;

import java.util.Collection;
import java.util.Map;

import cern.jet.random.Normal;
import cern.jet.random.StudentT;
import org.broadinstitute.sting.oneoffprojects.walkers.association.modules.casecontrol.*;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.Dichotomizable;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.Dichotomizer1D;
import org.broadinstitute.sting.oneoffprojects.walkers.association.statistics.StatisticalTest;
import org.broadinstitute.sting.utils.MannWhitneyU;
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
    final static Dichotomizer1D.Transform LOG_TRANSFORM = (new Dichotomizer1D()).new Transform() {
        private double EPSILON = 1e-20;
        @Override
        public double transform(double d) {
            return Math.log(d+EPSILON);
        }
    };

    final static Dichotomizer1D.Transform ARCSINE_TRANSFORM = (new Dichotomizer1D()).new Transform() {
        @Override
        public double transform(double d) {
            return Math.asin(1.0-d);
        }
    };

    public static int pToQ(double p) {
        return Math.min((int) Math.floor(QualityUtils.phredScaleErrorRate(p)),MAX_Q_VALUE);
    }

    public static Pair<Double,Pair<Double,Integer>> getTestValues(AssociationContext context) {
        if ( context instanceof StatisticalTest ) {
            Pair<Double,Double> statAndP = ((StatisticalTest) context).getStatisticAndPValue();
            return new Pair<Double,Pair<Double,Integer>>(statAndP.first,
                    new Pair<Double,Integer>(statAndP.second,pToQ(statAndP.second)));
        }

        return null;
    }

    public static String runTests(AssociationContext context) {
        if ( context instanceof StatisticalTest ) {
            Pair<Double,Pair<Double,Integer>> results = getTestValues(context);
            return String.format("%s: %.2f\tP: %.2e\tQ: %d",
                    ((StatisticalTest) context).getStatisticName() ,
                    results.first,results.second.first,results.second.second);
        }

        return null;
    }

    @Deprecated // this is used for testing only at the moment
    public static Pair<Double,Double> mannWhitneyUTest(ValueTest context) {
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

    public static Pair<Double,Double> getDichotomizedValues(AssociationContext context) {
        if ( context instanceof Dichotomizable ) {
            double raw = Dichotomizer1D.simpleGaussianDichotomy(((Dichotomizable)context));
            double log = Dichotomizer1D.simpleGaussianDichotomy(((Dichotomizable)context),LOG_TRANSFORM);
            return new Pair<Double,Double>(raw,log);
        }

        return new Pair<Double,Double>(Double.NaN,Double.NaN);
    }
}

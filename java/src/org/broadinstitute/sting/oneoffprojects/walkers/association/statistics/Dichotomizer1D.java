package org.broadinstitute.sting.oneoffprojects.walkers.association.statistics;

import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.collections.Pair;

import java.util.ArrayList;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: 4/4/11
 * Time: 7:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class Dichotomizer1D {

    private Dichotomizer1D() { } // no instantiation

    public abstract class Transform {
        public abstract double transform(double n);

        public Collection<Number> apply(Collection<Number> collection) {
            ArrayList<Number> tData = new ArrayList<Number>(collection.size());
            for ( Number n : collection ) {
                tData.add(transform(n.doubleValue()));
            }

            return tData;
        }
    }

    public static double simpleGaussianDichotomy(Collection<Number> setOne, Collection<Number> setTwo) {
        double meanOne = MathUtils.sum(setOne)/setOne.size();
        double meanTwo = MathUtils.sum(setTwo)/setTwo.size();
        double stdOne = Math.sqrt(MathUtils.variance(setOne, meanOne));
        double stdTwo = Math.sqrt(MathUtils.variance(setTwo, meanTwo));

        return 1.0 - (3.0*(stdOne+stdTwo))/Math.abs(meanOne-meanTwo);
    }

    public static double simpleGaussianDichotomy(Collection<Number> setOne, Collection<Number> setTwo, Transform transform) {
        return simpleGaussianDichotomy(transform.apply(setOne),transform.apply(setTwo));
    }

    public static double simpleGaussianDichotomy(Dichotomizable dichotomizable) {
        Pair<Collection<Number>,Collection<Number>> dichotomizedData = dichotomizable.getDichotomizedData();
        return simpleGaussianDichotomy(dichotomizedData.first,dichotomizedData.second);
    }

    public static double simpleGaussianDichotomy(Dichotomizable dichotomizable, Transform transform) {
        Pair<Collection<Number>,Collection<Number>> dichotomizedData = dichotomizable.getDichotomizedData();
        return simpleGaussianDichotomy(dichotomizedData.first,dichotomizedData.second, transform);
    }

}

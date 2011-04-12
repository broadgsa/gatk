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

    public Dichotomizer1D() { } // instantiation required for access to transform

    public abstract class Transform {
        public abstract double transform(double n);

        public Collection<Number> apply(Collection<Number> collection) {
            ArrayList<Number> tData = new ArrayList<Number>(collection.size());
            for ( Number n : collection ) {
                tData.add(transform(n.doubleValue()));
            }

            return tData;
        }

        public Transform() {

        }
    }

    /**
     * Tries to measure effect size by the separation of setOne and setTwo clusters, with a window spread factor of 3
     * @param setOne - collection of data from set one
     * @param setTwo - colleciton of data from set two
     * @return - the so-called "Z"-factor (effect size/spread)
     */
    public static double simpleGaussianDichotomy(Collection<Number> setOne, Collection<Number> setTwo) {
        double meanOne = MathUtils.average(setOne,true);
        double meanTwo = MathUtils.average(setTwo,true);
        double stdOne = Math.sqrt(MathUtils.variance(setOne, meanOne,true));
        double stdTwo = Math.sqrt(MathUtils.variance(setTwo, meanTwo,true));
        /*
        System.out.print("setOne: ");
        for ( Number n : setOne ) {
            System.out.printf(",%.2f",n.doubleValue());
        }
        System.out.print("\tsetTwo: ");
        for ( Number n : setTwo ) {
            System.out.printf(",%.2f",n.doubleValue());
        }

        System.out.printf("\tmn1: %.2f mn2: %.2f var1: %.2f var2: %.2f%n",meanOne,meanTwo,stdOne,stdTwo);
        */
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

    public static Pair<Double,Double> twoMeansDichotomy() { return null; }

}

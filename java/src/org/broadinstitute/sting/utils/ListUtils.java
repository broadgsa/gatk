package org.broadinstitute.sting.utils;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: andrew
 * Date: Jul 2, 2009
 * Time: 3:12:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class ListUtils {

    static Random rand = new Random(12321); //System.currentTimeMillis());

    /**
     * Returns n random indices drawn with replacement from the range 0..(k-1)
     *
     * @param n  the number of random indices to draw (with replacement)
     * @param k  the total number of indices allowed
     * @return a list of random indices ranging from 0 to (k-1) with possible duplicates
     */
    static public ArrayList<Integer> sampleIndicesWithReplacement(int n, int k) {

        ArrayList<Integer> chosen_balls = new ArrayList <Integer>();
        for (int i=0; i<n; i++) {
            //Integer chosen_ball = balls[rand.nextInt(k)];
            chosen_balls.add(rand.nextInt(k));
            //balls.remove(chosen_ball);
        }

        return chosen_balls;
    }

    /**
     * Returns n random indices drawn without replacement from the range 0..(k-1)
     *
     * @param n  the number of random indices to draw (without replacement)
     * @param k  the total number of indices allowed
     * @return a list of random indices ranging from 0 to (k-1) without duplicates
     */
    static public ArrayList<Integer> sampleIndicesWithoutReplacement(int n, int k) {
        ArrayList<Integer> chosen_balls = new ArrayList<Integer>(k);

        for (int i = 0; i < n; i++) {
            chosen_balls.add(i);
        }

        Collections.shuffle(chosen_balls, rand);

        //return (ArrayList<Integer>) chosen_balls.subList(0, k);
        return new ArrayList<Integer>(chosen_balls.subList(0, k));
    }

    /**
     * Given a list of indices into a list, return those elements of the list with the possibility of drawing list elements multiple times

     * @param indices  the list of indices for elements to extract
     * @param list     the list from which the elements should be extracted
     * @param <T>      the template type of the ArrayList
     * @return         a new ArrayList consisting of the elements at the specified indices
     */
    static public <T> ArrayList<T> sliceListByIndices(List<Integer> indices, List<T> list) {
        ArrayList<T> subset = new ArrayList<T>();

        for (int i : indices) {
            subset.add(list.get(i));
        }
        
        return subset;
    }


}

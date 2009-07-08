package org.broadinstitute.sting.utils;

import java.util.List;
import java.util.HashSet;
import java.util.ArrayList;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: andrew
 * Date: Jul 2, 2009
 * Time: 3:12:57 PM
 * To change this template use File | Settings | File Templates.
 */
public class ListUtils {

    static Random rand = new Random(12321); //System.currentTimeMillis());

    static public ArrayList<Integer> sampleIndicesWithReplacement(int n, int k) {
        // Returns n random indices drawn with replacement from the range 1..k

        ArrayList<Integer> chosen_balls = new ArrayList <Integer>();
        for (int i=0; i<n; i++) {
            //Integer chosen_ball = balls[rand.nextInt(k)];
            chosen_balls.add(rand.nextInt(k));
            //balls.remove(chosen_ball);
        }

        return chosen_balls;
    }

    static public <T> ArrayList<T> sliceListByIndices(List<Integer> indices, List<T> list) {
        // Given a list of indices into a list, return those elements of the list with the possibility
        // of drawing list elements multiple times

        ArrayList<T> subset = new ArrayList<T>();

        for (int i : indices) {
            subset.add(list.get(i));
        }
        
        return subset;
    }


}

package org.broadinstitute.sting.utils;

import java.util.*;

import net.sf.samtools.SAMRecord;

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
     * @param n the total number of indices sampled from
     * @param k the number of random indices to draw (with replacement)
     * @return a list of k random indices ranging from 0 to (n-1) with possible duplicates
     */
    static public ArrayList<Integer> sampleIndicesWithReplacement(int n, int k) {

        ArrayList<Integer> chosen_balls = new ArrayList <Integer>(k);
        for (int i=0; i< k; i++) {
            //Integer chosen_ball = balls[rand.nextInt(k)];
            chosen_balls.add(rand.nextInt(n));
            //balls.remove(chosen_ball);
        }

        return chosen_balls;
    }

    /**
     * Returns n random indices drawn without replacement from the range 0..(k-1)
     *
     * @param n the total number of indices sampled from
     * @param k the number of random indices to draw (without replacement)
     * @return a list of k random indices ranging from 0 to (n-1) without duplicates
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

    public static Comparable orderStatisticSearch(int orderStat, List<Comparable> list) {
        // this finds the order statistic of the list (kth largest element)
        // the list is assumed *not* to be sorted

        final Comparable x = list.get(orderStat);
        ListIterator iterator = list.listIterator();
        ArrayList lessThanX = new ArrayList();
        ArrayList equalToX = new ArrayList();
        ArrayList greaterThanX = new ArrayList();

        for(Comparable y : list) {
            if(x.compareTo(y) > 0) {
                lessThanX.add(y);
            } else if(x.compareTo(y) < 0) {
                greaterThanX.add(y);
            } else
                equalToX.add(y);
        }

        if(lessThanX.size() > orderStat)
            return orderStatisticSearch(orderStat, lessThanX);
        else if(lessThanX.size() + equalToX.size() >= orderStat)
            return orderStat;
        else
            return orderStatisticSearch(orderStat - lessThanX.size() - equalToX.size(), greaterThanX);

    }


    public static Object getMedian(List<Comparable> list) {
        return orderStatisticSearch((int) Math.ceil(list.size()/2), list);
    }

    public static byte getQScoreOrderStatistic(List<SAMRecord> reads, List<Integer> offsets, int k) {
        // version of the order statistic calculator for SAMRecord/Integer lists, where the
        // list index maps to a q-score only through the offset index
        // returns the kth-largest q-score.

        if( reads.size() == 0) {
            return 0;
        }
        
        ArrayList lessThanQReads = new ArrayList();
        ArrayList equalToQReads = new ArrayList();
        ArrayList greaterThanQReads = new ArrayList();
        ArrayList lessThanQOffsets = new ArrayList();
        ArrayList greaterThanQOffsets = new ArrayList();

        final byte qk = reads.get(k).getBaseQualities()[offsets.get(k)];

        for(int iter = 0; iter < reads.size(); iter ++) {
            SAMRecord read = reads.get(iter);
            int offset = offsets.get(iter);
            byte quality = read.getBaseQualities()[offset];

            if(quality < qk) {
                lessThanQReads.add(read);
                lessThanQOffsets.add(offset);
            } else if(quality > qk) {
                greaterThanQReads.add(read);
                greaterThanQOffsets.add(offset);
            } else {
                equalToQReads.add(reads.get(iter));
            }
        }

        if(lessThanQReads.size() > k)
            return getQScoreOrderStatistic(lessThanQReads, lessThanQOffsets, k);
        else if(equalToQReads.size() + lessThanQReads.size() >= k)
            return qk;
        else
            return getQScoreOrderStatistic(greaterThanQReads, greaterThanQOffsets, k - lessThanQReads.size() - equalToQReads.size());

    }


    public static byte getQScoreMedian(List<SAMRecord> reads, List<Integer> offsets) {
        return getQScoreOrderStatistic(reads, offsets, (int)Math.floor(reads.size()/2.));
    }

}

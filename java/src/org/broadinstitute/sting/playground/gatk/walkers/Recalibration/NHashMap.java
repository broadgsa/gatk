package org.broadinstitute.sting.playground.gatk.walkers.Recalibration;

import java.util.*;

/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Oct 30, 2009
 */

public class NHashMap<T> extends HashMap<List<? extends Comparable>, T> {

	private static final long serialVersionUID = 1L; //BUGBUG: what should I do here?
    private ArrayList<ArrayList<Comparable>> keyLists;

    public NHashMap() {
        super();
        keyLists = new ArrayList<ArrayList<Comparable>>();
    }

    public NHashMap( int initialCapacity, float loadingFactor ) {
        super( initialCapacity, loadingFactor );
        keyLists = new ArrayList<ArrayList<Comparable>>();
    }

    /*

    // This method is here only to help facilitate direct comparison of old and refactored recalibrator.
    // The old recalibrator prints out it's mappings in a sorted order but the refactored recalibrator doesn't need to.
    // Comment out this overriding method to improve performance
    public T put(List<? extends Comparable> key, T value) {

        ArrayList<Comparable> thisList;
        for( int iii = 0; iii < key.size(); iii++ ) {
            thisList = keyLists.get( iii );
            if( thisList == null ) {
                thisList = new ArrayList<Comparable>();
            }
            if( !thisList.contains( key.get( iii ) ) ) {
                thisList.add( key.get(iii ) );
            }
        }
        return super.put( key, value );
    }

    // This method is here only to help facilitate direct comparison of old and refactored recalibrator.
    // The old recalibrator prints out it's mappings in a sorted order but the refactored recalibrator doesn't need to.
    // Comment out this overriding method to improve performance
    public Set<Map.Entry<List<? extends Comparable>, T>> entrySet() {
        
        for( ArrayList<Comparable> list : keyLists ) {
            Collections.sort(list);
        }

        for( ArrayList<Comparable> list : keyLists ) {
            for( Comparable comp : list ) {
                // BUGBUG: unfinished
            }
        }
        
        return null;
    }

    */

	public static <T extends Comparable> List<T> makeList(T... args) {
        List<T> list = new ArrayList<T>();
        for (T arg : args)
        {
            list.add(arg);
        }
        return list;
    }
}






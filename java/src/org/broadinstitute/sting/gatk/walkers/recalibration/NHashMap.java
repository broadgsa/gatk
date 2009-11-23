package org.broadinstitute.sting.gatk.walkers.recalibration;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;

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
 *
 * A HashMap that maps a list of comparables to any object <T>.
 * There is functionality for the mappings to be given back to you in sorted order.
 */

public class NHashMap<T> extends HashMap<List<? extends Comparable>, T> {

	private static final long serialVersionUID = 1L; //BUGBUG: what should I do here?, added by Eclipse
    private ArrayList<ArrayList<Comparable>> keyLists;

    public NHashMap() {
        super();
        keyLists = null;
    }

    public NHashMap( int initialCapacity, float loadingFactor ) {
        super( initialCapacity, loadingFactor );
        keyLists = null;
    }


    // This method is here only to help facilitate direct comparison of old and refactored recalibrator.
    // The old recalibrator prints out its mappings in a sorted order but the refactored recalibrator doesn't need to.
    public T myPut(List<? extends Comparable> key, T value) {

        if( keyLists == null ) {
            keyLists = new ArrayList<ArrayList<Comparable>>();
            for( Comparable comp : key ) {
                keyLists.add( new ArrayList<Comparable>() );
            }
        }

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

    // This method is very ugly but is here only to help facilitate direct comparison of old and refactored recalibrator.
    // The old recalibrator prints out its mappings in a sorted order but the refactored recalibrator doesn't need to.
    @SuppressWarnings(value = "unchecked")
    public ArrayList<Pair<List<? extends Comparable>, T>> entrySetSorted4() {

        ArrayList<Pair<List<? extends Comparable>, T>> theSet = new ArrayList<Pair<List<? extends Comparable>, T>>();

        for( ArrayList<Comparable> list : keyLists ) {
            Collections.sort(list);
        }

        if( keyLists.size() != 4 ) {
            throw new StingException("Are you sure you want to be calling this ugly method? NHashMap.entrySetSorted4()");
        }

        ArrayList<Comparable> newKey = null;
        for( Comparable c0 : keyLists.get(0) ) {
            for( Comparable c1 : keyLists.get(1) ) {
                for( Comparable c2 : keyLists.get(2) ) {
                    for( Comparable c3 : keyLists.get(3) ) {
                        newKey = new ArrayList<Comparable>();
                        newKey.add(c0);
                        newKey.add(c1);
                        newKey.add(c2);
                        newKey.add(c3);
                        T value = this.get( newKey );
                        if( value!= null ) {
                            theSet.add(new Pair<List<? extends Comparable>,T>( newKey, value ) );
                        }
                    }
                }
            }
        }

        return theSet;
    }

	public static <T> List<T> makeList(T... args) {
        List<T> list = new ArrayList<T>();
        for (T arg : args)
        {
            list.add(arg);
        }
        return list;
    }
}






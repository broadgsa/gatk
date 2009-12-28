package org.broadinstitute.sting.utils;

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
 * A HashMap that maps a list of comparables to any Object <T>.
 * There is functionality for the mappings to be given back to you in sorted order.
 */

public class NHashMap<T> extends HashMap<List<? extends Comparable>, T> {

	private static final long serialVersionUID = 1L; // Added by Eclipse
    private ArrayList<ArrayList<Comparable>> keyLists;

    public NHashMap() {
        super();
        keyLists = null;
    }

    public NHashMap( int initialCapacity, float loadingFactor ) {
        super( initialCapacity, loadingFactor );
        keyLists = null;
    }


    // This method is here only to help facilitate outputting the mappings in sorted order
    public T sortedPut(List<? extends Comparable> key, T value) {

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

    public ArrayList<Pair<List<? extends Comparable>, T>> entrySetSorted() {

        ArrayList<Pair<List<? extends Comparable>, T>> theSet = new ArrayList<Pair<List<? extends Comparable>, T>>();

        for( ArrayList<Comparable> list : keyLists ) {
            Collections.sort(list);
        }

        int[] keyIndex = new int[ keyLists.size() ];
        int[] maxIndex = new int[ keyLists.size() ];
        for( int iii = 0; iii < keyLists.size(); iii++ ) {
            keyIndex[iii] = 0;
            maxIndex[iii] = keyLists.get(iii).size();
        }

        // Try all the possible keys in sorted order, add them to the output set if they are in the hashMap
        boolean triedAllKeys = false;
        ArrayList<Comparable> newKey = null;
        while( !triedAllKeys ) {
            newKey = new ArrayList<Comparable>();
            for( int iii = 0; iii < keyLists.size(); iii++ ) {
                newKey.add(keyLists.get(iii).get(keyIndex[iii]));
            }
            T value = this.get( newKey );
            if( value!= null ) {
                theSet.add(new Pair<List<? extends Comparable>,T>( newKey, value ) );
            }

            // Increment the keyIndex
            keyIndex[keyLists.size() - 1]++;
            for( int iii = keyLists.size() - 1; iii >= 0; iii-- ) {
                if( keyIndex[iii] >= maxIndex[iii] ) { // Carry it forward
                    keyIndex[iii] = 0;
                    if( iii > 0 ) {
                        keyIndex[iii-1]++;
                    } else {
                        triedAllKeys = true;
                        break;
                    }
                } else {
                    break;
                }
            }
        }
        return theSet;
    }

    // Used to make the key from a list of <T> objects
	public static <T> List<T> makeList(T... args) {
        List<T> list = new ArrayList<T>();
        for( T arg : args )
        {
            list.add(arg);
        }
        return list;
    }
}






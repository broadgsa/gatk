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
 * Date: Dec 29, 2009
 */

public class NestedHashMap{

    private final Map data = new HashMap<Object, Object>();
    private ArrayList<ArrayList<Comparable>> keyLists; // Used to output the mappings in sorted order

    public Object get( final Object... keys ) {
        Map map = this.data;
        for( int iii = 0; iii < keys.length; iii++ ) {
            if( iii == keys.length - 1 ) {
                return map.get(keys[iii]);
            }
            else {
                map = (Map) map.get(keys[iii]);
                if( map == null ) { return null; }
            }
        }

        return null;
    }

    public void put( final Object value, final Object... keys ) {

        if( keyLists == null ) {
            keyLists = new ArrayList<ArrayList<Comparable>>();
            for( Object obj : keys ) {
                keyLists.add( new ArrayList<Comparable>() );
            }
        }

        ArrayList<Comparable> thisList;
        for( int iii = 0; iii < keys.length; iii++ ) {
            thisList = keyLists.get( iii );
            if( thisList == null ) {
                thisList = new ArrayList<Comparable>();
            }
            if( !thisList.contains( (Comparable)keys[iii] ) ) {
                thisList.add( (Comparable)keys[iii] );
            }
        }


        Map map = this.data;
        for( int iii = 0; iii < keys.length; iii++ ) {
            if( iii == keys.length - 1 ) {
                map.put(keys[iii], value);
            }
            else {
                Map tmp = (Map) map.get(keys[iii]);
                if( tmp == null ) {
                    tmp = new HashMap();
                    map.put(keys[iii], tmp);
                }
                map = tmp;
            }
        }
    }

    public ArrayList<Pair<Object[], Object>> entrySetSorted() {

        ArrayList<Pair<Object[], Object>> theSet = new ArrayList<Pair<Object[], Object>>();

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
        // BUGBUG: make this more efficient
        boolean triedAllKeys = false;
        ArrayList<Object> newKey = null;
        while( !triedAllKeys ) {
            newKey = new ArrayList<Object>();
            for( int iii = 0; iii < keyLists.size(); iii++ ) {
                newKey.add(keyLists.get(iii).get(keyIndex[iii]));
            }
            Object value = this.get( newKey.toArray() );
            if( value!= null ) {
                theSet.add(new Pair<Object[], Object>( newKey.toArray(), value ) );
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

}

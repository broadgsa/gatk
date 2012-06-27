/*
 * Copyright (c) 2010 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.collections;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Dec 29, 2009
 */

public class NestedHashMap {

    public final Map data = new HashMap<Object, Object>();

    public Object get( final Object... keys ) {
        Map map = this.data;
        final int keysLength = keys.length;
        for( int iii = 0; iii < keysLength; iii++ ) {
            if( iii == keysLength - 1 ) {
                return map.get(keys[iii]);
            } else {
                map = (Map) map.get(keys[iii]);
                if( map == null ) { return null; }
            }
        }

        return null;
    }

    public synchronized void put( final Object value, final Object... keys ) { // WARNING! value comes before the keys!
        this.put(value, false, keys );
    }

    public synchronized Object put( final Object value, boolean keepOldBindingIfPresent, final Object... keys ) {
        Map map = this.data;
        final int keysLength = keys.length;
        for( int iii = 0; iii < keysLength; iii++ ) {
            if( iii == keysLength - 1 ) {
                if ( keepOldBindingIfPresent && map.containsKey(keys[iii]) ) {
                    // this code test is for parallel protection when you call put() multiple times in different threads
                    // to initialize the map.  It returns the already bound key[iii] -> value
                    return map.get(keys[iii]);
                } else {
                    // we are a new binding, put it in the map
                    map.put(keys[iii], value);
                    return value;
                }
            } else {
                Map tmp = (Map) map.get(keys[iii]);
                if( tmp == null ) {
                    tmp = new HashMap();
                    map.put(keys[iii], tmp);
                }
                map = tmp;
            }
        }

        return value; // todo -- should never reach this point
    }

    public List<Object> getAllValues() {
        List<Object> result = new ArrayList<Object>();
        fillAllValues(data, result);
        return result;
    }

    private void fillAllValues(final Map map, final List<Object> result) {
        for ( Object value : map.values() ) {
            if ( value == null )
                continue;
            if ( value instanceof Map )
                fillAllValues((Map)value, result);
            else
                result.add(value);
        }
    }

    public static class Leaf {
        public final List<Object> keys;
        public final Object value;

        public Leaf(final List<Object> keys, final Object value) {
            this.keys = keys;
            this.value = value;
        }
    }

    public List<Leaf> getAllLeaves() {
        List<Leaf> result = new ArrayList<Leaf>();
        List<Object> path = new ArrayList<Object>();
        fillAllLeaves(data, path, result);
        return result;
    }

    private void fillAllLeaves(final Map map, final List<Object> path, final List<Leaf> result) {
        for ( final Object key : map.keySet() ) {
            final Object value = map.get(key);
            if ( value == null )
                continue;
            final List<Object> newPath = new ArrayList<Object>(path);
            newPath.add(key);
            if ( value instanceof Map ) {
                fillAllLeaves((Map) value, newPath, result);
            } else {
                result.add(new Leaf(newPath, value));
            }
        }
    }
}

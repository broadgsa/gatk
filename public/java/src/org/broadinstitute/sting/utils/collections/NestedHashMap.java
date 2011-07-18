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

import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Dec 29, 2009
 */

public class NestedHashMap{

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
}

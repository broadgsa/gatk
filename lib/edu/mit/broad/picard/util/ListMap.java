package edu.mit.broad.picard.util;

import java.util.List;
import java.util.HashMap;
import java.util.ArrayList;

/**
 * A Map class that holds a list of entries under each key instead of a single entry, and
 * provides utility methods for adding an entry under a key.
 *
 * @author Tim Fennell
 */
public class ListMap<K,V> extends HashMap<K, List<V>> {
    /** Adds a single value to the list stored under a key. */
    public void add(K key, V value) {
        List<V> values = get(key);
        if (values == null) {
            values = new ArrayList<V>();
            put(key, values);
        }

        values.add(value);
    }
}

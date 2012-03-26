package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import java.util.Map;
import java.util.TreeMap;

/**
 * A final constant class representing the specific state configuration
 * for a VariantEvaluator instance.
 *
 * TODO optimizations to entirely remove the TreeMap and just store the HashMap for performance and use the tree for the sorted tostring function.
 */
public final class StateKey {
    /** High-performance cache of the toString operation for a constant class */
    private final String string;
    private final TreeMap<String, String> states;

    public StateKey(final Map<String, String> states) {
        this.states = new TreeMap<String, String>(states);
        this.string = formatString();
    }

    public StateKey(final StateKey toOverride, final String keyOverride, final String valueOverride) {
        if ( toOverride == null ) {
            this.states = new TreeMap<String, String>();
        } else {
            this.states = new TreeMap<String, String>(toOverride.states);
        }

        this.states.put(keyOverride, valueOverride);
        this.string = formatString();
    }

    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final StateKey stateKey = (StateKey) o;

        if (states != null ? !states.equals(stateKey.states) : stateKey.states != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        return states.hashCode();
    }

    @Override
    public String toString() {
        return string;
    }

    private final String formatString() {
        StringBuilder b = new StringBuilder();
        
        for ( Map.Entry<String, String> entry : states.entrySet() ) {
            b.append(String.format("%s:%s;", entry.getKey(), entry.getValue()));
        }

        return b.toString();
    }

    // TODO -- might be slow because of tree map
    public String get(final String key) {
        return states.get(key);
    }
}

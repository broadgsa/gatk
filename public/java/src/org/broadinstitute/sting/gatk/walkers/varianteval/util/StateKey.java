package org.broadinstitute.sting.gatk.walkers.varianteval.util;

import java.util.Map;
import java.util.TreeMap;

/**
 * A final constant class representing the specific state configuration
 * for a VariantEvaluator instance.
 *
 * The way this is currently implemented is by a map from the name of a VariantStratification to a
 * specific state string.  For example, the stratification Novelty has states all, known, novel. A
 * specific variant and comp would be tagged as "known" by the stratification, and this could be
 * represented here by the map (Novelty -> known).
 *
 * TODO -- PERFORMANCE PROBLEM -- MAD 03/27/12
 * TODO -- PERFORMANCE PROBLEM -- MAD 03/27/12
 * TODO -- PERFORMANCE PROBLEM -- MAD 03/27/12
 * TODO -- PERFORMANCE PROBLEM -- MAD 03/27/12
 * TODO -- PERFORMANCE PROBLEM -- MAD 03/27/12
 *
 * I've been staring at this state key code for a while.  It's just not right, and expensive to boot.
 * Here are my thoughts for future work.  The state key is both a key with specific state values for
 * every stratification.  For example, (known, sample1, ac=1).  This capability is used in some places,
 * such as below, to return a set of all states that should be updated given the eval and comp
 * VCs.  In principle there are a finite set of such combinations (the product of all states for all active
 * stratifications at initialization).  We could represent such keys as integers into the set of all combinations.
 *
 * Note that all of the code that manipulates these things is just terrible. It's all string manipulation and
 * HashMaps.  Since we are effectively always squaring off our VE analyses (i.e., we have a table with
 * all variable values for all stratification combinations) it doesn't make sense to allow so much dynamicism.  Instead
 * we should just upfront create a giant table indexed by integer keys, and manage data via a simple map from
 * specific strat state to this key.
 *
 * The reason this is so important is that >80% of the runtime of VE with VCFs with >1000 samples is spent in
 * the initializeStateKey function.  Instead, we should have code that looks like:
 *
 * init:
 * allStates <- initializeCombinationalStateSpace
 *
 * map:
 * for each eval / comp pair:
 *   for each relevantState based on eval / comp:
 *     allStates[relevantState].update(eval, comp)
 *
 *
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

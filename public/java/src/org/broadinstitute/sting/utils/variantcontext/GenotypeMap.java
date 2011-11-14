/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext;

import java.util.*;

/**
 *
 */
public class GenotypeMap implements Map<String, Genotype> {
    final TreeMap<String, Genotype> genotypes;
    boolean immutable = false;
    public final static GenotypeMap NO_GENOTYPES = new GenotypeMap();

    // ---------------------------------------------------------------------------
    //
    // private constructors -- you have to use static create methods to make these classes
    //
    // ---------------------------------------------------------------------------

    private GenotypeMap() {
        this(false);
    }

    private GenotypeMap(boolean immutable) {
        this(new TreeMap<String, Genotype>(), immutable);
    }

    private GenotypeMap(final TreeMap<String, Genotype> genotypes, final boolean immutable) {
        this.genotypes = genotypes;
        this.immutable = immutable;
    }

    // ---------------------------------------------------------------------------
    //
    // public static factory methods
    //
    // ---------------------------------------------------------------------------

    public static final GenotypeMap create() {
        return new GenotypeMap();
    }

    public static final GenotypeMap create(final int nGenotypes) {
        return new GenotypeMap();
    }

    public static final GenotypeMap create(final GenotypeMap genotypes) {
        return create(genotypes.values());
    }

    // todo -- differentiate between empty constructor and copy constructor
    // todo -- create constructor (Genotype ... genotypes)

    public static final GenotypeMap create(final Map<String, Genotype> genotypes) {
        return create(genotypes.values());
    }

    public static final GenotypeMap create(final Collection<Genotype> genotypes) {
        if ( genotypes == null )
            return null; // todo -- really should return an empty map
        else {
            GenotypeMap genotypeMap = new GenotypeMap().mutable();
            for ( final Genotype g : genotypes ) {
                if ( genotypeMap.containsKey(g.getSampleName() ) )
                    throw new IllegalArgumentException("Duplicate genotype added to VariantContext: " + g);
                genotypeMap.put(g.getSampleName(), g);
            }

            //return genotypeMap.immutable();  // todo enable when we have time to dive into mutability issue
            return genotypeMap;
        }
    }

    // ---------------------------------------------------------------------------
    //
    // Mutability methods
    //
    // ---------------------------------------------------------------------------

    public final GenotypeMap mutable() {
        immutable = false;
        return this;
    }

    public final GenotypeMap immutable() {
        immutable = true;
        return this;
    }

    public boolean isMutable() {
        return ! immutable;
    }

    public final void checkImmutability() {
        if ( immutable )
            throw new IllegalAccessError("GenotypeMap is currently immutable, but a mutator method was invoked on it");
    }

    // ---------------------------------------------------------------------------
    //
    // Map methods
    //
    // ---------------------------------------------------------------------------

    @Override
    public void clear() {
        checkImmutability();
        genotypes.clear();
    }

    @Override
    public int size() {
        return genotypes.size();
    }

    @Override
    public boolean isEmpty() {
        return genotypes.isEmpty();
    }

    @Override
    public boolean containsKey(final Object o) {
        return genotypes.containsKey(o);
    }

    @Override
    public boolean containsValue(final Object o) {
        return genotypes.containsValue(o);
    }

    @Override
    public Genotype get(final Object o) {
        return genotypes.get(o);
    }

    @Override
    public Genotype put(final String s, final Genotype genotype) {
        checkImmutability();
        return genotypes.put(s, genotype);
    }

    @Override
    public Genotype remove(final Object o) {
        checkImmutability();
        return genotypes.remove(o);
    }

    @Override
    public void putAll(final Map<? extends String, ? extends Genotype> map) {
        checkImmutability();
        genotypes.putAll(map);
    }

    @Override
    public Set<String> keySet() {
        return Collections.unmodifiableSet(genotypes.keySet());
    }

    @Override
    public Collection<Genotype> values() {
        return genotypes.values();
    }

    @Override
    public Set<Entry<String, Genotype>> entrySet() {
        return genotypes.entrySet();
    }
}

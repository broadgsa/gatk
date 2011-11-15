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
public class GenotypeCollection implements List<Genotype> {
    public final static GenotypeCollection NO_GENOTYPES = new GenotypeCollection();

    Map<String, Integer> sampleNameToOffset = null;
    boolean cacheIsInvalid = true;
    final ArrayList<Genotype> genotypes;
    boolean immutable = false;

    // ---------------------------------------------------------------------------
    //
    // private constructors -- you have to use static create methods to make these classes
    //
    // ---------------------------------------------------------------------------

    private GenotypeCollection() {
        this(10, false);
    }

    private GenotypeCollection(final int n, final boolean immutable) {
        this(new ArrayList<Genotype>(n), immutable);
    }

    private GenotypeCollection(final ArrayList<Genotype> genotypes, final boolean immutable) {
        this.genotypes = genotypes;
        this.immutable = immutable;
    }

    // ---------------------------------------------------------------------------
    //
    // public static factory methods
    //
    // ---------------------------------------------------------------------------

    public static final GenotypeCollection create() {
        return new GenotypeCollection();
    }

    public static final GenotypeCollection create(final int nGenotypes) {
        return new GenotypeCollection(nGenotypes, false);
    }

    public static final GenotypeCollection create(final ArrayList<Genotype> genotypes) {
        return genotypes == null ? NO_GENOTYPES : new GenotypeCollection(genotypes, false);
    }

    public static final GenotypeCollection create(final Genotype... genotypes) {
        return new GenotypeCollection(new ArrayList<Genotype>(Arrays.asList(genotypes)), false);
    }

    public static final GenotypeCollection copy(final GenotypeCollection toCopy) {
        return create(new ArrayList<Genotype>(toCopy.genotypes));
    }

    public static final GenotypeCollection copy(final Collection<Genotype> toCopy) {
        return toCopy == null ? NO_GENOTYPES : create(new ArrayList<Genotype>(toCopy));
    }

//    public static final GenotypeMap create(final Collection<Genotype> genotypes) {
//        if ( genotypes == null )
//            return null; // todo -- really should return an empty map
//        else {
//            GenotypeMap genotypeMap = new GenotypeMap(genotypes.size(), false);
//            for ( final Genotype g : genotypes ) {
//                if ( genotypeMap.containsKey(g.getSampleName() ) )
//                    throw new IllegalArgumentException("Duplicate genotype added to VariantContext: " + g);
//                genotypeMap.put(g.getSampleName(), g);
//            }
//
//            //return genotypeMap.immutable();  // todo enable when we have time to dive into mutability issue
//            return genotypeMap;
//        }
//    }

    // ---------------------------------------------------------------------------
    //
    // Mutability methods
    //
    // ---------------------------------------------------------------------------

    public final GenotypeCollection mutable() {
        immutable = false;
        return this;
    }

    public final GenotypeCollection immutable() {
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
    // caches
    //
    // ---------------------------------------------------------------------------

    private void invalidateCaches() {
        cacheIsInvalid = true;
        if ( sampleNameToOffset != null ) sampleNameToOffset.clear();
    }

    private void buildCache() {
        cacheIsInvalid = false;

        if ( sampleNameToOffset == null )
            sampleNameToOffset = new HashMap<String, Integer>(genotypes.size());

        for ( int i = 0; i < genotypes.size(); i++ )
            sampleNameToOffset.put(genotypes.get(i).getSampleName(), i);
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
    public boolean add(final Genotype genotype) {
        checkImmutability();
        invalidateCaches();
        return genotypes.add(genotype);
    }

    public boolean add(final Genotype ... genotype) {
        checkImmutability();
        invalidateCaches();
        return genotypes.addAll(Arrays.asList(genotype));
    }

    @Override
    public void add(final int i, final Genotype genotype) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean addAll(final Collection<? extends Genotype> genotypes) {
        checkImmutability();
        invalidateCaches();
        return this.genotypes.addAll(genotypes);
    }

    @Override
    public boolean addAll(final int i, final Collection<? extends Genotype> genotypes) {
        throw new UnsupportedOperationException();
    }

    @Override
    public boolean contains(final Object o) {
        return this.genotypes.contains(o);
    }

    @Override
    public boolean containsAll(final Collection<?> objects) {
        return this.genotypes.containsAll(objects);
    }

    @Override
    public Genotype get(final int i) {
        return genotypes.get(i);
    }

    public Genotype get(final String sampleName) {
        buildCache();
        Integer offset = sampleNameToOffset.get(sampleName);
        return offset == null ? null : genotypes.get(offset);
    }

    @Override
    public int indexOf(final Object o) {
        return genotypes.indexOf(o);
    }

    @Override
    public Iterator<Genotype> iterator() {
        return genotypes.iterator();
    }

    @Override
    public int lastIndexOf(final Object o) {
        return genotypes.lastIndexOf(o);
    }

    @Override
    public ListIterator<Genotype> listIterator() {
        // todo -- must be immutable
        return genotypes.listIterator();
    }

    @Override
    public ListIterator<Genotype> listIterator(final int i) {
        // todo -- must be immutable
        return genotypes.listIterator(i);
    }

    @Override
    public Genotype remove(final int i) {
        checkImmutability();
        invalidateCaches();
        return genotypes.remove(i);
    }

    @Override
    public boolean remove(final Object o) {
        checkImmutability();
        invalidateCaches();
        return genotypes.remove(o);
    }

    @Override
    public boolean removeAll(final Collection<?> objects) {
        checkImmutability();
        invalidateCaches();
        return genotypes.removeAll(objects);
    }

    @Override
    public boolean retainAll(final Collection<?> objects) {
        checkImmutability();
        invalidateCaches();
        return genotypes.retainAll(objects);
    }

    @Override
    public Genotype set(final int i, final Genotype genotype) {
        checkImmutability();
        invalidateCaches();
        return genotypes.set(i, genotype);
    }

    @Override
    public List<Genotype> subList(final int i, final int i1) {
        return genotypes.subList(i, i1);
    }

    @Override
    public Object[] toArray() {
        return genotypes.toArray();
    }

    @Override
    public <T> T[] toArray(final T[] ts) {
        return genotypes.toArray(ts);
    }

    public Iterable<Genotype> iterateInSampleNameOrder(final Iterable<String> sampleNamesInOrder) {
        return new Iterable<Genotype>() {
            @Override
            public Iterator<Genotype> iterator() {
                return new InOrderIterator(sampleNamesInOrder.iterator());
            }
        };
    }

    public Iterable<Genotype> iterateInSampleNameOrder() {
        return iterateInSampleNameOrder(getSampleNamesOrderedByName());
    }

    private final class InOrderIterator implements Iterator<Genotype> {
        final Iterator<String> sampleNamesInOrder;

        private InOrderIterator(final Iterator<String> sampleNamesInOrder) {
            this.sampleNamesInOrder = sampleNamesInOrder;
        }

        @Override
        public boolean hasNext() {
            return sampleNamesInOrder.hasNext();
        }

        @Override
        public Genotype next() {
            return get(sampleNamesInOrder.next());
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    public Set<String> getSampleNames() {
        buildCache();
        return sampleNameToOffset.keySet();
    }

    public Set<String> getSampleNamesOrderedByName() {
        return new TreeSet<String>(getSampleNames());
    }

    public boolean containsSample(final String sample) {
        buildCache();
        return sampleNameToOffset.containsKey(sample);
    }

    public boolean containsSamples(final Collection<String> samples) {
        buildCache();
        return getSampleNames().containsAll(samples);
    }

    public GenotypeCollection subsetToSamples( final Collection<String> samples ) {
        return subsetToSamples(new HashSet<String>(samples));
    }

    public GenotypeCollection subsetToSamples( final Set<String> samples ) {
        if ( samples.size() == genotypes.size() )
            return this;
        else if ( samples.isEmpty() )
            return NO_GENOTYPES;
        else {
            GenotypeCollection subset = create(samples.size());
            for ( final Genotype g : genotypes )
                if ( samples.contains(g.getSampleName()) )
                    subset.add(g);
            return subset;
        }
    }
}

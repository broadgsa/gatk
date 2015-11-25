/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.commandline;

import com.google.java.contract.Ensures;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.*;

/**
 * A RodBindingCollection represents a collection of RodBindings.
 *
 * The RodBindingCollection<T> is a formal GATK argument that is used to specify a file of RodBindings.
 *
 */
public final class RodBindingCollection<T extends Feature> {

    /** The Java class expected for this RodBinding.  Must correspond to the type emitted by Tribble */
    final private Class<T> type;

    private Collection<RodBinding<T>> rodBindings;

    public RodBindingCollection(final Class<T> type, final Collection<RodBinding<T>> rodBindings) {
        this.type = type;
        this.rodBindings = Collections.unmodifiableCollection(rodBindings);
    }

    /**
     * @return the collection of RodBindings
     */
    final public Collection<RodBinding<T>> getRodBindings() {
        return rodBindings;
    }

    /**
     * @return the string name of the tribble type, such as vcf, bed, etc.
     */
    @Ensures({"result != null"})
    final public Class<T> getType() {
        return type;
    }

    @Override
    public String toString() {
        return String.format("(RodBindingCollection %s)", getRodBindings());
    }

    /**
     * Utility method to help construct a RodBindingCollection of the given Feature type
     *
     * @param type         the Feature type
     * @param rodBindings  the rod bindings to put into the collection
     * @return a new RodBindingCollection object
     */
    public static Object createRodBindingCollectionOfType(final Class<? extends Feature> type, final Collection<RodBinding> rodBindings) {
        try {
            final Constructor ctor = RodBindingCollection.class.getConstructor(Class.class, Collection.class);
            return ctor.newInstance(type, rodBindings);
        } catch (final Exception e) {
            throw new IllegalStateException("Failed to create a RodBindingCollection for type " + type);
        }
    }
}

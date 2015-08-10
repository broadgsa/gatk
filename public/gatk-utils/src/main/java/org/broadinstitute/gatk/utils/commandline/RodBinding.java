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
import com.google.java.contract.Requires;
import htsjdk.tribble.Feature;

import java.util.*;

/**
 * A RodBinding represents a walker argument that gets bound to a ROD track.
 *
 * The RodBinding<T> is a formal GATK argument that bridges between a walker and
 * the RefMetaDataTracker to obtain data about this rod track at runtime.  The RodBinding
 * is explicitly typed with type of the Tribble.Feature expected to be produced by this
 * argument.  The GATK Engine takes care of initializing the binding and connecting it
 * to the RMD system.
 *
 * It is recommended that optional RodBindings be initialized to the value returned
 * by the static method makeUnbound().
 *
 * Note that this class is immutable.
 */
public final class RodBinding<T extends Feature> {
    protected final static String UNBOUND_VARIABLE_NAME = "";
    protected final static String UNBOUND_SOURCE = "UNBOUND";
    protected final static String UNBOUND_TRIBBLE_TYPE = "";

    /**
     * Create an unbound Rodbinding of type.  This is the correct programming
     * style for an optional RodBinding<T>
     *
     *     At Input()
     *     RodBinding<T> x = RodBinding.makeUnbound(T.class)
     *
     * The unbound binding is guaranteed to never match any binding.  It uniquely
     * returns false to isBound().
     *
     * @param type the Class type produced by this unbound object
     * @param <T> any class extending Tribble Feature
     * @return the UNBOUND RodBinding producing objects of type T
     */
    @Requires("type != null")
    protected final static <T extends Feature> RodBinding<T> makeUnbound(Class<T> type) {
        return new RodBinding<T>(type);
    }

    /** The name of this binding.  Often the name of the field itself, but can be overridden on cmdline */
    final private String name;
    /** where the data for this ROD is coming from.  A file or special value if coming from stdin */
    final private String source;
    /** the string name of the tribble type, such as vcf, bed, etc. */
    final private String tribbleType;
    /** The command line tags associated with this RodBinding */
    final private Tags tags;
    /** The Java class expected for this RodBinding.  Must correspond to the type emitted by Tribble */
    final private Class<T> type;
    /** True for all RodBindings except the special UNBOUND binding, which is the default for optional arguments */
    final private boolean bound;

    /**
     * The name counter.  This is how we create unique names for collections of RodBindings
     * on the command line.  If you have provide the GATK with -X file1 and -X file2 to a
     * RodBinding argument as List<RodBinding<T>> then each binding will receive automatically
     * the name of X and X2.
     */
    final private static Map<String, Integer> nameCounter = new HashMap<String, Integer>();

    /** for UnitTests */
    final public static void resetNameCounter() {
        nameCounter.clear();
    }

    @Requires("rawName != null")
    @Ensures("result != null")
    final private static synchronized String countedVariableName(final String rawName) {
        Integer count = nameCounter.get(rawName);
        if ( count == null ) {
            nameCounter.put(rawName, 1);
            return rawName;
        } else {
            nameCounter.put(rawName, count + 1);
            return rawName + (count + 1);
        }
    }

    @Requires({"type != null", "rawName != null", "source != null", "tribbleType != null", "tags != null"})
    public RodBinding(Class<T> type, final String rawName, final String source, final String tribbleType, final Tags tags) {
        this.type = type;
        this.name = countedVariableName(rawName);
        this.source = source;
        this.tribbleType = tribbleType;
        this.tags = tags;
        this.bound = true;
    }

    /**
     * For testing purposes only.  Creates a RodBinding sufficient for looking up associations to rawName
     * @param type
     * @param rawName
     */
    public RodBinding(Class<T> type, final String rawName) {
        this(type, rawName, "missing", type.getSimpleName(), new Tags());
    }

    /**
     * Make an unbound RodBinding<T>.  Only available for creating the globally unique UNBOUND object
     * @param type class this unbound RodBinding creates
     */
    @Requires({"type != null"})
    private RodBinding(Class<T> type) {
        this.type = type;
        this.name = UNBOUND_VARIABLE_NAME;  // special value can never be found in RefMetaDataTracker
        this.source = UNBOUND_SOURCE;
        this.tribbleType = UNBOUND_TRIBBLE_TYPE;
        this.tags = new Tags();
        this.bound = false;
    }


   /**
     * @return True for all RodBindings except the special UNBOUND binding, which is the default for optional arguments
     */
    final public boolean isBound() {
        return bound;
    }

    /**
     * @return The name of this binding.  Often the name of the field itself, but can be overridden on cmdline
     */
    @Ensures({"result != null"})
    final public String getName() {
        return name;
    }

    /**
     * @return the string name of the tribble type, such as vcf, bed, etc.
     */
    @Ensures({"result != null"})
    final public Class<T> getType() {
        return type;
    }

    /**
     * @return where the data for this ROD is coming from.  A file or special value if coming from stdin
     */
    @Ensures({"result != null"})
    final public String getSource() {
        return source;
    }

    /**
     * @return The command line tags associated with this RodBinding.  Will include the tags used to
     * determine the name and type of this RodBinding
     */
    @Ensures({"result != null"})
    final public Tags getTags() {
        return tags;
    }

    /**
     * @return The Java class expected for this RodBinding.  Must correspond to the type emited by Tribble
     */
    @Ensures({"result != null"})
    final public String getTribbleType() {
        return tribbleType;
    }

    @Override
    public String toString() {
        return String.format("(RodBinding name=%s source=%s)", getName(), getSource());
    }
}

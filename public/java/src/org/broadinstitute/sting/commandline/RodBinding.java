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

package org.broadinstitute.sting.commandline;

import org.broad.tribble.Feature;

import java.util.*;

/**
 * A RodBinding representing a walker argument that gets bound to a ROD track.
 *
 * There is no constraint on the type of the ROD bound.
 */
public class RodBinding<T extends Feature> {
    protected final static String UNBOUND_VARIABLE_NAME = "";
    protected final static String UNBOUND_SOURCE = "UNBOUND";
    protected final static String UNBOUND_TRIBBLE_TYPE = null;
    public final static <T extends Feature> RodBinding<T> makeUnbound(Class<T> type) {
        return new RodBinding<T>(type);
    }

    final private String name;
    final private String source;
    final private String tribbleType;
    final private Tags tags;
    final private Class<T> type;
    final private boolean bound;

    final private static Map<String, Integer> nameCounter = new HashMap<String, Integer>();

    final protected static void resetNameCounter() {
        nameCounter.clear();
    }

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

    public boolean isBound() {
        return bound;
    }

    public RodBinding(Class<T> type, final String rawName, final String source, final String tribbleType, final Tags tags) {
        this.type = type;
        this.name = countedVariableName(rawName);
        this.source = source;
        this.tribbleType = tribbleType;
        this.tags = tags;
        this.bound = true;
    }

    /**
     * Make an unbound RodBinding<T>
     * @param type
     */
    private RodBinding(Class<T> type) {
        this.type = type;
        this.name = UNBOUND_VARIABLE_NAME;  // special value can never be found in RefMetaDataTracker
        this.source = UNBOUND_SOURCE;
        this.tribbleType = UNBOUND_TRIBBLE_TYPE;
        this.tags = new Tags();
        this.bound = false;
    }

    public String getName() {
        return name;
    }
    public Class<T> getType() {
        return type;
    }
    public String getSource() {
        return source;
    }

    public Tags getTags() {
        return tags;
    }

    public String getTribbleType() {
        return tribbleType;
    }

    public String toString() {
        return String.format("(RodBinding name=%s source=%s)", getName(), getSource());
    }
}

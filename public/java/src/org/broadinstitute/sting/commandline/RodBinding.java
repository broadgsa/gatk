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
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;

import java.util.List;

/**
 * A RodBinding representing a walker argument that gets bound to a ROD track.
 *
 * There is no constraint on the type of the ROD bound.
 */
public class RodBinding<T extends Feature> {
    final String variableName;
    final String source;
    final ParsingEngine parser;
    final Class<T> type;

    protected RodBinding(Class<T> type, final String variableName, final String source, final ParsingEngine parser) {
        this.type = type;
        this.variableName = variableName;
        this.source = source;
        this.parser = parser;
    }

    public String getVariableName() {
        return variableName;
    }

    public String getSource() {
        return source;
    }

    public List<T> getValues(final RefMetaDataTracker tracker) {
        return tracker.getValues(variableName, type);
    }

//    public <T> List<T> getValues(final RefMetaDataTracker tracker, final Class<T> clazz) {
//        return tracker.getValues(variableName, clazz);
//    }

    public T getFirstValue(final RefMetaDataTracker tracker) {
        return tracker.getFirstValue(variableName, type);
    }

    public boolean hasValues(final RefMetaDataTracker tracker) {
        return tracker.hasValues(variableName);
    }

    public List<GATKFeature> getValuesAsGATKFeatures(final RefMetaDataTracker tracker) {
        return tracker.getValuesAsGATKFeatures(variableName);
    }

    public Tags getTags() {
        return parser.getTags(this);
    }

    public String toString() {
        return String.format("(RodBinding name=%s source=%s)", getVariableName(), getSource());
    }
}

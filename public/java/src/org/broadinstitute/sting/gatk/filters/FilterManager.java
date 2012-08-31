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

package org.broadinstitute.sting.gatk.filters;

import com.google.common.base.Function;
import com.google.common.collect.Collections2;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;

import java.util.Collection;
import java.util.List;

/**
 * Manage filters and filter options.  Any requests for basic filtering classes
 * should ultimately be made through this class.
 *
 * @author mhanna
 * @version 0.1
 */
public class FilterManager extends PluginManager<ReadFilter> {
    public FilterManager() {
        super(ReadFilter.class,"filter","Filter");
    }

    /**
     * Instantiate a filter of the given type.  Along the way, scream bloody murder if
     * the filter is not available.
     * @param filterType The type of the filter
     * @return The filter
     */
    public ReadFilter createFilterByType(Class<? extends ReadFilter> filterType) {
        return this.createByName(getName(filterType));
    }

    public Collection<Class<? extends ReadFilter>> getValues() {
        return this.getPlugins();
    }

    /**
     * Rather than use the default error message, print out a list of read filters as well.
     * @param pluginCategory - string, the category of the plugin (e.g. read filter)
     * @param pluginName - string, what we were trying to match (but failed to)
     * @return - A wall of text with the default message, followed by a listing of available read filters
     */
    @Override
    protected String formatErrorMessage(String pluginCategory, String pluginName) {
        List<Class<? extends ReadFilter>> availableFilters = this.getPluginsImplementing(ReadFilter.class);
        Collection<String> availableFilterNames = Collections2.transform(availableFilters, new Function<Class<? extends ReadFilter>,String>(){

            @Override
            public String apply(final Class<? extends ReadFilter> input) {
                return getName(input);
            }
        });

        return String.format("Read filter %s not found. Available read filters:%n%s.%n%n%s",pluginName,
                Utils.join(String.format(", "),availableFilterNames),
                "Please consult the GATK Documentation (http://www.broadinstitute.org/gatk/gatkdocs/) for more information.");
    }
}

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

package org.broadinstitute.gatk.engine.filters;

import org.broadinstitute.gatk.utils.classloader.PluginManager;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.GATKDocUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;

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


        return String.format("Read filter %s not found. Available read filters:%n%n%s%n%n%s",pluginName,
                userFriendlyListofReadFilters(availableFilters),
                "Please consult the GATK Documentation (" + HelpConstants.GATK_DOCS_URL + ") for more information.");
    }

    /**
     * Rather than use the default exception, return a MalformedReadFilterException.
     * @param errorMessage error message from formatErrorMessage()
     * @return - A MalformedReadFilterException with errorMessage
     */
    @Override
    protected UserException createMalformedArgumentException(final String errorMessage) {
        return new UserException.MalformedReadFilterException(errorMessage);
    }

    private String userFriendlyListofReadFilters(List<Class<? extends ReadFilter>> filters) {
        final String headName = "FilterName", headDoc = "Documentation";
        int longestNameLength = -1;
        for ( Class < ? extends ReadFilter> filter : filters ) {
            longestNameLength = Math.max(longestNameLength,this.getName(filter).length());
        }
        String format = "   %"+longestNameLength+"s        %s%n";

        StringBuilder listBuilder = new StringBuilder();
        listBuilder.append(String.format(format,headName,headDoc));
        for ( Class<? extends ReadFilter> filter : filters ) {
            String helpLink = GATKDocUtils.helpLinksToGATKDocs(filter);
            String filterName = this.getName(filter);
            listBuilder.append(String.format(format,filterName,helpLink));
        }

        return listBuilder.toString();
    }
}

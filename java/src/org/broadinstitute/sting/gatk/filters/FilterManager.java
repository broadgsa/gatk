package org.broadinstitute.sting.gatk.filters;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.PluginManager;

import net.sf.picard.filter.SamRecordFilter;

/**
 * Manage filters and filter options.  Any requests for basic filtering classes
 * should ultimately be made through this class.
 *
 * @author mhanna
 * @version 0.1
 */
public class FilterManager extends PluginManager<SamRecordFilter> {
    /**
     * our log, which we want to capture anything from this class
     */
    private static Logger logger = Logger.getLogger(FilterManager.class);

    public FilterManager() {
        super(SamRecordFilter.class,"filter","Filter");
    }

    /**
     * Instantiate a filter of the given type.  Along the way, scream bloody murder if
     * the filter is not available.
     * @param filterType
     * @return
     */
    public SamRecordFilter createFilterByType(Class<? extends SamRecordFilter> filterType) {
        return this.createByName(getName(filterType));
    }
}

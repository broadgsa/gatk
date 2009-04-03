package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.List;
import java.util.Collection;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Apr 3, 2009
 * Time: 3:05:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class RefMetaDataTracker {
    final HashMap<String, ReferenceOrderedDatum> map = new HashMap<String, ReferenceOrderedDatum>();
    protected static Logger logger = Logger.getLogger(RefMetaDataTracker.class);

    /**
     * Finds the reference meta data named name, if it exists, otherwise returns the defaultValue
     *
     * @param name
     * @param defaultValue
     * @return
     */
    public ReferenceOrderedDatum lookup(final String name, ReferenceOrderedDatum defaultValue) {
        //logger.debug(String.format("Lookup %s%n", name));
        if ( map.containsKey(name) )
            return map.get(name);
        else
            return defaultValue;
    }

    public Object lookup(final String name, Object defaultValue) {
        if ( map.containsKey(name) )
            return map.get(name);
        else
            return defaultValue;
    }

    public Object hasROD(final String name) {
       return map.containsKey(name);
    }

    public Collection<ReferenceOrderedDatum> getAllRods() {
        return map.values();
    }

    public void bind(final String name, ReferenceOrderedDatum rod) {
        //logger.debug(String.format("Binding %s to %s%n", name, rod));
        map.put(name, rod);
    }
}

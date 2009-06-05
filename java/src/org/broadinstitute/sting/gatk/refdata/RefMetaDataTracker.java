package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.List;
import java.util.Collection;

/**
 * This class represents the Reference Metadata available at a particular site in the genome.  It can be
 * used to conveniently lookup the RODs at this site, as well just getting a list of all of the RODs
 *
 * The standard interaction model is:
 *
 * Traversal system arrives at a site, which has a bunch of rods covering it
 * Traversal calls tracker.bind(name, rod) for each rod in rods
 * Traversal passes tracker to the walker
 * walker calls lookup(name, default) to obtain the rod values at this site, or default if none was
 *   bound at this site.
 *
 * User: mdepristo
 * Date: Apr 3, 2009
 * Time: 3:05:23 PM
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
        final String luName = canonicalName(name);
        if ( map.containsKey(luName) )
            return map.get(luName);
        else
            return defaultValue;
    }

    /**
     * @see this.lookup
     * @param name
     * @param defaultValue
     * @return
     */
    public Object lookup(final String name, Object defaultValue) {
        final String luName = canonicalName(name);
        if ( map.containsKey(luName) )
            return map.get(luName);
        else
            return defaultValue;
    }

    /**
     * Returns the canonical name of the rod name
     * @param name
     * @return
     */
    private final String canonicalName(final String name)
    {
        //return name; // .toLowerCase();
        return name.toLowerCase();
    }

    /**
     * Is there a binding at this site to a ROD with name?
     *
     * @param name
     * @return
     */
    public Object hasROD(final String name) {
       return map.containsKey(canonicalName(name));
    }

    /**
     * Get all of the RODs at the current site
     * 
     * @return
     */
    public Collection<ReferenceOrderedDatum> getAllRods() {
        return map.values();
    }

    /**
     * Binds the reference ordered datum ROD to name at this site.  Should be used only but the traversal
     * system to provide access to RODs in a structured way to the walkers.
     *
     * @param name
     * @param rod
     */
    public void bind(final String name, ReferenceOrderedDatum rod) {
        //logger.debug(String.format("Binding %s to %s", name, rod));
        map.put(canonicalName(name), rod);
    }
}

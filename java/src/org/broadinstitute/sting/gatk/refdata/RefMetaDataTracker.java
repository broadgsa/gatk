package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;

import java.util.*;

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
    final HashMap<String, RODRecordList<ReferenceOrderedDatum>> map = new HashMap<String, RODRecordList<ReferenceOrderedDatum>>();
    protected static Logger logger = Logger.getLogger(RefMetaDataTracker.class);

    /**
     * Finds the reference meta data named name, if it exists, otherwise returns the defaultValue.
     * This is a legacy method that works with "singleton" tracks, in which a single ROD record can be associated
     * with any given site. If track provides multiple records associated with a site, this method will return
     * the first one.
     * @param name
     * @param defaultValue
     * @return
     */
    @Deprecated
    public ReferenceOrderedDatum lookup(final String name, ReferenceOrderedDatum defaultValue) {
        //logger.debug(String.format("Lookup %s%n", name));
        final String luName = canonicalName(name);
        if ( map.containsKey(luName) ) {
            RODRecordList<ReferenceOrderedDatum> value = map.get(luName) ;
            if ( value != null ) {
                List<ReferenceOrderedDatum> l = value.getRecords();
                if ( l != null & l.size() > 0 ) return value.getRecords().get(0);
            }
        } 
        return defaultValue;
    }

    /**
     * Finds the reference metadata track named 'name' and returns all ROD records from that track associated
     * with the current site as a RODRecordList collection object. If no data track with specified name is available,
     * returns defaultValue wrapped as RODRecordList object. NOTE: if defaultValue is null, it will be wrapped up
     * with track name set to 'name' and location set to null; otherwise the wrapper object will have name and
     * location set to defaultValue.getName() and defaultValue.getLocation(), respectively (use caution,
     * defaultValue.getLocation() may be not equal to what RODRecordList's location would be expected to be otherwise:
     * for instance, on locus traversal, location is usually expected to be a single base we are currently looking at,
     * regardless of the presence of "extended" RODs overlapping with that location).
     * @param name
     * @param defaultValue
     * @return
     */
    public RODRecordList<ReferenceOrderedDatum> getTrackData(final String name, ReferenceOrderedDatum defaultValue) {
        //logger.debug(String.format("Lookup %s%n", name));
        final String luName = canonicalName(name);
        if ( map.containsKey(luName) )
            return map.get(luName);
        else {

            if ( defaultValue == null )
                return null;
            return new RODRecordList<ReferenceOrderedDatum>(defaultValue.getName(),
                                                            Collections.singletonList(defaultValue),
                                                            defaultValue.getLocation());
        }
    }
    /**
     * @see this.lookup
     * @param name
     * @param defaultValue
     * @return
     */
    @Deprecated
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
     * Is there a binding at this site to a ROD/track with the specified name?
     *
     * @param name the name of the rod
     * @return true if it has the rod
     */
    public boolean hasROD(final String name) {
       return map.containsKey(canonicalName(name));
    }

    /**
     * Get all of the RODs at the current site. The collection is "flattened": for any track that has multiple records
     * at the current site, they all will be added to the list as separate elements.
     * 
     * @return
     */
    public Collection<ReferenceOrderedDatum> getAllRods() {
        List<ReferenceOrderedDatum> l = new ArrayList<ReferenceOrderedDatum>();
        for ( RODRecordList<ReferenceOrderedDatum> rl : map.values() ) {
            if ( rl == null ) continue; // how do we get null value stored for a track? shouldn't the track be missing from the map alltogether?
            l.addAll(rl.getRecords());
        }
        return l;

    }

    /**
     * Get all of the ROD tracks at the current site. Each track is returned as a single compound
     * object (RODRecordList) that may contain multiple ROD records associated with the current site.
     *
     * @return
     */
    public Collection<RODRecordList<ReferenceOrderedDatum>> getBoundRodTracks() {
        LinkedList<RODRecordList<ReferenceOrderedDatum>> bound = new LinkedList<RODRecordList<ReferenceOrderedDatum>>();
        
        for ( RODRecordList<ReferenceOrderedDatum> value : map.values() ) {
             if ( value != null && value.size() != 0 ) bound.add(value);
        }

        return bound;
    }

    public int getNBoundRodTracks() {
        return getNBoundRodTracks(null);
    }

    public int getNBoundRodTracks(final String excludeIn ) {
        final String exclude = excludeIn == null ? null : canonicalName(excludeIn);

        int n = 0;
        for ( RODRecordList<ReferenceOrderedDatum> value : map.values() ) {
             if ( value != null && ! value.isEmpty() ) {
                 if ( exclude == null || ! value.getName().equals(exclude) )
                    n++;
             }
        }

        return n;
    }

    public Collection<ReferenceOrderedDatum> getBoundRodRecords() {
        LinkedList<ReferenceOrderedDatum> bound = new LinkedList<ReferenceOrderedDatum>();

        for ( RODRecordList<ReferenceOrderedDatum> valueList : map.values() ) {
            for ( ReferenceOrderedDatum value : valueList ) {
                if ( value != null )
                bound.add(value);
            }
        }

        return bound;
    }
    /**
     * Binds the list of reference ordered data records (RODs) to track name at this site.  Should be used only by the traversal
     * system to provide access to RODs in a structured way to the walkers.
     *
     * @param name
     * @param rod
     */
    public void bind(final String name, RODRecordList<ReferenceOrderedDatum>  rod) {
        //logger.debug(String.format("Binding %s to %s", name, rod));
        map.put(canonicalName(name), rod);
    }
/*
    public void bind(final String name, ReferenceOrderedDatum rod) {
        //logger.debug(String.format("Binding %s to %s", name, rod));
        map.put(canonicalName(name), rod);
    }
    */
}

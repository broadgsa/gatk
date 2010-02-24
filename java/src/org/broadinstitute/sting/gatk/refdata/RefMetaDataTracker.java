package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * This class represents the Reference Metadata available at a particular site in the genome.  It can be
 * used to conveniently lookup the RODs at this site, as well just getting a list of all of the RODs
 *
 * The standard interaction model is:
 *
 * Traversal system arrives at a site, which has a bunch of rods covering it
Genotype * Traversal calls tracker.bind(name, rod) for each rod in rods
 * Traversal passes tracker to the walker
 * walker calls lookup(name, default) to obtain the rod values at this site, or default if none was
 *   bound at this site.
 *
 * User: mdepristo
 * Date: Apr 3, 2009
 * Time: 3:05:23 PM
 */
public class RefMetaDataTracker {
    final HashMap<String, RODRecordList> map = new HashMap<String, RODRecordList>();
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
            RODRecordList value = map.get(luName) ;
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
    public RODRecordList getTrackData(final String name, ReferenceOrderedDatum defaultValue, boolean requireExactMatch) {
        //logger.debug(String.format("Lookup %s%n", name));

        final String luName = canonicalName(name);
        RODRecordList trackData = null;

        if ( requireExactMatch ) {
            if ( map.containsKey(luName) )
                trackData = map.get(luName);
        } else {
            for ( Map.Entry<String, RODRecordList> datum : map.entrySet() ) {
                final String rodName = datum.getKey();
                if ( rodName.startsWith(luName) ) {
                    if ( trackData == null ) trackData = new RODRecordList(name);
                    //System.out.printf("Adding bindings from %s to %s at %s%n", rodName, name, datum.getValue().getLocation());
                    trackData.add(datum.getValue(), true);
                }
            }
        }

        if ( trackData != null )
            return trackData;
        else if ( defaultValue == null )
            return null;
        else
            return new RODRecordList(defaultValue.getName(),
                    Collections.singletonList(defaultValue),
                    defaultValue.getLocation());
    }

    public RODRecordList getTrackData(final String name, ReferenceOrderedDatum defaultValue) {
        return getTrackData(name, defaultValue, true);
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
        for ( RODRecordList rl : map.values() ) {
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
    public Collection<RODRecordList> getBoundRodTracks() {
        LinkedList<RODRecordList> bound = new LinkedList<RODRecordList>();
        
        for ( RODRecordList value : map.values() ) {
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
        for ( RODRecordList value : map.values() ) {
             if ( value != null && ! value.isEmpty() ) {
                 if ( exclude == null || ! value.getName().equals(exclude) )
                    n++;
             }
        }

        return n;
    }

    public Collection<ReferenceOrderedDatum> getBoundRodRecords() {
        LinkedList<ReferenceOrderedDatum> bound = new LinkedList<ReferenceOrderedDatum>();

        for ( RODRecordList valueList : map.values() ) {
            for ( ReferenceOrderedDatum value : valueList ) {
                if ( value != null )
                bound.add(value);
            }
        }

        return bound;
    }


    /**
     * Converts all possible ROD tracks to VariantContexts objects, of all types, allowing any start and any number
     * of entries per ROD.
     */
    public Collection<VariantContext> getAllVariantContexts() {
        return getAllVariantContexts(null, null, false, false);
    }


    /**
     * Converts all possible ROD tracks to VariantContexts objects.  If allowedTypes != null, then only
     * VariantContexts in the allow set of types will be returned.  If requireStartsHere is true, then curLocation
     * must not be null, and only records whose start position is == to curLocation.getStart() will be returned.
     * If takeFirstOnly is true, then only a single VariantContext will be converted from any individual ROD.  Of course,
     * this single object must pass the allowed types and start here options if provided.  Note that the result
     * may return multiple VariantContexts with the same name if that particular track contained multiple RODs spanning
     * the current location.
     *
     * The name of each VariantContext corresponds to the ROD name.
     *
     * @param curLocation
     * @param allowedTypes
     * @param requireStartHere
     * @param takeFirstOnly
     * @return
     */
    public Collection<VariantContext> getAllVariantContexts(EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        List<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( RODRecordList rodList : getBoundRodTracks() ) {
            addVariantContexts(contexts, rodList, allowedTypes, curLocation, requireStartHere, takeFirstOnly);
        }

        return contexts;
    }

    /**
     * Gets the variant contexts associated with track name name
     *
     * see getVariantContexts for more information.
     *
     * @param name
     * @param curLocation
     * @param allowedTypes
     * @param requireStartHere
     * @param takeFirstOnly
     * @return
     */
    public Collection<VariantContext> getVariantContexts(String name, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        return getVariantContexts(Arrays.asList(name), allowedTypes, curLocation, requireStartHere, takeFirstOnly);
    }

    public Collection<VariantContext> getVariantContexts(Collection<String> names, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        Collection<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( String name : names ) {
            RODRecordList rodList = getTrackData(name, null);

            if ( rodList != null )
                addVariantContexts(contexts, rodList, allowedTypes, curLocation, requireStartHere, takeFirstOnly );
        }

        return contexts;
    }

    /**
     * Gets the variant context associated with name, and assumes the system only has a single bound track at this location.  Throws an exception if not.
     * see getVariantContexts for more information.
     *
     * @param name
     * @param curLocation
     * @param allowedTypes
     * @param requireStartHere
     * @return
     */
    public VariantContext getVariantContext(String name, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere ) {
        Collection<VariantContext> contexts = getVariantContexts(name, allowedTypes, curLocation, requireStartHere, false );

        if ( contexts.size() > 1 )
            throw new StingException("Requested a single VariantContext object for track " + name + " but multiple variants were present at position " + curLocation);
        else if ( contexts.size() == 0 )
            return null;
        else
            return contexts.iterator().next();
    }

    private void addVariantContexts(Collection<VariantContext> contexts, RODRecordList rodList, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        for ( ReferenceOrderedDatum rec : rodList.getRecords() ) {
            if ( VariantContextAdaptors.canBeConvertedToVariantContext(rec) ) {
                // ok, we might actually be able to turn this record in a variant context
                VariantContext vc = VariantContextAdaptors.toVariantContext(rodList.getName(), rec);

                if ( vc == null ) // sometimes the track has odd stuff in it that can't be converted 
                    continue;

                // now, let's decide if we want to keep it
                boolean goodType = allowedTypes == null || allowedTypes.contains(vc.getType());
                boolean goodPos = ! requireStartHere || rec.getLocation().getStart() == curLocation.getStart();

                if ( goodType && goodPos ) {  // ok, we are going to keep this thing
                    contexts.add(vc);

                    if ( takeFirstOnly )
                        // we only want the first passing instance, so break the loop over records in rodList
                        break;
                }
            }
        }
    }


    /**
     * Binds the list of reference ordered data records (RODs) to track name at this site.  Should be used only by the traversal
     * system to provide access to RODs in a structured way to the walkers.
     *
     * @param name
     * @param rod
     */
    public void bind(final String name, RODRecordList rod) {
        //logger.debug(String.format("Binding %s to %s", name, rod));
        map.put(canonicalName(name), rod);
    }
}

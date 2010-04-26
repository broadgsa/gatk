package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * This class represents the Reference Metadata available at a particular site in the genome.  It can be
 * used to conveniently lookup the RMDs at this site, as well just getting a list of all of the RMDs
 *
 * The standard interaction model is:
 *
 * Traversal system arrives at a site, which has a bunch of RMDs covering it
Genotype * Traversal calls tracker.bind(name, RMD) for each RMDs in RMDs
 * Traversal passes tracker to the walker
 * walker calls lookup(name, default) to obtain the RMDs values at this site, or default if none was
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
     * get all the reference meta data associated with a track name.
     * @param name the name of the track we're looking for
     * @return a list of objects, representing the underlying objects that the tracks produce.  I.e. for a
     *         dbSNP RMD this will be a RodDbSNP, etc.
     *
     * Important: The list returned by this function is guaranteed not to be null, but may be empty!
     */
    public List<Object> getReferenceMetaData(final String name) {
        RODRecordList list = getTrackDataByName(name, true);
        List<Object> objects = new ArrayList<Object>();
        if (list == null) return objects;
        for (GATKFeature feature : list)
            objects.add(feature.getUnderlyingObject());
        return objects;
    }

    /**
     * get all the reference meta data associated with a track name.
     * @param name the name of the track we're looking for
     * @param requireExactMatch do we require an exact match for the name (true) or do we require only that the name starts with
     *        the passed in parameter (false).
     * @return a list of objects, representing the underlying objects that the tracks produce.  I.e. for a
     *         dbSNP rod this will be a RodDbSNP, etc.
     *
     * Important: The list returned by this function is guaranteed not to be null, but may be empty!
     */
    public List<Object> getReferenceMetaData(final String name, boolean requireExactMatch) {
        RODRecordList list = getTrackDataByName(name, requireExactMatch);
        List<Object> objects = new ArrayList<Object>();
        if (list == null) return objects;
        for (GATKFeature feature : list)
            objects.add(feature.getUnderlyingObject());
        return objects;
    }

    /**
     * get a singleton record, given the name and a type.  This function will return the first record at the current position seen,
     * and emit a logger warning if there were more than one option.
     *
     * WARNING: this method is deprecated, since we now suppport more than one RMD at a single position for all tracks.  If there are
     * are multiple RMD objects at this location, there is no contract for which object this method will pick, and which object gets
     * picked may change from time to time!  BE WARNED!
     * 
     * @param name the name of the track
     * @param clazz the underlying type to return
     * @param <T> the type to parameterize on, matching the clazz argument
     * @return a record of type T, or null if no record is present.
     */
    @Deprecated
    public <T> T lookup(final String name, Class<T> clazz) {
        RODRecordList objects = getTrackDataByName(name, true);

        // if emtpy or null return null;
        if (objects == null || objects.size() < 1) return null;

        if (objects.size() > 1)
            logger.info("lookup is choosing the first record from " + (objects.size() - 1) + " options");

        Object obj = objects.get(0).getUnderlyingObject();
        if (!(clazz.isAssignableFrom(obj.getClass())))
            throw new StingException("Unable to case track named " + name + " to type of " + clazz.toString()
                                     + " it's of type " + obj.getClass());

        return (T)obj;
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
     * Get all of the RMDs at the current site. The collection is "flattened": for any track that has multiple records
     * at the current site, they all will be added to the list as separate elements.
     *
     * @return collection of all rods
     */
    public Collection<GATKFeature> getAllRods() {
        List<GATKFeature> l = new ArrayList<GATKFeature>();
        for ( RODRecordList rl : map.values() ) {
            if ( rl == null ) continue; // how do we get null value stored for a track? shouldn't the track be missing from the map alltogether?
            l.addAll(rl);
        }
        return l;

    }

    /**
     * Get all of the RMD tracks at the current site. Each track is returned as a single compound
     * object (RODRecordList) that may contain multiple RMD records associated with the current site.
     *
     * @return collection of all tracks
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


    /**
     * Binds the list of reference ordered data records (RMDs) to track name at this site.  Should be used only by the traversal
     * system to provide access to RMDs in a structured way to the walkers.
     *
     * @param name the name of the track
     * @param rod the collection of RMD data
     */
    public void bind(final String name, RODRecordList rod) {
        //logger.debug(String.format("Binding %s to %s", name, rod));
        map.put(canonicalName(name), rod);
    }


    /**
     * Converts all possible ROD tracks to VariantContexts objects, of all types, allowing any start and any number
     * of entries per ROD.
     * The name of each VariantContext corresponds to the ROD name.
     *
     * @param ref                reference context
     * @return variant context
     */
    public Collection<VariantContext> getAllVariantContexts(ReferenceContext ref) {
        return getAllVariantContexts(ref, null, null, false, false);
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
     * @param ref                reference context
     * @param allowedTypes       allowed types
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @param takeFirstOnly      do we take the first rod only?
     * @return variant context
     */
    public Collection<VariantContext> getAllVariantContexts(ReferenceContext ref, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        List<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( RODRecordList rodList : getBoundRodTracks() ) {
            addVariantContexts(contexts, rodList, ref, allowedTypes, curLocation, requireStartHere, takeFirstOnly);
        }

        return contexts;
    }

    /**
     * Gets the variant contexts associated with track name name
     *
     * see getVariantContexts for more information.
     *
     * @param name               name
     * @param curLocation        location
     * @param allowedTypes       allowed types
     * @param requireStartHere   do we require the rod to start at this location?
     * @param takeFirstOnly      do we take the first rod only?
     * @return variant context
     */
    public Collection<VariantContext> getVariantContexts(String name, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        return getVariantContexts(null, Arrays.asList(name), allowedTypes, curLocation, requireStartHere, takeFirstOnly);
    }

    public Collection<VariantContext> getVariantContexts(ReferenceContext ref, String name, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        return getVariantContexts(ref, Arrays.asList(name), allowedTypes, curLocation, requireStartHere, takeFirstOnly);
    }

    public Collection<VariantContext> getVariantContexts(Collection<String> names, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        return getVariantContexts(null, names, allowedTypes, curLocation, requireStartHere, takeFirstOnly);
    }

    public Collection<VariantContext> getVariantContexts(ReferenceContext ref, Collection<String> names, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        Collection<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( String name : names ) {
            RODRecordList rodList = getTrackDataByName(name,true); // require that the name is an exact match

            if ( rodList != null )
                addVariantContexts(contexts, rodList, ref, allowedTypes, curLocation, requireStartHere, takeFirstOnly );
        }

        return contexts;
    }

    /**
     * Gets the variant context associated with name, and assumes the system only has a single bound track at this location.  Throws an exception if not.
     * see getVariantContexts for more information.
     *
     * @param name               name
     * @param curLocation        location
     * @param allowedTypes       allowed types
     * @param requireStartHere   do we require the rod to start at this location?
     * @return variant context
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


    private void addVariantContexts(Collection<VariantContext> contexts, RODRecordList rodList, ReferenceContext ref, EnumSet<VariantContext.Type> allowedTypes, GenomeLoc curLocation, boolean requireStartHere, boolean takeFirstOnly ) {
        for ( GATKFeature rec : rodList ) {
            if ( VariantContextAdaptors.canBeConvertedToVariantContext(rec.getUnderlyingObject()) ) {
                // ok, we might actually be able to turn this record in a variant context
                VariantContext vc = VariantContextAdaptors.toVariantContext(rodList.getName(), rec.getUnderlyingObject(), ref);

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
     * Finds the reference metadata track named 'name' and returns all ROD records from that track associated
     * with the current site as a RODRecordList collection object. If no data track with specified name is available,
     * returns defaultValue wrapped as RODRecordList object. NOTE: if defaultValue is null, it will be wrapped up
     * with track name set to 'name' and location set to null; otherwise the wrapper object will have name and
     * location set to defaultValue.getName() and defaultValue.getLocation(), respectively (use caution,
     * defaultValue.getLocation() may be not equal to what RODRecordList's location would be expected to be otherwise:
     * for instance, on locus traversal, location is usually expected to be a single base we are currently looking at,
     * regardless of the presence of "extended" RODs overlapping with that location).
     * @param name                track name
     * @param requireExactMatch   do we require an exact match of the rod name?
     * @return track data for the given rod
     */
    private RODRecordList getTrackDataByName(final String name, boolean requireExactMatch) {
        //logger.debug(String.format("Lookup %s%n", name));

        final String luName = canonicalName(name);
        RODRecordList trackData = null;

        if ( requireExactMatch ) {
            if ( map.containsKey(luName) )
                trackData = map.get(luName);
        } else {
            for ( Map.Entry<String, RODRecordList> datum : map.entrySet() ) {
                final String rodName = datum.getKey();
                if ( datum.getValue() != null && rodName.startsWith(luName) ) {
                    if ( trackData == null ) trackData = new RODRecordListImpl(name);
                    //System.out.printf("Adding bindings from %s to %s at %s%n", rodName, name, datum.getValue().getLocation());
                    ((RODRecordListImpl)trackData).add(datum.getValue(), true);
                }
            }
        }
        return trackData;
    }

    /**
     * Returns the canonical name of the rod name (lowercases it)
     * @param name the name of the rod
     * @return canonical name of the rod
     */
    private final String canonicalName(final String name) {
        return name.toLowerCase();
    }
}

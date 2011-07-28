package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

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
    final Map<String, RODRecordList> map;
    protected static Logger logger = Logger.getLogger(RefMetaDataTracker.class);

    public RefMetaDataTracker(int nBindings) {
        if ( nBindings == 0 )
            map = Collections.emptyMap();
        else
            map = new HashMap<String, RODRecordList>(nBindings);
    }

    /**
     * No-assumption version of getValues(name, class).  Returns Objects.
     */
    public List<Object> getValues(final String name) {
        return getValues(name, Object.class);
    }

    /**
     * get all the reference meta data associated with a track name.
     * @param name the name of the track we're looking for
     * @param clazz the expected class of the elements bound to rod name
     * @return a list of objects, representing the underlying objects that the tracks produce.  I.e. for a
     *         dbSNP RMD this will be a RodDbSNP, etc.
     *
     * Important: The list returned by this function is guaranteed not to be null, but may be empty!
     */
    public <T> List<T> getValues(final String name, final Class<T> clazz) {
        RODRecordList list = getTrackDataByName(name);

        if (list == null)
            return Collections.emptyList();
        else {
            List<T> objects = new ArrayList<T>();
            for (GATKFeature feature : list) {
                final Object obj = feature.getUnderlyingObject();
                if (!(clazz.isAssignableFrom(obj.getClass())))
                    throw new UserException.CommandLineException("Unable to case track named " + name + " to type of " + clazz.toString()
                            + " it's of type " + obj.getClass());
                objects.add((T)obj);
            }
            return objects;
        }
    }

    /**
     * get a singleton record, given the name and a type.  This function will return the first record at the
     * current position seen.  The object is cast into a type clazz, or thoses an error if this isn't possible.
     *
     * * WARNING: we now suppport more than one RMD at a single position for all tracks.  If there are
     * are multiple RMD objects at this location, there is no contract for which object this method will pick, and which object gets
     * picked may change from time to time!  BE WARNED!
     * 
     * @param name the name of the track
     * @param clazz the underlying type to return
     * @param <T> the type to parameterize on, matching the clazz argument
     * @return a record of type T, or null if no record is present.
     */
    public <T> T getFirstValue(final String name, final Class<T> clazz) {
        RODRecordList objects = getTrackDataByName(name);

        // if empty or null return null;
        if (objects == null || objects.size() < 1) return null;

        Object obj = objects.get(0).getUnderlyingObject();
        if (!(clazz.isAssignableFrom(obj.getClass())))
            throw new UserException.CommandLineException("Unable to case track named " + name + " to type of " + clazz.toString()
                    + " it's of type " + obj.getClass());
        else
            return (T)obj;
    }

    /**
     * Is there a binding at this site to a ROD/track with the specified name?
     *
     * @param name the name of the rod
     * @return true if it has the rod
     */
    public boolean hasValues(final String name) {
        return map.containsKey(canonicalName(name));
    }


    /**
     * Get all of the RMDs at the current site. The collection is "flattened": for any track that has multiple records
     * at the current site, they all will be added to the list as separate elements.
     *
     * @return collection of all rods
     */
    public Collection<GATKFeature> getAllValuesAsGATKFeatures() {
        List<GATKFeature> l = new ArrayList<GATKFeature>();
        for ( RODRecordList rl : map.values() ) {
            if ( rl != null )
                l.addAll(rl);
        }
        return l;
    }

        /**
     * get all the GATK features associated with a specific track name
     * @param name the name of the track we're looking for
     * @return a list of GATKFeatures for the target rmd
     *
     * Important: The list returned by this function is guaranteed not to be null, but may be empty!
     */
    public List<GATKFeature> getValuesAsGATKFeatures(final String name) {
        List<GATKFeature> feat = getTrackDataByName(name);
        return (feat == null) ? new ArrayList<GATKFeature>() : feat; // to satisfy the above requirement that we don't return null
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

    /**
     * The number of tracks with at least one value bound here
     * @return
     */
    public int getNumberOfTracksWithValue() {
        int n = 0;
        for ( RODRecordList value : map.values() ) {
            if ( value != null && ! value.isEmpty() ) {
                n++;
            }
        }
        return n;
    }


    /**
     * Binds the list of reference ordered data records (RMDs) to track name at this site.  Should be used only by the traversal
     * system to provide access to RMDs in a structured way to the walkers.
     *
     * DO NOT USE THIS FUNCTION UNLESS YOU ARE THE GATK ENGINE
     *
     * @param name the name of the track
     * @param rod the collection of RMD data
     */
    public void bind(final String name, RODRecordList rod) {
        //logger.debug(String.format("Binding %s to %s", name, rod));
        map.put(canonicalName(name), rod);
    }

    // ------------------------------------------------------------------------------------------
    //
    //
    // VariantContext helpers
    //
    //
    // ------------------------------------------------------------------------------------------

    /**
     * Converts all possible ROD tracks to VariantContexts objects, of all types, allowing any start and any number
     * of entries per ROD.
     * The name of each VariantContext corresponds to the ROD name.
     *
     * @param ref                reference context
     * @return variant context
     */
    public Collection<VariantContext> getAllVariantContexts(final ReferenceContext ref) {
        return getAllVariantContexts(ref, null, false, false);
    }

    /**
     * Returns all of the variant contexts that start at the current location
     * @param ref
     * @param curLocation
     * @return
     */
    public Collection<VariantContext> getAllVariantContexts(final ReferenceContext ref,
                                                            final GenomeLoc curLocation) {
        return getAllVariantContexts(ref, curLocation, true, false);
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
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @param takeFirstOnly      do we take the first rod only?
     * @return variant context
     */
    public Collection<VariantContext> getAllVariantContexts(final ReferenceContext ref,
                                                            final GenomeLoc curLocation,
                                                            final boolean requireStartHere,
                                                            final boolean takeFirstOnly ) {
        List<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( RODRecordList rodList : getBoundRodTracks() ) {
            addVariantContexts(contexts, rodList, ref, curLocation, requireStartHere, takeFirstOnly);
        }

        return contexts;
    }

    /**
     * Gets the variant contexts associated with track name name
     *
     * see getVariantContexts for more information.
     *
     * @param ref                ReferenceContext to enable conversion to variant context
     * @param name               name
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @param takeFirstOnly      do we take the first rod only?
     * @return variant context
     */
    public Collection<VariantContext> getVariantContexts(final ReferenceContext ref,
                                                         final String name,
                                                         final GenomeLoc curLocation,
                                                         final boolean requireStartHere,
                                                         final boolean takeFirstOnly ) {
        return getVariantContexts(ref, Arrays.asList(name), curLocation, requireStartHere, takeFirstOnly);
    }

    public Collection<VariantContext> getVariantContexts(final ReferenceContext ref,
                                                         final Collection<String> names,
                                                         final GenomeLoc curLocation,
                                                         final boolean requireStartHere,
                                                         final boolean takeFirstOnly ) {
        Collection<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( String name : names ) {
            RODRecordList rodList = getTrackDataByName(name); // require that the name is an exact match

            if ( rodList != null )
                addVariantContexts(contexts, rodList, ref, curLocation, requireStartHere, takeFirstOnly );
        }

        return contexts;
    }

    /**
     * Gets the variant context associated with name, and assumes the system only has a single bound track at this location.  Throws an exception if not.
     * see getVariantContexts for more information.
     *
     * @param name               name
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @return variant context
     */
    public VariantContext getVariantContext(final ReferenceContext ref,
                                            final String name,
                                            final GenomeLoc curLocation,
                                            final boolean requireStartHere ) {
        Collection<VariantContext> contexts = getVariantContexts(ref, name, curLocation, requireStartHere, false );

        if ( contexts.size() > 1 )
            throw new ReviewedStingException("Requested a single VariantContext object for track " + name + " but multiple variants were present at position " + curLocation);
        else if ( contexts.size() == 0 )
            return null;
        else
            return contexts.iterator().next();
    }

    /**
     * Very simple accessor that gets the first (and only!) VC associated with name at the current location, or
     * null if there's no binding here.
     * 
     * @param ref
     * @param name
     * @param curLocation
     * @return
     */
    public VariantContext getVariantContext(final ReferenceContext ref,
                                            final String name,
                                            final GenomeLoc curLocation) {
        return getVariantContext(ref, name, curLocation, true);
    }

    private void addVariantContexts(final Collection<VariantContext> contexts,
                                    final RODRecordList rodList,
                                    final ReferenceContext ref,
                                    final GenomeLoc curLocation,
                                    final boolean requireStartHere,
                                    final boolean takeFirstOnly ) {
        for ( GATKFeature rec : rodList ) {
            if ( VariantContextAdaptors.canBeConvertedToVariantContext(rec.getUnderlyingObject()) ) {
                // ok, we might actually be able to turn this record in a variant context
                final VariantContext vc = VariantContextAdaptors.toVariantContext(rodList.getName(), rec.getUnderlyingObject(), ref);

                if ( vc == null ) // sometimes the track has odd stuff in it that can't be converted
                    continue;

                if ( ! requireStartHere || rec.getLocation().getStart() == curLocation.getStart() ) {  // ok, we are going to keep this thing
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
     * @return track data for the given rod
     */
    private RODRecordList getTrackDataByName(final String name) {
        final String luName = canonicalName(name);
        return map.get(luName);
    }

    /**
     * Returns the canonical name of the rod name (lowercases it)
     * @param name the name of the rod
     * @return canonical name of the rod
     */
    private final String canonicalName(final String name) {
        // todo -- remove me after switch to RodBinding syntax
        return name.toLowerCase();
    }
}

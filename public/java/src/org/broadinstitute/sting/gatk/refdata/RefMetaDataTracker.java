package org.broadinstitute.sting.gatk.refdata;

import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.Reference;
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
    // TODO: this should be a list, not a map, actually

    private final static RODRecordList EMPTY_ROD_RECORD_LIST = new RODRecordListImpl("EMPTY");

    final Map<String, RODRecordList> map;
    final ReferenceContext ref;
    final protected static Logger logger = Logger.getLogger(RefMetaDataTracker.class);

    // ------------------------------------------------------------------------------------------
    //
    //
    // Special ENGINE interaction functions
    //
    //
    // ------------------------------------------------------------------------------------------

    public RefMetaDataTracker(final Collection<RODRecordList> allBindings, final ReferenceContext ref) {
        this.ref = ref;
        if ( allBindings.isEmpty() )
            map = Collections.emptyMap();
        else {
            map = new HashMap<String, RODRecordList>(allBindings.size());
            for ( RODRecordList rod : allBindings ) {
                //logger.debug(String.format("Binding %s to %s", name, rod));
                if ( rod != null )
                    map.put(canonicalName(rod.getName()), maybeConvertToVariantContext(rod));
            }
        }
    }

    /**
     * A private converter that transforms a RODRecordList of objects of type X into
     * a list of VariantContexts, if possible.
     *
     * TODO: should be removed when Features like dbsnp and hapmap produce VCs directly
     *
     * @param bindings
     * @return
     */
    private final RODRecordList maybeConvertToVariantContext(RODRecordList bindings) {
        List<GATKFeature> values = new ArrayList<GATKFeature>(bindings.size());

        for ( GATKFeature rec : bindings ) {
            if ( VariantContextAdaptors.canBeConvertedToVariantContext(rec.getUnderlyingObject()) ) {
                final VariantContext vc = VariantContextAdaptors.toVariantContext(bindings.getName(), rec.getUnderlyingObject(), ref);
                if ( vc != null ) // it's possible that the conversion failed, but we continue along anyway
                    values.add(new GATKFeature.TribbleGATKFeature(ref.getGenomeLocParser(), vc, rec.getName()));
            } else
                values.add(rec);
        }

        return new RODRecordListImpl(bindings.getName(), values, bindings.getLocation());
    }

    // ------------------------------------------------------------------------------------------
    //
    //
    // Generic accessors
    //
    //
    // ------------------------------------------------------------------------------------------

    public <T extends Feature> List<T> getValues(Class<T> type) {
        return addValues(map.keySet(), type, new ArrayList<T>(), null, false, false);
    }
    public <T extends Feature> List<T> getValues(Class<T> type, final GenomeLoc onlyAtThisLoc) {
        return addValues(map.keySet(), type, new ArrayList<T>(), onlyAtThisLoc, true, false);
    }
    public <T extends Feature> List<T> getValues(Class<T> type, final String name) {
        return addValues(name, type, new ArrayList<T>(), getTrackDataByName(name), null, false, false);
    }
    public <T extends Feature> List<T> getValues(Class<T> type, final String name, final GenomeLoc onlyAtThisLoc) {
        return addValues(name, type, new ArrayList<T>(), getTrackDataByName(name), onlyAtThisLoc, true, false);
    }
    public <T extends Feature> List<T> getValues(Class<T> type, final Collection<String> names) {
        return addValues(names, type, new ArrayList<T>(), null, false, false);
    }
    public <T extends Feature> List<T> getValues(Class<T> type, final Collection<String> names, final GenomeLoc onlyAtThisLoc) {
        return addValues(names, type, new ArrayList<T>(), onlyAtThisLoc, true, false);
    }

    public <T extends Feature> T getFirstValue(Class<T> type) {
        return safeGetFirst(getValues(type));
    }
    public <T extends Feature> T getFirstValue(Class<T> type, final GenomeLoc onlyAtThisLoc) {
        return safeGetFirst(getValues(type, onlyAtThisLoc));
    }
    public <T extends Feature> T getFirstValue(Class<T> type, final String name) {
        return safeGetFirst(getValues(type, name));
    }
    public <T extends Feature> T getFirstValue(Class<T> type, final String name, final GenomeLoc onlyAtThisLoc) {
        return safeGetFirst(getValues(type, name, onlyAtThisLoc));
    }
    public <T extends Feature> T getFirstValue(Class<T> type, final Collection<String> names) {
        return safeGetFirst(getValues(type, names));
    }
    public <T extends Feature> T getFirstValue(Class<T> type, final Collection<String> names, final GenomeLoc onlyAtThisLoc) {
        return safeGetFirst(getValues(type, names, onlyAtThisLoc));
    }

    //
    // ROD binding accessors
    //
    public <T extends Feature> List<T> getValues(RodBinding<T> rodBinding) {
        return getValues(rodBinding.getType(), rodBinding.getVariableName());
    }
    public <T extends Feature> List<T> getValues(RodBinding<T> rodBinding, final GenomeLoc onlyAtThisLoc) {
        return getValues(rodBinding.getType(), rodBinding.getVariableName(), onlyAtThisLoc);
    }

    public <T extends Feature> T getFirstValue(RodBinding<T> rodBinding) {
        return getFirstValue(rodBinding.getType(), rodBinding.getVariableName());
    }
    public <T extends Feature> T getFirstValue(RodBinding<T> rodBinding, final GenomeLoc onlyAtThisLoc) {
        return getFirstValue(rodBinding.getType(), rodBinding.getVariableName(), onlyAtThisLoc);
    }

    public boolean hasValues(RodBinding rodBinding) {
        return hasValues(rodBinding.getVariableName());
    }

    public List<GATKFeature> getValuesAsGATKFeatures(RodBinding rodBinding) {
        return getValuesAsGATKFeatures(rodBinding.getVariableName());
    }

    /**
     * Helper function for getFirst() operations that takes a list of <T> and
     * returns the first element, or null if no such element exists.
     *
     * TODO: determine specific behavior for l.size() > 1.  Do we turn first or an error?
     * TODO: right now we return the first.  Should be clearer
     *
     * @param l
     * @param <T>
     * @return
     */
    final private <T extends Feature> T safeGetFirst(List<T> l) {
        // todo: should we be warning people here?  Throwing an error?
        return l.isEmpty() ? null : l.get(0);
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
    public List<GATKFeature> getAllValuesAsGATKFeatures() {
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
        return getTrackDataByName(name);
    }

    /**
     * Get all of the RMD tracks at the current site. Each track is returned as a single compound
     * object (RODRecordList) that may contain multiple RMD records associated with the current site.
     *
     * @return List of all tracks
     */
    public List<RODRecordList> getBoundRodTracks() {
        LinkedList<RODRecordList> bound = new LinkedList<RODRecordList>();

        for ( RODRecordList value : map.values() ) {
            if ( value.size() != 0 ) bound.add(value);
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
            if ( ! value.isEmpty() ) {
                n++;
            }
        }
        return n;
    }

    // ------------------------------------------------------------------------------------------
    //
    //
    // old style Generic accessors
    //
    // TODO -- DELETE ME
    //
    //
    // ------------------------------------------------------------------------------------------

    /**
     * No-assumption version of getValues(name, class).  Returns Objects.
     */
    @Deprecated
    public List<Object> getValues(final String name) {
        return (List<Object>)(List)getValues(name, Feature.class);
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
    @Deprecated
    public <T extends Feature> List<T> getValues(final String name, final Class<T> clazz) {
        RODRecordList list = getTrackDataByName(name);

        if (list.isEmpty())
            return Collections.emptyList();
        else {
            return addValues(name, clazz, new ArrayList<T>(), list, list.getLocation(), false, false);
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
     * @param <T extends Feature> the type to parameterize on, matching the clazz argument
     * @return a record of type T, or null if no record is present.
     */
    @Deprecated
    public <T extends Feature> T getFirstValue(final String name, final Class<T> clazz) {
        RODRecordList objects = getTrackDataByName(name);

        if (objects.isEmpty()) return null;

        Object obj = objects.get(0).getUnderlyingObject();
        if (!(clazz.isAssignableFrom(obj.getClass())))
            throw new UserException.CommandLineException("Unable to case track named " + name + " to type of " + clazz.toString()
                    + " it's of type " + obj.getClass());
        else
            return (T)obj;
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
     * @return variant context
     */
    @Deprecated
    public List<VariantContext> getAllVariantContexts() {
        return getAllVariantContexts(null, false, false);
    }

    /**
     * Returns all of the variant contexts that start at the current location
     *
     * @param curLocation
     * @return
     */
    @Deprecated
    public List<VariantContext> getAllVariantContexts(final GenomeLoc curLocation) {
        return getAllVariantContexts(curLocation, true, false);
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
     *
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @param takeFirstOnly      do we take the first rod only?
     * @return variant context
     */
    @Deprecated
    public List<VariantContext> getAllVariantContexts(final GenomeLoc curLocation,
                                                      final boolean requireStartHere,
                                                      final boolean takeFirstOnly) {
        List<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( RODRecordList rodList : getBoundRodTracks() ) {
            addVariantContexts(contexts, rodList, curLocation, requireStartHere, takeFirstOnly);
        }

        return contexts;
    }

    /**
     * Gets the variant contexts associated with track name name
     *
     * see getVariantContexts for more information.
     *
     *
     * @param name               name
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @param takeFirstOnly      do we take the first rod only?
     * @return variant context
     */
    @Deprecated
    public List<VariantContext> getVariantContexts(final String name,
                                                   final GenomeLoc curLocation,
                                                   final boolean requireStartHere,
                                                   final boolean takeFirstOnly) {
        return getVariantContexts(Arrays.asList(name), curLocation, requireStartHere, takeFirstOnly);
    }

    @Deprecated
    public List<VariantContext> getVariantContexts(final Collection<String> names,
                                                   final GenomeLoc curLocation,
                                                   final boolean requireStartHere,
                                                   final boolean takeFirstOnly) {
        List<VariantContext> contexts = new ArrayList<VariantContext>();

        for ( String name : names ) {
            RODRecordList rodList = getTrackDataByName(name); // require that the name is an exact match
            addVariantContexts(contexts, rodList, curLocation, requireStartHere, takeFirstOnly );
        }

        return contexts;
    }

    /**
     * Gets the variant context associated with name, and assumes the system only has a single bound track at this location.  Throws an exception if not.
     * see getVariantContexts for more information.
     *
     *
     * @param name               name
     * @param curLocation        location
     * @param requireStartHere   do we require the rod to start at this location?
     * @return variant context
     */
    @Deprecated
    public VariantContext getVariantContext(final String name,
                                            final GenomeLoc curLocation,
                                            final boolean requireStartHere) {
        List<VariantContext> contexts = getVariantContexts(name, curLocation, requireStartHere, false );

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
     *
     * @param name
     * @param curLocation
     * @return
     */
    @Deprecated
    public VariantContext getVariantContext(final String name,
                                            final GenomeLoc curLocation) {
        return getVariantContext(name, curLocation, true);
    }

    @Deprecated
    private void addVariantContexts(final List<VariantContext> contexts,
                                    final RODRecordList rodList,
                                    final GenomeLoc curLocation,
                                    final boolean requireStartHere,
                                    final boolean takeFirstOnly ) {
        addValues("xxx", VariantContext.class, contexts, rodList, curLocation, requireStartHere, takeFirstOnly);
    }

    private <T extends Feature> List<T> addValues(final Collection<String> names,
                                  final Class<T> type,
                                  final List<T> values,
                                  final GenomeLoc curLocation,
                                  final boolean requireStartHere,
                                  final boolean takeFirstOnly ) {
        for ( String name : names ) {
            RODRecordList rodList = getTrackDataByName(name); // require that the name is an exact match
            addValues(name, type, values, rodList, curLocation, requireStartHere, takeFirstOnly );
        }

        return values;
    }

    private <T extends Feature> List<T> addValues(final String name,
                                  final Class<T> type,
                                  final List<T> values,
                                  final RODRecordList rodList,
                                  final GenomeLoc curLocation,
                                  final boolean requireStartHere,
                                  final boolean takeFirstOnly ) {
        for ( GATKFeature rec : rodList ) {
            if ( ! requireStartHere || rec.getLocation().getStart() == curLocation.getStart() ) {  // ok, we are going to keep this thing
                Object obj = rec.getUnderlyingObject();
                if (!(type.isAssignableFrom(obj.getClass())))
                    throw new UserException.CommandLineException("Unable to cast track named " + name + " to type of " + type.toString()
                            + " it's of type " + obj.getClass());

                values.add((T)obj);

                if ( takeFirstOnly )
                    // we only want the first passing instance, so break the loop over records in rodList
                    break;
            }
        }

        return values;
    }

    /**
     * Finds the reference metadata track named 'name' and returns all ROD records from that track associated
     * with the current site as a RODRecordList List object. If no data track with specified name is available,
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
        RODRecordList l = map.get(luName);
        return l == null ? EMPTY_ROD_RECORD_LIST : l;
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

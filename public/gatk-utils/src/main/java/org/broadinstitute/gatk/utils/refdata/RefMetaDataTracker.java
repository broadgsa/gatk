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

package org.broadinstitute.gatk.utils.refdata;

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.*;

/**
 * This class represents the Reference Metadata available at a particular site in the genome.  It can be
 * used to conveniently lookup the RMDs at this site, as well just getting a list of all of the RMDs
 *
 * The standard interaction model is:
 *
 * Traversal system arrives at a site, which has a bunch of RMDs covering it
 * Traversal passes creates a tracker and passes it to the walker
 * walker calls get(rodBinding) to obtain the RMDs values at this site for the track
 * associated with rodBinding.
 *
 * Note that this is an immutable class.  Once created the underlying data structures
 * cannot be modified
 *
 * User: mdepristo
 * Date: Apr 3, 2009
 * Time: 3:05:23 PM
 */
public class RefMetaDataTracker {
    // TODO: this should be a list, not a bindings, actually
    private final static RODRecordList EMPTY_ROD_RECORD_LIST = new RODRecordListImpl("EMPTY");

    final Map<String, RODRecordList> bindings;
    final protected static Logger logger = Logger.getLogger(RefMetaDataTracker.class);
    public final static RefMetaDataTracker EMPTY_TRACKER = new RefMetaDataTracker();

    // ------------------------------------------------------------------------------------------
    //
    //
    // Special ENGINE interaction functions
    //
    //
    // ------------------------------------------------------------------------------------------

    /**
     * Create an tracker with no bindings
     */
    public RefMetaDataTracker() {
        bindings = Collections.emptyMap();
    }

    public RefMetaDataTracker(final Collection<RODRecordList> allBindings) {
        // set up the bindings
        if ( allBindings.isEmpty() )
            bindings = Collections.emptyMap();
        else {
            final Map<String, RODRecordList> tmap = new HashMap<String, RODRecordList>(allBindings.size());
            for ( RODRecordList rod : allBindings ) {
                if ( rod != null && ! rod.isEmpty() )
                    tmap.put(canonicalName(rod.getName()), rod);
            }

            // ensure that no one modifies the bindings itself
            bindings = Collections.unmodifiableMap(tmap);
        }
    }

    // ------------------------------------------------------------------------------------------
    //
    //
    // Generic accessors
    //
    //
    // ------------------------------------------------------------------------------------------

    /**
     * Gets all of the Tribble features spanning this locus, returning them as a list of specific
     * type T extending Feature.  This function looks across all tracks to find the Features, so
     * if you have two tracks A and B each containing 1 Feature, then getValues will return
     * a list containing both features.
     *
     * Note that this function assumes that all of the bound features are instances of or
     * subclasses of T.  A ClassCastException will occur if this isn't the case.  If you want
     * to get all Features without any danger of such an exception use the root Tribble
     * interface Feature.
     *
     * @param type The type of the underlying objects bound here
     * @param <T> as above
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"type != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getValues(final Class<T> type) {
        return addValues(bindings.keySet(), type, new ArrayList<T>(), null, false, false);
    }

    /**
     * Provides the same functionality as @link #getValues(Class<T>) but will only include
     * Features that start as the GenomeLoc provide onlyAtThisLoc.
     *
     * @param type The type of the underlying objects bound here
     * @param onlyAtThisLoc
     * @param <T> as above
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"type != null", "onlyAtThisLoc != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getValues(final Class<T> type, final GenomeLoc onlyAtThisLoc) {
        return addValues(bindings.keySet(), type, new ArrayList<T>(), onlyAtThisLoc, true, false);
    }

    /**
     * Uses the same logic as @link #getValues(Class) but arbitrary select one of the resulting
     * elements of the list to return.  That is, if there would be two elements in the result of
     * @link #getValues(Class), one of these two is selected, and which one it will be isn't
     * specified.  Consequently, this method is only really safe if (1) you absolutely know
     * that only one binding will meet the constraints of @link #getValues(Class) or (2)
     * you truly don't care which of the multiple bindings available you are going to examine.
     *
     * If there are no bindings here, getFirstValue() return null
     *
     * @param type The type of the underlying objects bound here
     * @param <T> as above
     * @return A random single element the RODs bound here, or null if none are bound.
     */
    @Requires({"type != null"})
    public <T extends Feature> T getFirstValue(final Class<T> type) {
        return safeGetFirst(getValues(type));
    }

    /**
     * Uses the same logic as @link #getValue(Class,GenomeLoc) to determine the list
     * of eligible Features and @link #getFirstValue(Class) to select a single
     * element from the interval list.
     *
     * @param type The type of the underlying objects bound here
     * @param <T> as above
     * @param onlyAtThisLoc only Features starting at this site are considered
     * @return A random single element the RODs bound here starting at onlyAtThisLoc, or null if none are bound.
     */
    @Requires({"type != null", "onlyAtThisLoc != null"})
    public <T extends Feature> T getFirstValue(final Class<T> type, final GenomeLoc onlyAtThisLoc) {
        return safeGetFirst(getValues(type, onlyAtThisLoc));
    }

    /**
     * Same logic as @link #getFirstValue(RodBinding, boolean) but prioritizes records from prioritizeThisLoc if available
     *
     * @param rodBindings Only Features coming from the tracks associated with one of rodBindings are fetched
     * @param <T> The Tribble Feature type of the rodBinding, and consequently the type of the resulting list of Features
     * @param prioritizeThisLoc only Features starting at this site are considered
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"rodBindings != null", "prioritizeThisLoc != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getPrioritizedValue(final Collection<RodBinding<T>> rodBindings, final GenomeLoc prioritizeThisLoc) {
        final List<T> results = new ArrayList<>();

        for ( final RodBinding<T> rodBinding : rodBindings ) {

            // if there's a value at the prioritized location, take it
            T value = getFirstValue(rodBinding, prioritizeThisLoc);

            // otherwise, grab any one
            if ( value == null )
                value = getFirstValue(rodBinding);

            // add if not null
            if ( value != null )
                results.add(value);
        }

        return results;
    }

    /**
     * Gets all of the Tribble features bound to RodBinding spanning this locus, returning them as
     * a list of specific type T extending Feature.
     *
     * Note that this function assumes that all of the bound features are instances of or
     * subclasses of T.  A ClassCastException will occur if this isn't the case.
     *
     * @param rodBinding Only Features coming from the track associated with this rodBinding are fetched
     * @param <T> The Tribble Feature type of the rodBinding, and consequently the type of the resulting list of Features
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"rodBinding != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getValues(final RodBinding<T> rodBinding) {
        return addValues(rodBinding.getName(), rodBinding.getType(), new ArrayList<T>(1), getTrackDataByName(rodBinding), null, false, false);
    }

    /**
     * Gets all of the Tribble features bound to any RodBinding in rodBindings,
     * spanning this locus, returning them as a list of specific type T extending Feature.
     *
     * Note that this function assumes that all of the bound features are instances of or
     * subclasses of T.  A ClassCastException will occur if this isn't the case.
     *
     * @param rodBindings Only Features coming from the tracks associated with one of rodBindings are fetched
     * @param <T> The Tribble Feature type of the rodBinding, and consequently the type of the resulting list of Features
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"rodBindings != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getValues(final Collection<RodBinding<T>> rodBindings) {
        List<T> results = new ArrayList<T>(1);
        for ( RodBinding<T> rodBinding : rodBindings )
            results.addAll(getValues(rodBinding));
        return results;
    }

    /**
     * The same logic as @link #getValues(RodBinding) but enforces that each Feature start at onlyAtThisLoc
     *
     * @param rodBinding Only Features coming from the track associated with this rodBinding are fetched
     * @param <T> The Tribble Feature type of the rodBinding, and consequently the type of the resulting list of Features
     * @param onlyAtThisLoc only Features starting at this site are considered
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"rodBinding != null", "onlyAtThisLoc != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getValues(final RodBinding<T> rodBinding, final GenomeLoc onlyAtThisLoc) {
        return addValues(rodBinding.getName(), rodBinding.getType(), new ArrayList<T>(1), getTrackDataByName(rodBinding), onlyAtThisLoc, true, false);
    }

    /**
     * The same logic as @link #getValues(List) but enforces that each Feature start at onlyAtThisLoc
     *
     * @param rodBindings Only Features coming from the tracks associated with one of rodBindings are fetched
     * @param <T> The Tribble Feature type of the rodBinding, and consequently the type of the resulting list of Features
     * @param onlyAtThisLoc only Features starting at this site are considered
     * @return A freshly allocated list of all of the bindings, or an empty list if none are bound.
     */
    @Requires({"rodBindings != null", "onlyAtThisLoc != null"})
    @Ensures("result != null")
    public <T extends Feature> List<T> getValues(final Collection<RodBinding<T>> rodBindings, final GenomeLoc onlyAtThisLoc) {
        List<T> results = new ArrayList<T>(1);
        for ( RodBinding<T> rodBinding : rodBindings )
            results.addAll(getValues(rodBinding, onlyAtThisLoc));
        return results;
    }

    /**
     * Uses the same logic as @getValues(RodBinding) to determine the list
     * of eligible Features and select a single element from the resulting set
     * of eligible features.
     *
     * @param rodBinding Only Features coming from the track associated with this rodBinding are fetched
     * @param <T> as above
     * @return A random single element the eligible Features found, or null if none are bound.
     */
    @Requires({"rodBinding != null"})
    public <T extends Feature> T getFirstValue(final RodBinding<T> rodBinding) {
        return safeGetFirst(addValues(rodBinding.getName(), rodBinding.getType(), null, getTrackDataByName(rodBinding), null, false, true));
    }

    /**
     * Uses the same logic as @getValues(RodBinding, GenomeLoc) to determine the list
     * of eligible Features and select a single element from the resulting set
     * of eligible features.
     *
     * @param rodBinding Only Features coming from the track associated with this rodBinding are fetched
     * @param <T> as above
     * @param onlyAtThisLoc only Features starting at this site are considered
     * @return A random single element the eligible Features found, or null if none are bound.
     */
    @Requires({"rodBinding != null", "onlyAtThisLoc != null"})
    public <T extends Feature> T getFirstValue(final RodBinding<T> rodBinding, final GenomeLoc onlyAtThisLoc) {
        return safeGetFirst(addValues(rodBinding.getName(), rodBinding.getType(), null, getTrackDataByName(rodBinding), onlyAtThisLoc, true, true));
    }

    /**
     * Uses the same logic as @getValues(List) to determine the list
     * of eligible Features and select a single element from the resulting set
     * of eligible features.
     *
     * @param rodBindings Only Features coming from the tracks associated with these rodBindings are fetched
     * @param <T> as above
     * @return A random single element the eligible Features found, or null if none are bound.
     */
    @Requires({"rodBindings != null"})
    public <T extends Feature> T getFirstValue(final Collection<RodBinding<T>> rodBindings) {
        for ( RodBinding<T> rodBinding : rodBindings ) {
            T val = getFirstValue(rodBinding);
            if ( val != null )
                return val;
        }
        return null;
    }

    /**
     * Uses the same logic as @getValues(RodBinding,GenomeLoc) to determine the list
     * of eligible Features and select a single element from the resulting set
     * of eligible features.
     *
     * @param rodBindings Only Features coming from the tracks associated with these rodBindings are fetched
     * @param <T> as above
     * @param onlyAtThisLoc only Features starting at this site are considered
     * @return A random single element the eligible Features found, or null if none are bound.
     */
    @Requires({"rodBindings != null", "onlyAtThisLoc != null"})
    public <T extends Feature> T getFirstValue(final Collection<RodBinding<T>> rodBindings, final GenomeLoc onlyAtThisLoc) {
        for ( RodBinding<T> rodBinding : rodBindings ) {
            T val = getFirstValue(rodBinding, onlyAtThisLoc);
            if ( val != null )
                return val;
        }
        return null;
    }

    /**
     * Is there a binding at this site to a ROD/track with the specified name?
     *
     * @param rodBinding the rod binding we want to know about
     * @return true if any Features are bound in this tracker to rodBinding
     */
    @Requires({"rodBinding != null"})
    public boolean hasValues(final RodBinding rodBinding) {
        return bindings.containsKey(canonicalName(rodBinding.getName()));
    }

    /**
     * Get all of the RMD tracks at the current site. Each track is returned as a single compound
     * object (RODRecordList) that may contain multiple RMD records associated with the current site.
     *
     * @return List of all tracks
     */
    public List<RODRecordList> getBoundRodTracks() {
        return new ArrayList<RODRecordList>(bindings.values());
    }

    /**
     * The number of tracks with at least one value bound here
     * @return the number of tracks with at least one bound Feature
     */
    public int getNTracksWithBoundFeatures() {
        return bindings.size();
    }

    // ------------------------------------------------------------------------------------------
    // Protected accessors using strings for unit testing
    // ------------------------------------------------------------------------------------------

    protected boolean hasValues(final String name) {
        return bindings.containsKey(canonicalName(name));
    }

    protected <T extends Feature> List<T> getValues(final Class<T> type, final String name) {
        return addValues(name, type, new ArrayList<T>(), getTrackDataByName(name), null, false, false);
    }

    protected <T extends Feature> List<T> getValues(final Class<T> type, final String name, final GenomeLoc onlyAtThisLoc) {
        return addValues(name, type, new ArrayList<T>(), getTrackDataByName(name), onlyAtThisLoc, true, false);
    }

    protected <T extends Feature> T getFirstValue(final Class<T> type, final String name) {
        return safeGetFirst(getValues(type, name));
    }

    protected <T extends Feature> T getFirstValue(final Class<T> type, final String name, final GenomeLoc onlyAtThisLoc) {
        return safeGetFirst(getValues(type, name, onlyAtThisLoc));
    }

    // ------------------------------------------------------------------------------------------
    //
    //
    // Private utility functions
    //
    //
    // ------------------------------------------------------------------------------------------

    /**
     * Helper function for getFirst() operations that takes a list of <T> and
     * returns the first element, or null if no such element exists.
     *
     * @param l
     * @param <T>
     * @return
     */
    @Requires({"l != null"})
    private <T extends Feature> T safeGetFirst(final List<T> l) {
        return l.isEmpty() ? null : l.get(0);
    }

    private <T extends Feature> List<T> addValues(final Collection<String> names,
                                                  final Class<T> type,
                                                  List<T> values,
                                                  final GenomeLoc curLocation,
                                                  final boolean requireStartHere,
                                                  final boolean takeFirstOnly ) {
        for ( String name : names ) {
            RODRecordList rodList = getTrackDataByName(name); // require that the name is an exact match
            values = addValues(name, type, values, rodList, curLocation, requireStartHere, takeFirstOnly );
            if ( takeFirstOnly && ! values.isEmpty() )
                break;
        }

        return values;
    }



    private <T extends Feature> List<T> addValues(final String name,
                                                  final Class<T> type,
                                                  List<T> values,
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

                T objT = (T)obj;
                if ( takeFirstOnly ) {
                    if ( values == null )
                        values = Arrays.asList(objT);
                    else
                        values.add(objT);

                    break;
                } else {
                    if ( values == null )
                        values = new ArrayList<T>();
                    values.add(objT);
                }
            }
        }

        return values == null ? Collections.<T>emptyList() : values;
    }

    /**
     * Finds the reference metadata track named 'name' and returns all ROD records from that track associated
     * with the current site as a RODRecordList List object. If no data track with specified name is available,
     * returns defaultValue wrapped as RODRecordList object. NOTE: if defaultValue is null, it will be wrapped up
     * with track name set to 'name' and location set to null; otherwise the wrapper object will have name and
     * location set to defaultValue.getID() and defaultValue.getLocation(), respectively (use caution,
     * defaultValue.getLocation() may be not equal to what RODRecordList's location would be expected to be otherwise:
     * for instance, on locus traversal, location is usually expected to be a single base we are currently looking at,
     * regardless of the presence of "extended" RODs overlapping with that location).
     * @param name                track name
     * @return track data for the given rod
     */
    private RODRecordList getTrackDataByName(final String name) {
        final String luName = canonicalName(name);
        RODRecordList l = bindings.get(luName);
        return l == null ? EMPTY_ROD_RECORD_LIST : l;
    }

    private RODRecordList getTrackDataByName(final RodBinding binding) {
        return getTrackDataByName(binding.getName());
    }

    /**
     * Returns the canonical name of the rod name (lowercases it)
     * @param name the name of the rod
     * @return canonical name of the rod
     */
    private String canonicalName(final String name) {
        // todo -- remove me after switch to RodBinding syntax
        return name.toLowerCase();
    }
}

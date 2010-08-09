/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.refdata.tracks;

import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.util.*;

/**
 * Find the available track builders, and create the requisite tracks from the command line.
 *
 * In Tribble RMD tracks have two classes:
 *  - a Feature that is the model/view for the data
 *  - a Codec that is the controller to generate the Feature.
 *
 * In this class, the track types are the Codecs.  The track record types are the Features.
 */
public class RMDTrackManager extends PluginManager<RMDTrackBuilder> {
    // the input strings we use to create RODs from
    List<RMDTriplet> inputs = new ArrayList<RMDTriplet>();

    // create an active mapping of builder instances, and a map of the name -> class for convenience
    /** the tracks that are available to us, associated with their builder */
    Map<String, RMDTrackBuilder> availableTrackBuilders;
    /** the classes names, with their class description (think the Controller Codecs) */
    Map<String, Class> availableTrackTypes;
    /** the available track record types (think the Model/View Features) */
    Map<String, Class> availableTrackRecordTypes;

    /** Create a new track plugin manager. */
    public RMDTrackManager() {
        super(RMDTrackBuilder.class, "TrackBuilders", null);
    }

    /**
     * find the associated reference meta data
     *
     * @param triplets the triplets of strings from the -B command line option
     *
     * @return a list of RMDTracks, one for each -B option
     */
    public List<RMDTrack> getReferenceMetaDataSources(List<String> triplets) {
        initializeTrackTypes();
        initializeTriplets(triplets);
        // try and make the tracks given their requests
        return createRequestedTrackObjects();
    }


    /**
     * Returns a collection of track names that match the record type.
     * @param trackRecordType the record type specified in the @RMD annotation
     * @return a collection of available track record type names that match the record type
     */
    public Collection<String> getTrackRecordTypeNames(Class trackRecordType) {
        initializeTrackTypes();
        initializeTrackRecordTypes();
        Set<String> names = new TreeSet<String>();
        for (Map.Entry<String, Class> availableTrackRecordType: availableTrackRecordTypes.entrySet()) {
            if (trackRecordType.isAssignableFrom(availableTrackRecordType.getValue()))
                names.add(availableTrackRecordType.getKey());
        }
        return names;
    }

    /**
     * initialize our lists of triplets
     * @param triplets the input to the GATK, as a list of strings passed in through the -B options
     */
    private void initializeTriplets(List<String> triplets) {
        // NOTE: Method acts as a static.  Once the inputs have been passed once they are locked in.
        if (inputs.size() > 0 || triplets.size() == 0)
            return;

        for (String value: triplets) {
            String[] split = value.split(",");
            if (split.length != 3) throw new IllegalArgumentException(value + " is not a valid reference metadata track description");
            inputs.add(new RMDTriplet(split[0], split[1], split[2]));
        }
    }

    /**
     * initialize our lists of tracks and builders
     */
    private void initializeTrackTypes() {
        if (availableTrackBuilders != null && availableTrackTypes != null)
            return;

        // create an active mapping of builder instances, and a map of the name -> class for convenience
        availableTrackBuilders = new HashMap<String, RMDTrackBuilder>();
        availableTrackTypes = new HashMap<String, Class>();
        createBuilderObjects();
    }

    /**
     * create the builder objects from the retrieved list
     */
    private void createBuilderObjects() {
        // create a track builder instance for each track builder, and find out what tracks we can make
        for (String builderName : this.pluginsByName.keySet()) {
            RMDTrackBuilder builder = this.createByName(builderName);
            Map<String, Class> mapping = builder.getAvailableTrackNamesAndTypes();
            for (String name : mapping.keySet()) {
                availableTrackBuilders.put(name.toUpperCase(), builder);
                availableTrackTypes.put(name.toUpperCase(), mapping.get(name));
            }
        }
    }

    /**
     * initialize our list of track record types
     */
    private void initializeTrackRecordTypes() {
        if (availableTrackRecordTypes != null)
            return;

        availableTrackRecordTypes = new HashMap<String, Class>();
        for (RMDTrackBuilder builder : availableTrackBuilders.values()) {
            Map<String, Class> mapping = builder.getAvailableTrackNamesAndRecordTypes();
            for (String name : mapping.keySet()) {
                availableTrackRecordTypes.put(name.toUpperCase(), mapping.get(name));
            }
        }
    }

    /**
     * create the requested track objects
     *
     * @return a list of the tracks, one for each of the requested input tracks
     */
    private List<RMDTrack> createRequestedTrackObjects() {
        // create of live instances of the tracks
        List<RMDTrack> tracks = new ArrayList<RMDTrack>();

        // create instances of each of the requested types
        for (RMDTriplet trip : inputs) {
            RMDTrackBuilder b = availableTrackBuilders.get(trip.getType().toUpperCase());
            if (b == null) throw new StingException("Unable to find track for " + trip.getType());
            tracks.add(b.createInstanceOfTrack(availableTrackTypes.get(trip.getType().toUpperCase()), trip.getName(), new File(trip.getFile())));
        }
        return tracks;
    }
}

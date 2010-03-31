/*
 * Copyright (c) 2010.  The Broad Institute
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
 * THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.refdata.tracks;

import org.broadinstitute.sting.gatk.refdata.tracks.builders.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.refdata.utils.RMDTriplet;
import org.broadinstitute.sting.utils.PluginManager;
import org.broadinstitute.sting.utils.StingException;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class RMDTrackManager
 *         <p/>
 *         Find the available track builders, and create the requisite tracks from the command line.
 */
public class RMDTrackManager extends PluginManager<RMDTrackBuilder> {
    // the input strings we use to create RODs from
    List<RMDTriplet> inputs = new ArrayList<RMDTriplet>();

    // create an active mapping of builder instances, and a map of the name -> class for convenience
    Map<String, RMDTrackBuilder> availableTracks;
    Map<String, Class> availableTrackClasses;

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
        if (availableTracks == null || availableTrackClasses == null) initialize(triplets);
        // try and make the tracks given their requests
        return createRequestedTrackObjects(availableTracks, availableTrackClasses);
    }

    /**
     * initialize our lists of tracks and builders
     * @param triplets the input to the GATK, as a list of strings passed in through the -B options
     */
    private void initialize(List<String> triplets) {        
        for (String value: triplets) {
            String[] split = value.split(",");
            if (split.length != 3) throw new IllegalArgumentException(value + " is not a valid reference metadata track description");
            inputs.add(new RMDTriplet(split[0], split[1], split[2]));
        }

        // create an active mapping of builder instances, and a map of the name -> class for convenience
        availableTracks = new HashMap<String, RMDTrackBuilder>();
        availableTrackClasses = new HashMap<String, Class>();
        createBuilderObjects();


    }

    /**
     * create the builder objects from the retrieved list
     */
    private void createBuilderObjects() {
        // create a track builder instance for each track builder, and find out what tracks we can make
        for (String builderName : this.pluginsByName.keySet()) {
            RMDTrackBuilder builder = this.createByName(builderName);
            for (String name : builder.getAvailableTrackNamesAndTypes().keySet()) {
                availableTracks.put(name.toUpperCase(), builder);
                availableTrackClasses.put(name.toUpperCase(), builder.getAvailableTrackNamesAndTypes().get(name));
            }
        }
    }

    /**
     * create the requested track objects
     *
     * @param availableTracks       the tracks that are available to us, associated with their builder
     * @param availableTrackClasses the classes names, with their class description
     *
     * @return a list of the tracks, one for each of the requested input tracks
     */
    private List<RMDTrack> createRequestedTrackObjects(Map<String, RMDTrackBuilder> availableTracks, Map<String, Class> availableTrackClasses) {
        // create of live instances of the tracks
        List<RMDTrack> tracks = new ArrayList<RMDTrack>();

        // create instances of each of the requested types
        for (RMDTriplet trip : inputs) {
            RMDTrackBuilder b = availableTracks.get(trip.getType().toUpperCase());
            if (b == null) throw new StingException("Unable to find track for " + trip.getType());
            tracks.add(b.createInstanceOfTrack(availableTrackClasses.get(trip.getType().toUpperCase()), trip.getName(), new File(trip.getFile())));
        }
        return tracks;
    }
}


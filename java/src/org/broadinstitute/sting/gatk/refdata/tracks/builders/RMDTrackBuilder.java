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

package org.broadinstitute.sting.gatk.refdata.tracks.builders;

import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException;

import java.io.File;
import java.util.Map;



/**
 * @author aaron
 *         <p/>
 *         Interface RMDTrackBuilder
 *         <p/>
 *         The basic interface for finding and parsing RMDTracks.  Track builders present an interface that allows
 *         the track manager to find and create tracks of the specified type.
 */
public interface RMDTrackBuilder {

    /** @return a list of all available tracks types we currently have access to create */
    public Map<String, Class> getAvailableTrackNamesAndTypes();

    /**
     * create a RMDTrack of the specified type
     *
     * @param targetClass the target class of track
     * @param name        what to call the track
     * @param inputFile   the input file
     *
     * @return an instance of the track
     * @throws RMDTrackCreationException if we don't know of the target class or we couldn't create it
     */
    public RMDTrack createInstanceOfTrack(Class targetClass, String name, File inputFile) throws RMDTrackCreationException;
}

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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException;
import org.broadinstitute.sting.gatk.refdata.tracks.RODRMDTrack;
import org.broadinstitute.sting.oneoffprojects.refdata.HapmapVCFROD;

import java.io.File;
import java.util.HashMap;
import java.util.Map;


/**
 * @author aaron
 *         <p/>
 *         Class RODTrackBuilder
 *         <p/>
 *        the builder for tracks of the current ROD system, a holdover until Tribble supports binary and multi-line formats
 */
public class RODTrackBuilder implements RMDTrackBuilder {

    /** our log, which we want to capture anything from this class */
    private static Logger logger = Logger.getLogger(ReferenceOrderedData.class);

    /**
     * the bindings from track name to the ROD class we use
     */
    private static HashMap<String, Class<? extends ReferenceOrderedDatum>> Types = new HashMap<String, Class<? extends ReferenceOrderedDatum>>();

    static {
        // All known ROD types
        Types.put("GELI", rodGELI.class);
        Types.put("RefSeq", rodRefSeq.class);
        Types.put("Table", TabularROD.class);
        Types.put("HapMap", HapMapROD.class);
        Types.put("Intervals", IntervalRod.class);
        Types.put("GLF", RodGLF.class);
        Types.put("PicardDbSNP", rodPicardDbSNP.class);
        Types.put("HapmapVCF", HapmapVCFROD.class);
        Types.put("Beagle", BeagleROD.class);
        Types.put("Plink", PlinkRod.class);
        Types.put("Bed", RodBed.class);
    }

    /**
        * create a RMDTrack of the specified type
        *
        * @param targetClass the target class of track
        * @param name        what to call the track
        * @param inputFile   the input file
        *
        * @return an instance of the track
        * @throws org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackCreationException
        *          if we don't know of the target class or we couldn't create it
        */
       //@Override
       public RMDTrack createInstanceOfTrack(Class targetClass, String name, File inputFile) throws RMDTrackCreationException {
           return new RODRMDTrack(targetClass, name, inputFile, createROD(name,targetClass,inputFile));
       }

    /** @return a map of all available tracks we currently have access to create */
    public Map<String, Class> getAvailableTrackNamesAndTypes() {
        Map<String, Class> ret = new HashMap<String, Class>();
        for (String name : Types.keySet())
            ret.put(name, Types.get(name));
        return ret;
    }

/**
     * Helpful function that parses a single triplet of <name> <type> <file> and returns the corresponding ROD with
     * <name>, of type <type> that reads its input from <file>.
     *
     * @param trackName the name of the track to create
     * @param type the type of the track to create
     * @param fileName the filename to create the track from
     * @return a reference ordered data track
     */
    public ReferenceOrderedData createROD(final String trackName, Class type, File fileName) {

        // Create the ROD
        ReferenceOrderedData<?> rod = new ReferenceOrderedData<ReferenceOrderedDatum>(trackName, fileName, type );
        logger.info(String.format("Created binding from %s to %s of type %s", trackName, fileName, type));
        return rod;
    }

}

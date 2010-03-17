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

    public static HashMap<String, ReferenceOrderedData.RODBinding> Types = new HashMap<String, ReferenceOrderedData.RODBinding>();

    public static void addModule(final String name, final Class<? extends ReferenceOrderedDatum> rodType) {
        final String boundName = name.toLowerCase();
        if (Types.containsKey(boundName)) {
            throw new RuntimeException(String.format("GATK BUG: adding ROD module %s that is already bound", boundName));
        }
        logger.info(String.format("* Adding rod class %s", name));
        Types.put(boundName, new ReferenceOrderedData.RODBinding(name, rodType));
    }

    static {
        // All known ROD types
        addModule("GFF", RodGenotypeChipAsGFF.class);
        //addModule("dbSNP", rodDbSNP.class);
        addModule("HapMapAlleleFrequencies", HapMapAlleleFrequenciesROD.class);
        addModule("SAMPileup", rodSAMPileup.class);
        addModule("GELI", rodGELI.class);
        addModule("RefSeq", rodRefSeq.class);
        addModule("Table", TabularROD.class);
        addModule("PooledEM", PooledEMSNPROD.class);
        addModule("CleanedOutSNP", CleanedOutSNPROD.class);
        addModule("Sequenom", SequenomROD.class);
        addModule("SangerSNP", SangerSNPROD.class);
        addModule("SimpleIndel", SimpleIndelROD.class);
        addModule("PointIndel", PointIndelROD.class);
        addModule("HapMapGenotype", HapMapGenotypeROD.class);
        addModule("Intervals", IntervalRod.class);
        addModule("Variants", RodGeliText.class);
        addModule("GLF", RodGLF.class);
        addModule("VCF", RodVCF.class);
        addModule("PicardDbSNP", rodPicardDbSNP.class);
        addModule("HapmapVCF", HapmapVCFROD.class);
        addModule("Beagle", BeagleROD.class);
        addModule("Plink", PlinkRod.class);
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
           return new RODRMDTrack(targetClass, name, inputFile, ReferenceOrderedData.parse1Binding(name,targetClass.getName(),inputFile.getAbsolutePath()));
       }

/** @return a map of all available tracks we currently have access to create */
    //@Override
    public Map<String, Class> getAvailableTrackNamesAndTypes() {
        Map<String, Class> ret = new HashMap<String, Class>();
        for (ReferenceOrderedData.RODBinding binding: Types.values())
            ret.put(binding.name, binding.type);
        return ret;
    }    
}

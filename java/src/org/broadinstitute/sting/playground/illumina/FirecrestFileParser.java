/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package org.broadinstitute.sting.playground.illumina;

import edu.mit.broad.picard.util.PasteParser;
import edu.mit.broad.picard.util.FormatUtil;
import edu.mit.broad.picard.util.BasicTextFileParser;
import edu.mit.broad.picard.PicardException;

import java.io.File;

/**
 * Class to parse the data in an Illumina Firecrest directory and return an iterator over that data, in order
 * by tile.
 *
 * @author Kiran Garimella
 */
public class FirecrestFileParser extends AbstractFirecrestFileParser {

    private BasicTextFileParser parser;
    private final FormatUtil formatter = new FormatUtil();
    private final File[] intensityFiles;
    private int cycle_start = -1;
    private int cycle_stop = -1;

    /**
     * Constructor
     *
     * @param firecrestDirectory  directory where the Firecrest files can be located
     * @param lane              the lane to parse
     */
    public FirecrestFileParser(final File firecrestDirectory, final int lane) {
        super(firecrestDirectory, lane);
        intensityFiles = getFilesMatchingRegexp("s_" + lane + "_\\d{4}_int.txt(.gz)?");
    }

    public FirecrestFileParser(final File firecrestDirectory, final int lane, final int cycle_start, final int cycle_stop) {
        super(firecrestDirectory, lane);
        intensityFiles = getFilesMatchingRegexp("s_" + lane + "_\\d{4}_int.txt(.gz)?");
        this.cycle_start = cycle_start;
        this.cycle_stop = cycle_stop;
    }

    @Override
    public boolean isValidFirecrestDirectory() {
        return (intensityFiles.length > 0);
    }

    /**
     * Sorts the relevant files in the firecrestDirectory.  Does some basic sanity checking to ensure that some files
     * are found and that they are the expected multiple for paired-end or not.
     *
     */
    @Override
    protected void prepareToIterate() {
        // Some basic sanity checking on file counts
        if (intensityFiles.length == 0) {
            throw new PicardException("No Firecrest 1.3 intensity files found in " + firecrestDirectory.getAbsolutePath() + " for lane " + lane);
        }

        // Sort each set of reads and create a text parser for it
        parser = makeParserForTextFiles(true, intensityFiles);
    }

    /**
     * Parses the next line from the parser and constructs a FirecrestReadData object from it
     * The first 4 fields are position information for the read, and the remaining value are
     * the intensities data.
     *
     * @return  a fully populated FirecrestReadData object
     */
    protected FirecrestReadData readNext() {
        if (!parser.hasNext()) {
                return null;
        }
        final String[] data = parser.next();
        final int lane = formatter.parseInt(data[0]);
        final int tile = formatter.parseInt(data[1]);
        final int x = formatter.parseInt(data[2]);
        final int y = formatter.parseInt(data[3]);

        int intensityOffset = 4;
        int numIntensities = (data.length - 4)/4;
        FourIntensity[] intensities = new FourIntensity[numIntensities];

        for (int cycle = 0, index = intensityOffset; cycle < numIntensities; cycle++) {
            float[] fIntensities = new float[4];
            for (int channel = 0; channel < 4; channel++, index++) {
                fIntensities[channel] = formatter.parseFloat(data[index]);
            }

            intensities[cycle] = new FourIntensity(fIntensities);
        }

        FourIntensity[] intensities2 = new FourIntensity[(cycle_start > 0 && cycle_stop > 0 && cycle_stop > cycle_start) ? numIntensities : (cycle_stop - cycle_start)];
        for (int cycle = 0, offset = (cycle_start >= 0 ? cycle_start : 0); cycle < intensities2.length; cycle++) {
            intensities2[cycle] = intensities[offset + cycle];
        }

        return new FirecrestReadData(lane, tile, x, y, intensities2);
    }

    /**
     * Closes the underlying PasteParser
     */
    @Override
    public void close() {
        if (parser != null) {
            parser.close();
        }
    }

}

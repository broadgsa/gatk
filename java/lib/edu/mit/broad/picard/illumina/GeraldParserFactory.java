/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.picard.illumina;

import edu.mit.broad.picard.io.IoUtil;

import java.io.File;

/**
 * Given a Gerald directory, create a GeraldParser for one end or both ends as appropriate.
 */
public class GeraldParserFactory {

    // A Map of squashed reference chunk to reference genome sequence/chromosome. The chunk is represented as
    // a mapping (sequence=chunk file, startPos=offset into chunk file).
    private final SquashedCoordinateMap geraldToArachne;
    private final File geraldDir;
    private final int lane;

    public GeraldParserFactory(final File geraldDir, final int lane, final File squashedMapFile) {
        this.geraldDir = geraldDir;
        this.lane = lane;
        geraldToArachne = new SquashedCoordinateMap(squashedMapFile);
    }

    /** Attempts to determine if an analysis on a lane is PE or single. */
    public boolean isPairedRun() {
        if (new File(geraldDir, "s_" + lane + "_1_eland_query.txt").exists()) return true;
        else if (new File(geraldDir, "s_" + lane + "_eland_query.txt").exists()) return false;

        throw new IllegalStateException("Could not determine if gerald run is PE or fragment.");
    }

    private String makeLanePrefix(final Integer readNumber) {
        return "s_" + lane + "_" + (readNumber == null ? "" : readNumber + "_");

    }

    /**
     * @param readNumber 1 == first end of pair; 2 == second end of pair; null == unpaired
     * @return a GeraldParser for the given end
     */
    public GeraldParser makeParser(final Integer readNumber) {
        final File elandExtendedFile = new File(geraldDir, makeLanePrefix(readNumber) + "eland_extended.txt");
        final File exportFile = new File(geraldDir, makeLanePrefix(readNumber) + "export.txt");
        IoUtil.assertFileIsReadable(elandExtendedFile);
        IoUtil.assertFileIsReadable(exportFile);
        return new GeraldParser(geraldToArachne, elandExtendedFile, exportFile);
    }

}

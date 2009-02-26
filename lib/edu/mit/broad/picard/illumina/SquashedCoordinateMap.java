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

import edu.mit.broad.sam.util.CoordMath;
import edu.mit.broad.picard.cmdline.CommandLineUtils;

import java.util.Map;
import java.util.HashMap;
import java.io.File;
import java.io.BufferedReader;
import java.io.IOException;

public class SquashedCoordinateMap {
    private final Map<SimpleMapping, String> geraldToArachne = new HashMap<SimpleMapping, String>();
    private long genomeSize;

    public SquashedCoordinateMap(final File squashedMapFile) {
        try {
            final BufferedReader in = CommandLineUtils.getReader(squashedMapFile);
            String line;
            genomeSize = 0;

            while ((line = in.readLine()) != null) {
                final String[] fields = CommandLineUtils.SPACE_SPLITTER.split(line);
                final String arachneIndex = fields[0].trim().intern();
                final String squashedRefIndex = fields[1].trim().intern();
                final long squashedStart = Long.parseLong(fields[2]);
                final long length = Long.parseLong(fields[3]);
                final String sequenceName = fields[4];

                final SimpleMapping mapping = new SimpleMapping(squashedRefIndex, squashedStart, 
                        CoordMath.getEnd(squashedStart, length), sequenceName);
                geraldToArachne.put(mapping, arachneIndex);

                genomeSize += length;
            }

            in.close();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /* Converts a read's mapping from Gerald's vretarded space to arachne index + coords. */
    public void convertToArachneCoords(final SimpleMapping read) {
        if (this.geraldToArachne == null || this.geraldToArachne.isEmpty()) {
            throw new IllegalStateException("Cannot invoke convertToArachneCoords before parseSquashedMapFile");
        }

        for (final Map.Entry<SimpleMapping,String> entry : this.geraldToArachne.entrySet()) {
            final SimpleMapping chunk = entry.getKey();
            if (chunk.intersects(read)) {
                read.setArachneIndex(entry.getValue());
                read.setStartPos( read.getStartPos() - chunk.getStartPos() );
                read.setEndPos(   read.getEndPos()   - chunk.getStartPos() );
                read.setSequenceName(chunk.getSequenceName());
                return;
            }
        }

        throw new RuntimeException("Could not convert read: " + read);
    }

    long getGenomeSize() {
        return genomeSize;
    }
}

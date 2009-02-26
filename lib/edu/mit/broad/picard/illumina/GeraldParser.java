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

import edu.mit.broad.picard.util.PasteParser;
import edu.mit.broad.picard.util.TabbedTextFileParser;
import edu.mit.broad.picard.PicardException;
import edu.mit.broad.sam.util.CloseableIterator;

import java.io.File;
import java.util.Iterator;
import java.util.Arrays;
import java.util.regex.Pattern;
import java.text.ParsePosition;
import java.text.NumberFormat;

/**
 * Parse the pair of files (eland_extended.txt and export.txt) that correspond to an end of a Gerald run for a lane.
 */
public class GeraldParser implements Iterable<GeraldParser.GeraldAlignment>, CloseableIterator<GeraldParser.GeraldAlignment> {
    private static final int EXPECTED_ELAND_FIELDS = 4;
    // Regex used to split apart multiple alignments in the eland output
    private static final Pattern ALIGN_SPLITTER = Pattern.compile("\\,+");

    // export.txt constants
    private static final int PASSING_FILTER_COLUMN = 21;
    private static final int QUALITIES_COLUMN = 9;
    private static final int REQUIRED_EXPORT_COLUMNS = PASSING_FILTER_COLUMN + 1;

    private final NumberFormat integerFormat = NumberFormat.getIntegerInstance();

    private final SquashedCoordinateMap geraldToArachne;
    private final PasteParser pasteParser;
    private final File elandExtended;
    private final File export;
    private boolean iteratorCalled = false;
    private final byte[] solexaToPhredQualityConverter = new SolexaQualityConverter().getSolexaToPhredConversionTable();

    /**
     * @param geraldToArachne for converting btw Gerald coordinate and genomic coordinate
     */
    public GeraldParser(final SquashedCoordinateMap geraldToArachne, final File elandExtended, final File export) {
        this.geraldToArachne = geraldToArachne;
        this.elandExtended = elandExtended;
        this.export = export;
        final TabbedTextFileParser[] parsers = {
                new TabbedTextFileParser(false, elandExtended),
                new TabbedTextFileParser(false, export)
        };
        pasteParser = new PasteParser(parsers);
    }

    public Iterator<GeraldAlignment> iterator() {
        if (iteratorCalled) {
            throw new IllegalStateException("iterator() cannot be called more than once on a GeraldParser instance.");
        }
        iteratorCalled = true;
        return this;
    }

    public void close() {
        pasteParser.close();
    }

    public boolean hasNext() {
        return pasteParser.hasNext();
    }

    public GeraldAlignment next() {
        final GeraldAlignment ret = new GeraldAlignment();
        final String[][] fields = pasteParser.next();

        // Parse eland_extended.txt fields
        final String[] elandExtendedFields = fields[0];
        if (elandExtendedFields.length < EXPECTED_ELAND_FIELDS) {
            throw new PicardException("Not enough fields in file: " + elandExtended);
        }

        ret.readName   = elandExtendedFields[0].substring(1);
        ret.readBases = elandExtendedFields[1];
        ret.readLength = ret.readBases.length();
        final String[] alignCounts = elandExtendedFields[2].split(":");
        if (alignCounts.length == 3) {
            ret.zeroMismatchPlacements = Short.parseShort(alignCounts[0]);
            ret.oneMismatchPlacements  = Short.parseShort(alignCounts[1]);
            ret.twoMismatchPlacements  = Short.parseShort(alignCounts[2]);
        }

        final String[] alignments = ALIGN_SPLITTER.split(elandExtendedFields[3]);
        if (alignments.length == 1 && !"-".equals(alignments[0])) {
            final int lastDot   = alignments[0].lastIndexOf(".");
            final int colon     = alignments[0].indexOf(':');

            final String tmp = alignments[0].substring(colon + 1);
            final ParsePosition pos = new ParsePosition(0);
            final long start = integerFormat.parse(tmp, pos).longValue();
            if (pos.getIndex() == 0) {
                throw new RuntimeException("Problem parsing eland extended alignment record: " + Arrays.toString(elandExtendedFields));
            }

            final SimpleMapping m = new SimpleMapping(alignments[0].substring(lastDot+1, colon).trim(),
                    start, start + ret.readLength - 1, null);
            geraldToArachne.convertToArachneCoords(m);
            ret.primaryChrom = m.getSequenceName();
            ret.primaryStart = m.getStartPos();
            ret.primaryStop  = m.getEndPos();
            ret.orientation  = tmp.substring(pos.getIndex(), pos.getIndex() + 1);
            ret.mismatchString = tmp.substring(pos.getIndex() + 1);

            // Count the mismatches in the alignment
            for (int i=pos.getIndex(); i<tmp.length(); ++i) {
                final char ch = tmp.charAt(i);
                if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T') {
                    ret.primaryMismatches += 1;
                }
            }
        }

        final String[] exportFields = fields[1];
        // Parse export.txt fields
        if (exportFields.length < REQUIRED_EXPORT_COLUMNS) {
            throw new RuntimeException("Not enough columns in _export.txt file " + export);
        }
        if (exportFields[PASSING_FILTER_COLUMN].equals("Y")) {
            ret.passingFilter = true;
        } else if (exportFields[PASSING_FILTER_COLUMN].equals("N")) {
            ret.passingFilter = false;
        } else {
            throw new RuntimeException("Strange value for PF column in _export.txt file " +  export + ": '" +
                    exportFields[PASSING_FILTER_COLUMN] + "'.");
        }
        ret.phredQualities = exportFields[QUALITIES_COLUMN].getBytes();
        decodeSolexaQualitiesToPhred(ret.phredQualities);



        return ret;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    /** Decodes an array of solexa quality chars into SOLEXA numeric space.
     * Decode in place in order to avoid extra object allocation */
    private void decodeSolexaQualitiesToPhred(final byte[] solexaQuals) {
        for (int i=0; i<solexaQuals.length; ++i) {
            solexaQuals[i] = solexaToPhredQualityConverter[solexaQuals[i]];
        }

    }

    public class GeraldAlignment {
        // From eland_extended.txt
        private String readName = null;
        private String readBases = null;
        private int readLength = 0;
        private short zeroMismatchPlacements = 0;
        private short oneMismatchPlacements = 0;
        private short twoMismatchPlacements = 0;
        private String primaryChrom = null;
        private long primaryStart = 0;
        private long primaryStop = 0;
        private String orientation = null;
        private short primaryMismatches = 0;
        private String mismatchString = null;

        // from export.txt
        private boolean passingFilter;
        private byte[] phredQualities;

        public String getMismatchString() {
            return mismatchString;
        }

        public short getOneMismatchPlacements() {
            return oneMismatchPlacements;
        }

        public String getOrientation() {
            return orientation;
        }

        public boolean isPassingFilter() {
            return passingFilter;
        }

        public byte[] getPhredQualities() {
            return phredQualities;
        }

        public String getPrimaryChrom() {
            return primaryChrom;
        }

        public short getPrimaryMismatches() {
            return primaryMismatches;
        }

        public long getPrimaryStart() {
            return primaryStart;
        }

        public long getPrimaryStop() {
            return primaryStop;
        }

        public String getReadBases() {
            return readBases;
        }

        public int getReadLength() {
            return readLength;
        }

        public String getReadName() {
            return readName;
        }

        public short getTwoMismatchPlacements() {
            return twoMismatchPlacements;
        }

        public short getZeroMismatchPlacements() {
            return zeroMismatchPlacements;
        }
    }
}

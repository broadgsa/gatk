/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/
package edu.mit.broad.sam;

import edu.mit.broad.sam.util.LineReader;
import edu.mit.broad.sam.util.RuntimeIOException;
import edu.mit.broad.sam.util.StringUtil;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * This is actually two classes in one (not sure if that is a good idea) -- a parser
 * for a SAM text header, and a generator of SAM text header.
 */
public class SAMTextHeaderCodec {
    private static final String HEADER_LINE_START = "@";

    // These attributes are populated when parsing or generating
    private SAMFileHeader mFileHeader;

    // These attributes are populated when parsing text
    private String mCurrentLine;
    private LineReader mReader;
    private File mFile;
    private List<SAMSequenceRecord> sequences;
    private List<SAMReadGroupRecord> readGroups;

    // These attributes are populated when generating text
    private BufferedWriter writer;

    private static final String TAG_KEY_VALUE_SEPARATOR = ":";
    private static final String FIELD_SEPARATOR = "\t";

    public SAMTextHeaderCodec() {
    }

    /**
     * Reads text and converts to a SAMFileHeader object.  Note that one line past
     * the header must be read in order to determine the end of the header.  This line can be
     * obtained after parseTextHeader() has returned by calling getCurrentLine()
     * @param reader Where to get header text from.
     * @param file Name of the input file, for error messages.  May be null.
     * @return complete header object.
     */
    public SAMFileHeader decode(final LineReader reader, final File file) {
        mFileHeader = new SAMFileHeader();
        mReader = reader;
        mFile = file;
        sequences = new ArrayList<SAMSequenceRecord>();
        readGroups = new ArrayList<SAMReadGroupRecord>();

        while (advanceLine() != null) {
            if (!mCurrentLine.startsWith(HEADER_LINE_START)) {
                break;
            }
            final ParsedHeaderLine parsedHeaderLine = new ParsedHeaderLine(mCurrentLine);
            switch (parsedHeaderLine.getHeaderRecordType()) {

                case HD:
                    parseHDLine(parsedHeaderLine);
                    break;
                case PG:
                    parsePGLine(parsedHeaderLine);
                    break;
                case RG:
                    parseRGLine(parsedHeaderLine);
                    break;
                case SQ:
                    parseSQLine(parsedHeaderLine);
                    break;
                default:
                    throw new IllegalStateException("Unrecognized header record type: " +
                            parsedHeaderLine.getHeaderRecordType());
            }
        }
        mFileHeader.setSequences(sequences);
        mFileHeader.setReadGroups(readGroups);
        return mFileHeader;
    }

    private String advanceLine() {
        mCurrentLine = mReader.readLine();
        return mCurrentLine;
    }

    private void parsePGLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.PG.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMProgramRecord.PROGRAM_GROUP_ID_TAG);
        final SAMProgramRecord programRecord = new SAMProgramRecord(parsedHeaderLine.removeValue(SAMProgramRecord.PROGRAM_GROUP_ID_TAG));
        for (final Map.Entry<String, String> entry : parsedHeaderLine.mKeyValuePairs.entrySet()) {
            programRecord.setAttribute(entry.getKey(), entry.getValue());
        }
        mFileHeader.addProgramRecord(programRecord);
    }

    private void parseRGLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.RG.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMReadGroupRecord.READ_GROUP_ID_TAG);
        parsedHeaderLine.requireTag(SAMReadGroupRecord.READ_GROUP_SAMPLE_TAG);
        final SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord(parsedHeaderLine.removeValue(SAMReadGroupRecord.READ_GROUP_ID_TAG));
        for (final Map.Entry<String, String> entry : parsedHeaderLine.mKeyValuePairs.entrySet()) {
            samReadGroupRecord.setAttribute(entry.getKey(), entry.getValue());
        }

        // Convert non-String attributes to the appropriate types
        final String predictedMedianInsertSize =
                (String)samReadGroupRecord.getAttribute(SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG);
        if (predictedMedianInsertSize != null) {
            try {
                samReadGroupRecord.setAttribute(SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG,
                    Integer.parseInt(predictedMedianInsertSize));
            } catch (NumberFormatException e) {
                throw new SAMFormatException(SAMReadGroupRecord.PREDICTED_MEDIAN_INSERT_SIZE_TAG +
                        " is not numeric: " + predictedMedianInsertSize, e);
            }
        }

/*
TODO: Need an ISO 6801 date parser
        String dateRunProduced = (String)samReadGroupRecord.getAttribute(SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG);
        if (dateRunProduced != null) {
            try {
                Date date = dateParser.parse(dateRunProduced);
                samReadGroupRecord.setAttribute(SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG, date);
            } catch (ParseException e) {
                throw new SAMFormatException(SAMReadGroupRecord.DATE_RUN_PRODUCED_TAG + " cannot be parsed as a date: " +
                        dateRunProduced, e);
            }
        }
*/

        readGroups.add(samReadGroupRecord);
    }

    private void parseSQLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.SQ.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMSequenceRecord.SEQUENCE_NAME_TAG);
        parsedHeaderLine.requireTag(SAMSequenceRecord.SEQUENCE_LENGTH_TAG);
        final SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(parsedHeaderLine.removeValue(SAMSequenceRecord.SEQUENCE_NAME_TAG));
        samSequenceRecord.setSequenceLength(Integer.parseInt(parsedHeaderLine.removeValue(SAMSequenceRecord.SEQUENCE_LENGTH_TAG)));
        for (final Map.Entry<String, String> entry : parsedHeaderLine.mKeyValuePairs.entrySet()) {
            samSequenceRecord.setAttribute(entry.getKey(), entry.getValue());
        }
        sequences.add(samSequenceRecord);
    }

    private void parseHDLine(final ParsedHeaderLine parsedHeaderLine) {
        assert(HeaderRecordType.HD.equals(parsedHeaderLine.getHeaderRecordType()));
        parsedHeaderLine.requireTag(SAMFileHeader.VERSION_TAG);
        for (final Map.Entry<String, String> entry : parsedHeaderLine.mKeyValuePairs.entrySet()) {
            mFileHeader.setAttribute(entry.getKey(), entry.getValue());
        }
    }

    private RuntimeException reportErrorParsingLine(final String reason) {
        String fileMessage = "";
        if (mFile != null) {
            fileMessage = "File " + mFile + "; ";
        }
        return new SAMFormatException("Error parsing text SAM file. " + reason + "; " + fileMessage +
                "Line " + mReader.getLineNumber() + "\nLine: " + mCurrentLine);
    }

    private enum HeaderRecordType {
        HD, SQ, RG, PG
    }

    private class ParsedHeaderLine {
        private final HeaderRecordType mHeaderRecordType;
        private final Map<String, String> mKeyValuePairs = new HashMap<String, String>();

        ParsedHeaderLine(final String line) {
            assert(line.startsWith(HEADER_LINE_START));
            final String[] fields = line.split(FIELD_SEPARATOR);
            try {
                mHeaderRecordType = HeaderRecordType.valueOf(fields[0].substring(1));
            } catch (IllegalArgumentException e) {
                throw reportErrorParsingLine("Unrecognized header record type");
            }
            for (int i = 1; i < fields.length; ++i) {
                final String[] keyAndValue = fields[i].split(TAG_KEY_VALUE_SEPARATOR, 2);
                if (keyAndValue.length != 2) {
                    throw reportErrorParsingLine("Problem parsing " + HEADER_LINE_START + mHeaderRecordType +
                            " key:value pair");
                }
                mKeyValuePairs.put(keyAndValue[0], keyAndValue[1]);
            }
        }

        void requireTag(final String tag) {
            if (!mKeyValuePairs.containsKey(tag)) {
                throw reportErrorParsingLine(HEADER_LINE_START + mHeaderRecordType + " line missing " + tag + " tag");
            }
        }

        public HeaderRecordType getHeaderRecordType() {
            return mHeaderRecordType;
        }

        boolean containsKey(final String key) {
            return mKeyValuePairs.containsKey(key);
        }

        String getValue(final String key) {
            return mKeyValuePairs.get(key);
        }

        String removeValue(final String key) {
            final String ret = mKeyValuePairs.get(key);
            mKeyValuePairs.remove(key);
            return ret;
        }

    }

    /**
     * After parsing the text header, this object has gobbled one line too many.  Call this to get that line.
     * @return the first non-header line, or null if there isn't one.
     */
    public String getCurrentLine() {
        return mCurrentLine;
    }

    /**
     *
     * @param writer where to write the header text
     * @param header object to be converted to text.
     */
    public void encode(final Writer writer, final SAMFileHeader header) {
        mFileHeader = header;
        this.writer = new BufferedWriter(writer);
        writeHDLine();
        for (final SAMSequenceRecord sequenceRecord: header.getSequences()) {
            writeSQLine(sequenceRecord);
        }

        for (final SAMReadGroupRecord readGroup : header.getReadGroups()) {
            writeRGLine(readGroup);
        }
        for (final SAMProgramRecord programRecord : header.getProgramRecords()) {
            writePGLine(programRecord);
        }
        try {
            this.writer.flush();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private void println(final String s) {
        try {
            writer.append(s);
            writer.append("\n");
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private void writePGLine(SAMProgramRecord programRecord) {
        if (programRecord == null) {
            return;
        }
        final String[] fields = new String[2 + programRecord.getAttributes().size()];
        fields[0] = HEADER_LINE_START + HeaderRecordType.PG;
        fields[1] = SAMProgramRecord.PROGRAM_GROUP_ID_TAG + TAG_KEY_VALUE_SEPARATOR + programRecord.getProgramGroupId();
        int i = 2;
        for (final Map.Entry<String, String> entry: programRecord.getAttributes()) {
            fields[i++] = entry.getKey() + TAG_KEY_VALUE_SEPARATOR + entry.getValue();
        }
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    private void writeRGLine(final SAMReadGroupRecord readGroup) {
        final String[] fields = new String[2 + readGroup.getAttributes().size()];
        fields[0] = HEADER_LINE_START + HeaderRecordType.RG;
        fields[1] = SAMReadGroupRecord.READ_GROUP_ID_TAG + TAG_KEY_VALUE_SEPARATOR + readGroup.getReadGroupId();
        int i = 2;
        for (final Map.Entry<String, Object> entry: readGroup.getAttributes()) {
            fields[i++] = entry.getKey() + TAG_KEY_VALUE_SEPARATOR + entry.getValue().toString();
        }
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    private void writeHDLine() {
        final String[] fields = new String[1 + mFileHeader.getAttributes().size()];
        fields[0] = HEADER_LINE_START + HeaderRecordType.HD;
        int i = 1;
        for (final Map.Entry<String, Object> entry: mFileHeader.getAttributes()) {
            fields[i++] = entry.getKey() + TAG_KEY_VALUE_SEPARATOR + entry.getValue().toString();
        }
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

    private void writeSQLine(final SAMSequenceRecord sequenceRecord) {
        final int numAttributes =sequenceRecord.getAttributes() != null ? sequenceRecord.getAttributes().size() : 0;
        final String[] fields = new String[3 + numAttributes];
        fields[0] = HEADER_LINE_START + HeaderRecordType.SQ;
        fields[1] = SAMSequenceRecord.SEQUENCE_NAME_TAG + TAG_KEY_VALUE_SEPARATOR + sequenceRecord.getSequenceName();
        fields[2] = SAMSequenceRecord.SEQUENCE_LENGTH_TAG + TAG_KEY_VALUE_SEPARATOR + Integer.toString(sequenceRecord.getSequenceLength());
        int i = 3;
        if (sequenceRecord.getAttributes() != null) {
            for (final Map.Entry<String, Object> entry: sequenceRecord.getAttributes()) {
                fields[i++] = entry.getKey() + TAG_KEY_VALUE_SEPARATOR + entry.getValue().toString();
            }
        }
        println(StringUtil.join(FIELD_SEPARATOR, fields));
    }

}

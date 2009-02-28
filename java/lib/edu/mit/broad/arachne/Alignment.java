/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2008 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package edu.mit.broad.arachne;


/**
 * This class represents an arachne LookAlign alignment (or other related data structures).
 */
public class Alignment {

    private static final char TAB = '\t';

    private int mASequenceId;
    private int mASequenceLength;
    private int mAStart;
    private int mAEnd;
    private int mBSequenceId;
    private int mBSequenceLength;
    private int mBStart;
    private int mBEnd;
    private char mOrientation;
    private int[] mAlignmentBlocks;


    public Alignment() {
    }

    public int getASequenceId() {
        return mASequenceId;
    }

    public void setASequenceId(int value) {
        mASequenceId = value;
    }

    public int getASequenceLength() {
        return mASequenceLength;
    }

    public void setASequenceLength(int value) {
        mASequenceLength = value;
    }

    public int getAStart() {
        return mAStart;
    }

    public void setAStart(int value) {
        mAStart = value;
    }

    public int getAEnd() {
        return mAEnd;
    }

    public void setAEnd(int value) {
        mAEnd = value;
    }

    public int getBSequenceId() {
        return mBSequenceId;
    }

    public void setBSequenceId(int value) {
        mBSequenceId = value;
    }

    public int getBSequenceLength() {
        return mBSequenceLength;
    }

    public void setBSequenceLength(int value) {
        mBSequenceLength = value;
    }

    public int getBStart() {
        return mBStart;
    }

    public void setBStart(int value) {
        mBStart = value;
    }

    public int getBEnd() {
        return mBEnd;
    }

    public void setBEnd(int value) {
        mBEnd = value;
    }

    public char getOrientation() {
        return mOrientation;
    }

    public void setOrientation(char value) {
        mOrientation = value;
    }

    public int[] getAlignmentBlocks() {
        return mAlignmentBlocks;
    }

    public void setAlignmentBlocks(int[] value) {
        mAlignmentBlocks = value;
    }

    public static Alignment parse(String text) {

        if (text == null) {
            return null;
        }

        String[] fields = text.trim().split("\t");
        if (fields.length == 0) {
            return null;
        }

        if (!fields[0].equals("QUERY")) {
            throw new IllegalArgumentException("Invalid alignment: " + text);
        }
        if (fields.length < 14) {
            throw new IllegalArgumentException("Invalid alignment: " + text);
        }

        int seqAId = parseIntField(fields[1]);
        int seqAStart = parseIntField(fields[2]);
        int seqAEnd = parseIntField(fields[3]);
        int seqALength = parseIntField(fields[4]);
        int orientation = parseIntField(fields[5]);
        int seqBId = parseIntField(fields[6]);
        int seqBStart = parseIntField(fields[7]);
        int seqBEnd = parseIntField(fields[8]);
        int seqBLength = parseIntField(fields[9]);
        int blockCount = parseIntField(fields[10]);

        if (seqAStart < 0 || seqAEnd <= 0 || seqALength <= 0 ||
            seqAStart >= seqALength || seqAEnd > seqALength || seqAStart >= seqAEnd) {
            throw new IllegalArgumentException("Invalid alignment: " + text);
        }
        if (seqBStart < 0 || seqBEnd <= 0 || seqBLength <= 0 ||
            seqBStart >= seqBLength || seqBEnd > seqBLength || seqBStart >= seqBEnd) {
            throw new IllegalArgumentException("Invalid alignment: " + text);
        }
        if (orientation < 0 || orientation > 1) {
            throw new IllegalArgumentException("Invalid alignment: " + text);
        }
        if (fields.length != (11 + 3*blockCount)) {
            throw new IllegalArgumentException("Invalid alignment: " + text);
        }

        int[] alignmentBlocks = new int[3*blockCount];
        for (int i = 0; i < 3*blockCount; i++) {
            alignmentBlocks[i] = parseIntField(fields[11 + i]);
        }

        Alignment alignment = new Alignment();
        alignment.setASequenceId(seqAId);
        alignment.setASequenceLength(seqALength);
        alignment.setAStart(seqAStart+1);
        alignment.setAEnd(seqAEnd);
        alignment.setBSequenceId(seqBId);
        alignment.setBSequenceLength(seqBLength);
        alignment.setBStart(seqBStart+1);
        alignment.setBEnd(seqBEnd);
        alignment.setOrientation((orientation == 0) ? '+' : '-');
        alignment.setAlignmentBlocks(alignmentBlocks);
        return alignment;
    }

    private static int parseIntField(String text) {
        try {
            return Integer.parseInt(text);
        } catch (NumberFormatException exc) {
            throw new IllegalArgumentException("Illegal alignment field: " + text);
        }
    }

    public String arachneFormat() {
        StringBuilder builder = new StringBuilder();
        builder.append("QUERY");
        builder.append(TAB);
        builder.append(mASequenceId);
        builder.append(TAB);
        builder.append(mAStart-1); // zero based
        builder.append(TAB);
        builder.append(mAEnd);
        builder.append(TAB);
        builder.append(mASequenceLength);
        builder.append(TAB);
        builder.append(mOrientation == '+' ? 0 : 1);
        builder.append(TAB);
        builder.append(mBSequenceId);
        builder.append(TAB);
        builder.append(mBStart-1); // zero based
        builder.append(TAB);
        builder.append(mBEnd);
        builder.append(TAB);
        builder.append(mBSequenceLength);
        builder.append(TAB);
        builder.append(mAlignmentBlocks.length / 3);
        for (int i = 0; i < mAlignmentBlocks.length; i++) {
            builder.append(TAB);
            builder.append(mAlignmentBlocks[i]);
        }
        return builder.toString();
    }

    public String format() {
        StringBuilder builder = new StringBuilder();
        builder.append("Alignment");
        builder.append(' ');
        builder.append(mASequenceId);
        builder.append(' ');
        builder.append(mAStart);
        builder.append(' ');
        builder.append(mAEnd);
        builder.append(' ');
        builder.append(mOrientation);
        builder.append(' ');
        builder.append(mBSequenceId);
        builder.append(' ');
        builder.append(mBStart);
        builder.append(' ');
        builder.append(mBEnd);
        builder.append(' ');
        builder.append(mAlignmentBlocks.length / 3);
        for (int i = 0; i < mAlignmentBlocks.length; i++) {
            builder.append(' ');
            builder.append(mAlignmentBlocks[i]);
        }
        return builder.toString();
    }        
}

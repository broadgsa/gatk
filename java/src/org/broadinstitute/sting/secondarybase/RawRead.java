package org.broadinstitute.sting.secondarybase;

/**
 * RawRead represents lane and tile coordinates, raw intensities, read bases, and quality scores
 *
 * @author Kiran Garimella
 */
public class RawRead {
    private byte lane;
    private short tile;
    private short x;
    private short y;
    
    private byte[] sequence;
    private byte[] quals;
    private short[][] intensities;

    /**
     * Construct a raw read from the output of a PasteParser (in the order of int, seq, prb).
     * Takes data within specified cycle ranges.
     *
     * @param pastedReadString  the 3x(fragment length) output array from the PasteParser.
     * @param cycleBegin        the start cycle for the read (0-based, inclusive)
     * @param cycleEnd          the end cycle for the read (0-based, inclusive)
     */
    public RawRead(String[][] pastedReadString, int cycleBegin, int cycleEnd) {
        lane = Byte.valueOf(pastedReadString[0][0]);
        tile = Short.valueOf(pastedReadString[0][1]);
        x = Short.valueOf(pastedReadString[0][2]);
        y = Short.valueOf(pastedReadString[0][3]);

        sequence = pastedReadString[1][4].substring(cycleBegin, cycleEnd).getBytes();

        quals = new byte[sequence.length];
        intensities = new short[sequence.length][4];

        for (int cycle = 0; cycle < sequence.length; cycle++) {
            byte maxQual = -50;

            for (int fullReadIndex = 4*cycle; fullReadIndex < 4*cycle + 4; fullReadIndex++) {
                byte qual = Byte.valueOf(pastedReadString[2][fullReadIndex]);

                if (qual > maxQual) { maxQual = qual; }
            }

            quals[cycle] = maxQual >= 0 ? maxQual : 0;

            for (int fullReadIndex = 4*cycle + 4, channel = 0; fullReadIndex < 4*cycle + 8; fullReadIndex++, channel++) {
                double doubleChannelIntensity = Double.valueOf(pastedReadString[0][fullReadIndex]);
                short shortChannelIntensity = (short) doubleChannelIntensity;

                intensities[cycle][channel] = shortChannelIntensity;
            }
        }
    }

    /**
     * Get lane number of read
     *
     * @return lane number of read
     */
    public byte getLane() { return lane; }

    /**
     * Get tile number of read
     *
     * @return tile number of read
     */
    public int getTile() { return tile; }

    /**
     * Get x-coordinate of read
     *
     * @return x-coordinate of read
     */
    public int getXCoordinate() { return x; }

    /**
     * Get y-coordinate of read
     *
     * @return y-coordinate of read
     */
    public int getYCoordinate() { return y; }

    /**
     * Get read key (lane:tile:x:y)
     *
     * @return read key (lane:tile:x:y)
     */
    public String getReadKey() { return String.format("%d:%d:%d:%d", lane, tile, x, y); }

    /**
     * Get the read sequence between the cycles specified in the constructor as a byte array
     *
     * @return read sequence
     */
    public byte[] getSequence() { return sequence; }

    /**
     * Set the read sequence from a byte array
     *
     * @param sequence  the read sequence in byte array form
     */
    public void setSequence(byte[] sequence) { this.sequence = sequence; }

    /**
     * Get the read sequence as a string
     *
     * @return  the read sequence in string form
     */
    public String getSequenceAsString() {
        return new String(getSequence());
    }

    /**
     * Get the quals
     *
     * @return a byte array of quals
     */
    public byte[] getQuals() { return quals; }

    /**
     * Set the quals
     *
     * @param quals  a byte array of quals
     */
    public void setQuals(byte[] quals) { this.quals = quals; }

    /**
     * Get the raw read intensities
     *
     * @return the (readLength)x(numChannels) array of raw intensities
     */
    public short[][] getIntensities() { return intensities; }

    /**
     * Set the raw intensities
     *
     * @param intensities  the (readLength)x(numChannels) array of raw intensities
     */
    public void setIntensities(short[][] intensities) { this.intensities = intensities; }

    /**
     * Get the read length
     *
     * @return the read length
     */
    public int getReadLength() { return sequence.length; }
}

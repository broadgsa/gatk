package org.broadinstitute.sting.secondarybase;

/**
 * RawRead represents lane and tile coordinates, raw intensities, read bases, and quality scores
 *
 * @author Kiran Garimella
 */
public class RawRead implements Comparable<RawRead> {
    private byte lane;
    private short tile;
    private short x;
    private short y;
    
    private byte[] sequence;
    private byte[] quals;
    private short[][] intensities;

    /**
     * Blank constructor.
     */
    public RawRead() {}

    /**
     * Construct a raw read from the output of a PasteParser (in the order of int, seq, prb).
     * Takes data from entire read.
     *
     * @param pastedReadString  the 3x(fragment length) output array from the PasteParser.
     */
    public RawRead(String[][] pastedReadString) {
        loadRange(pastedReadString, 0, pastedReadString[1][4].length() - 1);
    }

    /**
     * Construct a raw read from the output of a PasteParser (in the order of int, seq, prb).
     * Takes data within specified cycle ranges.
     *
     * @param pastedReadString  the 3x(fragment length) output array from the PasteParser.
     * @param cycleBegin        the start cycle for the read (0-based, inclusive)
     * @param cycleEnd          the end cycle for the read (0-based, inclusive)
     */
    public RawRead(String[][] pastedReadString, int cycleBegin, int cycleEnd) {
        loadRange(pastedReadString, cycleBegin, cycleEnd);
    }

    /**
     * Does the actual parsing of the PasteParser output.
     *
     * @param pastedReadString  the 3x(fragment length) output array from the PasteParser.
     * @param cycleBegin        the start cycle for the read (0-based, inclusive)
     * @param cycleEnd          the end cycle for the read (0-based, inclusive)
     */
    private void loadRange(String pastedReadString[][], int cycleBegin, int cycleEnd) {
        lane = Byte.valueOf(pastedReadString[0][0]);
        tile = Short.valueOf(pastedReadString[0][1]);
        x = Short.valueOf(pastedReadString[0][2]);
        y = Short.valueOf(pastedReadString[0][3]);

        sequence = pastedReadString[1][4].substring(cycleBegin, cycleEnd + 1).getBytes();

        quals = new byte[sequence.length];
        //intensities = new short[sequence.length][4];
        intensities = new short[4][sequence.length];

        for (int cycle = cycleBegin, offset = 0; cycle <= cycleEnd; cycle++, offset++) {
            byte maxQual = -50;

            for (int fullReadIndex = 4*cycle; fullReadIndex < 4*cycle + 4; fullReadIndex++) {
                byte qual = Byte.valueOf(pastedReadString[2][fullReadIndex]);

                if (qual > maxQual) { maxQual = qual; }
            }

            quals[offset] = maxQual >= 0 ? maxQual : 0;

            for (int fullReadIndex = 4*cycle + 4, channel = 0; fullReadIndex < 4*cycle + 8; fullReadIndex++, channel++) {
                double doubleChannelIntensity = Double.valueOf(pastedReadString[0][fullReadIndex]);
                short shortChannelIntensity = (short) doubleChannelIntensity;

                //intensities[offset][channel] = shortChannelIntensity;
                intensities[channel][offset] = shortChannelIntensity;
            }
        }
    }

    /**
     * Get lane number of read.
     *
     * @return lane number of read
     */
    public byte getLane() { return lane; }

    /**
     * Set lane number of read.
     *
     * @param lane  lane number of read
     */
    public void setLane(byte lane) { this.lane = lane; }

    /**
     * Get tile number of read.
     *
     * @return tile number of read
     */
    public short getTile() { return tile; }

    /**
     * Set tile number of read.
     *
     * @param tile tile number of read
     */
    public void setTile(short tile) { this.tile = tile; }

    /**
     * Get x-coordinate of read.
     *
     * @return x-coordinate of read
     */
    public int getXCoordinate() { return x; }

    /**
     * Set x-coordinate of read.
     *
     * @param x  x-coordinate of read
     */
    public void setXCoordinate(short x) { this.x = x; }

    /**
     * Get y-coordinate of read.
     *
     * @return y-coordinate of read
     */
    public int getYCoordinate() { return y; }

    /**
     * Set y-coordinate of read.
     *
     * @param y  y-coordinate of read
     */
    public void setYCoordinate(short y) { this.y = y; }

    /**
     * Get read key (lane:tile:x:y).
     *
     * @return read key (lane:tile:x:y)
     */
    public String getReadKey() { return String.format("%d:%d:%d:%d", lane, tile, x, y); }

    /**
     * Get the read sequence between the cycles specified in the constructor as a byte array.
     *
     * @return read sequence
     */
    public byte[] getSequence() { return sequence; }

    /**
     * Set the read sequence from a byte array.
     *
     * @param sequence  the read sequence in byte array form
     */
    public void setSequence(byte[] sequence) { this.sequence = sequence; }

    /**
     * Get the read sequence as a String.
     *
     * @return  the read sequence in String form
     */
    public String getSequenceAsString() { return new String(getSequence()); }

    /**
     * Get the quals.
     *
     * @return a byte array of quals
     */
    public byte[] getQuals() { return quals; }

    /**
     * Set the quals.
     *
     * @param quals  a byte array of quals
     */
    public void setQuals(byte[] quals) { this.quals = quals; }

    /**
     * Get the raw read intensities.
     *
     * @return the (numChannels)x(readLength) array of raw intensities
     */
    public short[][] getIntensities() { return intensities; }

    /**
     * Set the raw intensities.
     *
     * @param intensities  the (numChannels)x(readLength) array of raw intensities
     */
    public void setIntensities(short[][] intensities) { this.intensities = intensities; }

    /**
     * Get the read length.
     *
     * @return the read length
     */
    public int getReadLength() { return sequence.length; }

    /**
     * Return the sum of the quality scores for this RawRead.
     *
     * @return  the sum of the quality scores
     */
    public int getQualityScoreSum() {
        int qualSum = 0;
        for ( byte qual : quals ) {
            qualSum  += (int) qual;
        }

        return qualSum;
    }

    /**
     * Compare two RawRead objects by summing their quality scores.  The one with lower aggregate quality is the "lesser" RawRead.
     *
     * @param rawRead  the other RawRead
     * @return -1, 0, or 1 if the RawRead on which this method is called is the lesser one, is equal to the comparison RawRead, or is greater than the comparison RawRead, respectively.
     */
    public int compareTo(RawRead rawRead) {
        int qualSum1 = this.getQualityScoreSum();
        int qualSum2 = rawRead.getQualityScoreSum();

        if      (qualSum1 < qualSum2) { return -1; }
        else if (qualSum1 > qualSum2) { return  1; }

        return 0;
    }

    public RawRead getSubset(int cycleStart, int cycleStop) {
        RawRead subRead = new RawRead();

        subRead.setLane(lane);
        subRead.setTile(tile);
        subRead.setXCoordinate(x);
        subRead.setYCoordinate(y);

        byte[] newSequence = new byte[cycleStop - cycleStart + 1];
        byte[] newQuals = new byte[cycleStop - cycleStart + 1];
        //short[][] newIntensities = new short[cycleStop - cycleStart + 1][4];
        short[][] newIntensities = new short[4][cycleStop - cycleStart + 1];

        for (int cycle = cycleStart, offset = 0; cycle <= cycleStop; cycle++, offset++) {
            newSequence[offset] = sequence[cycle];
            newQuals[offset] = quals[cycle];
            //newIntensities[offset] = intensities[cycle];

            for (int channel = 0; channel < 4; channel++) {
                newIntensities[channel][offset] = intensities[channel][cycle];
            }
        }

        subRead.setSequence(newSequence);
        subRead.setQuals(newQuals);
        subRead.setIntensities(newIntensities);

        return subRead;
    }
}

package org.broadinstitute.sting.secondarybase;

import org.broadinstitute.sting.utils.BaseUtils;

public class RawRead {
    private byte lane;
    private short tile;
    private short x;
    private short y;
    
    private byte[] sequence;
    private byte[] quals;
    private short[][] intensities;

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

    public byte getLane() { return lane; }
    public int getTile() { return tile; }
    public int getXCoordinate() { return x; }
    public int getYCoordinate() { return y; }

    public String getReadKey() { return String.format("%d:%d:%d:%d", lane, tile, x, y); }

    public byte[] getSequence() { return sequence; }
    public void setSequence(byte[] sequence) { this.sequence = sequence; }

    public String getSequenceAsString() {
        return new String(getSequence());
    }

    public byte[] getQuals() { return quals; }
    public void setQuals(byte[] quals) { this.quals = quals; }

    public short[][] getIntensities() { return intensities; }
    public void setIntensities(short[][] intensities) { this.intensities = intensities; }

    public int getReadLength() { return sequence.length; }
}

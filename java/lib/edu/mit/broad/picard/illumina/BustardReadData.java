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

/**
 * Holds all the Bustard-level data we need (so far) about an individual read.
 *
 * @author Kathleen Tibbetts
 */
public class BustardReadData {

    private static final String PADDING ="00000";

    final private String machineName;
    final private int runNumber;
    final private int laneNumber;
    final private int tileNumber;
    final private String firstReadSequence;
    final private String firstReadQualities;
    final private String secondReadSequence;
    final private String secondReadQualities;
    final private boolean pf;
    final private double intensities[][];
    final private int xCoordinate;
    final private int yCoordinate;
    private final SolexaQualityConverter converter = new SolexaQualityConverter();


    /**
     * Constructor that takes everything to populate this object
     * 
     * @param machineName
     * @param runNumber
     * @param laneNumber
     * @param tileNumber
     * @param firstReadSequence
     * @param firstReadQualities
     * @param secondReadSequence
     * @param secondReadQualities
     * @param pf
     * @param intensities
     * @param xCoordinate
     * @param yCoordinate
     */
    public BustardReadData(String machineName, int runNumber, int laneNumber, int tileNumber,
                           String firstReadSequence, String firstReadQualities,
                           String secondReadSequence, String secondReadQualities, 
                           boolean pf, double[][] intensities, int xCoordinate, int yCoordinate ) {

        this.machineName = machineName;
        this.runNumber = runNumber;
        this.laneNumber = laneNumber;
        this.tileNumber = tileNumber;
        this.firstReadSequence = firstReadSequence;
        this.firstReadQualities = firstReadQualities;
        this.secondReadSequence = secondReadSequence;
        this.secondReadQualities = secondReadQualities;
        this.pf = pf;
        this.intensities = intensities;
        this.xCoordinate = xCoordinate;
        this.yCoordinate = yCoordinate;
    }

    // TODO: Finalize read name -- ask Tim
    /**
     * Composes a name for this read from its values
     *
     * @return the read name
     */
    public String getReadName() {
        return this.machineName + ":" + this.laneNumber + ":" + this.tileNumber +
                ":" + this.xCoordinate + ":" + this.yCoordinate;
    }

    /**
     * Gets Phred-style qualitites for the first read
     *
     * @return  the String of qualities
     */
    public String getFirstReadPhredQualities() {
        return decodeSolexaQualitiesToPhred(getFirstReadQualities());
    }

    /**
     * Gets Phred-style qualitites for the second read
     *
     * @return  the String of qualities
     */
    public String getSecondReadPhredQualities() {
        return decodeSolexaQualitiesToPhred(getSecondReadQualities());
    }

    /**
     * Converts a string of Solexa qualities to a Phred-style quality String
     *
     * @param qualities the Solexa qualities to decode
     * @return  the String of Phred qualities
     */
    private String decodeSolexaQualitiesToPhred(String qualities) {
        StringBuilder sb = new StringBuilder();
        for (char c : qualities.toCharArray()) {
            // Quality char is phred score + 33
            sb.append((char)(converter.solexaToPhred((byte)c)+33));
        }
        return sb.toString();
    }

    public String getMachineName() { return machineName; }
    public int getRunNumber() { return runNumber; }
    public int getLaneNumber() { return laneNumber; }
    public int getTileNumber() { return tileNumber; }
    public String getFirstReadSequence() { return firstReadSequence; }
    public String getFirstReadQualities() { return firstReadQualities; }
    public String getSecondReadSequence() { return secondReadSequence; }
    public String getSecondReadQualities() { return secondReadQualities; }
    public double[][] getIntensities() { return intensities; }
    public boolean isPf() { return pf; }
    public int getXCoordinate() { return xCoordinate; }
    public int getYCoordinate() { return yCoordinate; }

}

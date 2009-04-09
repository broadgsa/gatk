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

/**
 * Holds all the Firecrest-level data we need (so far) about an individual read.
 *
 * @author Kiran Garimella
 */
public class FirecrestReadData {
    final private int laneNumber;
    final private int tileNumber;
    final private int xCoordinate;
    final private int yCoordinate;
    final private double[][] intensities;


    /**
     * Constructor that takes everything to populate this object
     * 
     * @param laneNumber
     * @param tileNumber
     * @param xCoordinate
     * @param yCoordinate
     * @param intensities
     */
    public FirecrestReadData(int laneNumber, int tileNumber, int xCoordinate, int yCoordinate, double[][] intensities) {
        this.laneNumber = laneNumber;
        this.tileNumber = tileNumber;
        this.xCoordinate = xCoordinate;
        this.yCoordinate = yCoordinate;
        this.intensities = intensities;
    }

    /**
     * Composes a name for this read from its values.
     *
     * @return the read name
     */
    public String getReadName() {
        return this.laneNumber + ":" + this.tileNumber + ":" + this.xCoordinate + ":" + this.yCoordinate + "#0";
    }

    public int getLaneNumber() { return laneNumber; }
    public int getTileNumber() { return tileNumber; }
    public int getXCoordinate() { return xCoordinate; }
    public int getYCoordinate() { return yCoordinate; }
    public double[][] getIntensities() { return intensities; }

}

/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.genotype;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Data object for Genotype Likelihoods. One object represents one row in a GELI file.
 *
 * @author Doug Voet
 */
public class GenotypeLikelihoods {
    /** this is a guess at how much memory an instance of this object occupies */
    public static final int OBJECT_SIZE_BYTES = 150;
    
    public static final int AA_GENOTYPE = 0;
    public static final int AC_GENOTYPE = 1;
    public static final int AG_GENOTYPE = 2;
    public static final int AT_GENOTYPE = 3;
    public static final int CC_GENOTYPE = 4;
    public static final int CG_GENOTYPE = 5;
    public static final int CT_GENOTYPE = 6;
    public static final int GG_GENOTYPE = 7;
    public static final int GT_GENOTYPE = 8;
    public static final int TT_GENOTYPE = 9;
    
    private static final char[][] GENOTYPES = {
        "AA".toCharArray(),
        "AC".toCharArray(),
        "AG".toCharArray(),
        "AT".toCharArray(),
        "CC".toCharArray(),
        "CG".toCharArray(),
        "CT".toCharArray(),
        "GG".toCharArray(),
        "GT".toCharArray(),
        "TT".toCharArray()
    };
    
    /** compares first by reference index then by position */
    public static class GenotypeLikelihoodsComparator implements Comparator<GenotypeLikelihoods> {
        @Override
        public int compare(GenotypeLikelihoods thing1, GenotypeLikelihoods thing2) {
            long refCompare = thing1.referenceIndex - thing2.referenceIndex;
            if (refCompare == 0) {
                long posCompare = thing1.position - thing2.position;
                return (int) posCompare;
            } else {
                return (int) refCompare;
            }
        }
    }


    private long referenceIndex;
    private long position;
    private byte referenceBase;
    private int numReads;
    private short maxMappingQuality;
    private float[] likelihoods = new float[10];
    private byte bestLikelihoodIndex = -1; // stored as byte to reduce memory footprint
    private byte secondBestLikelihoodIndex = -1; // stored as byte to reduce memory footprint
    
    public static int getLikelihoodIndex(char[] genotype) {
        char first = Character.isLowerCase(genotype[0]) ? Character.toUpperCase(genotype[0]) : genotype[0];
        char second = Character.isLowerCase(genotype[1]) ? Character.toUpperCase(genotype[1]) : genotype[1];
        if (first > second) {
            char temp = first;
            first = second;
            second = temp;
        }
        for (int i=0; i<GENOTYPES.length; i++) {
            if (first == GENOTYPES[i][0] && second == GENOTYPES[i][1]) {
                return i;
            }
        }
        throw new IllegalArgumentException("Unknown genotype string [" + new String(genotype) + 
                "], any pair of ACTG case insensitive is acceptable");
    }
    
    public float getLikelihood(int genotype) {
        return likelihoods[genotype];
    }
    
    public void setLikelihood(int genotype, float value) {
        likelihoods[genotype] = value;
    }
    
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("referenc ").append(referenceIndex).append(":").append(position);
        builder.append(", ref base ").append((char) referenceBase);
        builder.append(", #reads ").append(numReads);
        builder.append(", quality ").append(maxMappingQuality);
        builder.append(" [");
        for (int i=0; i<likelihoods.length; i++) {
            builder.append(GENOTYPES[i]).append(":").append(likelihoods[i]).append(" ");
        }
        builder.append("]");
        return builder.toString();
    }
    
    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + Arrays.hashCode(likelihoods);
        result = prime * result + maxMappingQuality;
        result = prime * result + numReads;
        result = prime * result + (int) (position ^ (position >>> 32));
        result = prime * result + referenceBase;
        result = prime * result + (int) (referenceIndex ^ (referenceIndex >>> 32));
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        GenotypeLikelihoods other = (GenotypeLikelihoods) obj;
        if (!Arrays.equals(likelihoods, other.likelihoods))
            return false;
        if (maxMappingQuality != other.maxMappingQuality)
            return false;
        if (numReads != other.numReads)
            return false;
        if (position != other.position)
            return false;
        if (referenceBase != other.referenceBase)
            return false;
        if (referenceIndex != other.referenceIndex)
            return false;
        return true;
    }

    public long getReferenceIndex() { return referenceIndex; }
    public void setReferenceIndex(long sequenceIndex) { this.referenceIndex = sequenceIndex; }
    public long getPosition() { return position; }
    public void setPosition(long position) { this.position = position; }
    public byte getReferenceBase() { return referenceBase; }
    public void setReferenceBase(byte referenceBase) { this.referenceBase = referenceBase; }
    public int getNumReads() { return numReads; }
    public void setNumReads(int numReads) { this.numReads = numReads; }
    public short getMaxMappingQuality() { return maxMappingQuality; }
    public void setMaxMappingQuality(short maxMappingQuality) { this.maxMappingQuality = maxMappingQuality; }
    float[] getLikelihoods() { return likelihoods; }
    public int getBestLikelihoodIndex() { return bestLikelihoodIndex; }
    public void setBestLikelihoodIndex(int bestLikelihoodIndex) { this.bestLikelihoodIndex = (byte) bestLikelihoodIndex; }
    public int getSecondBestLikelihoodIndex() { return secondBestLikelihoodIndex; }
    public void setSecondBestLikelihoodIndex(int secondBestLikelihoodIndex) { this.secondBestLikelihoodIndex = (byte) secondBestLikelihoodIndex; }
}

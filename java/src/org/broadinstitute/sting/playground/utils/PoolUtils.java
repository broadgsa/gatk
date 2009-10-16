package org.broadinstitute.sting.playground.utils;

import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;

import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.BaseUtils;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 27, 2009
 * Time: 12:31:08 PM
 * To change this template use File | Settings | File Templates.
 */
public class PoolUtils {

    private PoolUtils() {
    }

    public static final int BASE_A_OFFSET = 0;
    public static final int BASE_C_OFFSET = 1;
    public static final int BASE_G_OFFSET = 2;
    public static final int BASE_T_OFFSET = 3;
    public static final int BASE_INDEXED_ARRAY_SIZE = 4;

    public static Pair<List<SAMRecord>, List<Integer>> splitReadsByIndels( List<SAMRecord> reads, List<Integer> offsets, boolean returnBases ) {

        List<SAMRecord> baseReads = new ArrayList<SAMRecord>();
        List<SAMRecord> indelReads = new ArrayList<SAMRecord>();
        List<Integer> baseOffsets = new ArrayList<Integer>();
        List<Integer> indelOffsets = new ArrayList<Integer>();


        for ( int r = 0; r < reads.size(); r ++ ) {
            SAMRecord read = reads.get(r);
            int offset = offsets.get(r);
            if (BaseUtils.isRegularBase( (char) read.getReadBases()[offset] ) ) {
                baseReads.add(read);
                baseOffsets.add(offset);
            } else {
                indelReads.add(read);
                indelOffsets.add(offset);
            }
        }

        if (returnBases) {
            return new Pair<List<SAMRecord>,List<Integer>>(baseReads,baseOffsets);
        } else {
            return new Pair<List<SAMRecord>,List<Integer>>(indelReads,indelOffsets);
        }
    }

    public static ReadOffsetQuad splitReadsByReadDirection(List<SAMRecord> reads, List<Integer> offsets) {
        ArrayList<SAMRecord> forwardReads;
        ArrayList<SAMRecord> reverseReads;
        ArrayList<Integer> forwardOffsets;
        ArrayList<Integer> reverseOffsets;

        if (reads == null) {
            forwardReads = null;
            reverseReads = null;
            forwardOffsets = null;
            reverseOffsets = null;
        } else {
            forwardReads = new ArrayList();
            reverseReads = new ArrayList();
            forwardOffsets = new ArrayList();
            reverseOffsets = new ArrayList();

            for (int readNo = 0; readNo < reads.size(); readNo++) {
                if (reads.get(readNo).getReadNegativeStrandFlag()) {
                    forwardReads.add(reads.get(readNo));
                    forwardOffsets.add(offsets.get(readNo));
                } else {
                    reverseReads.add(reads.get(readNo));
                    reverseOffsets.add(offsets.get(readNo));
                }
            }
        }

        return new ReadOffsetQuad(forwardReads,forwardOffsets,reverseReads,reverseOffsets);
    }

    public static Pair<List<SAMRecord>[], List<Integer>[]> splitReadsByBase(List<SAMRecord> reads, List<Integer> offsets) {
        ArrayList<SAMRecord>[] readsByBase;
        ArrayList<Integer>[] offsetsByBase;
        if (reads == null) {
            readsByBase = null;
            offsetsByBase = null;
        } else {
            readsByBase = new ArrayList[4];
            offsetsByBase = new ArrayList[4];
            for (int readNum = 0; readNum < reads.size(); readNum++) {
                switch (reads.get(readNum).getReadBases()[offsets.get(readNum)]) {
                    case 'A':
                    case 'a':
                        readsByBase[BASE_A_OFFSET].add(reads.get(readNum));
                        offsetsByBase[BASE_A_OFFSET].add(offsets.get(readNum));
                        break;
                    case 'C':
                    case 'c':
                        readsByBase[BASE_C_OFFSET].add(reads.get(readNum));
                        offsetsByBase[BASE_C_OFFSET].add(offsets.get(readNum));
                        break;
                    case 'G':
                    case 'g':
                        readsByBase[BASE_G_OFFSET].add(reads.get(readNum));
                        offsetsByBase[BASE_G_OFFSET].add(offsets.get(readNum));
                        break;
                    case 'T':
                    case 't':
                        readsByBase[BASE_T_OFFSET].add(reads.get(readNum));
                        offsetsByBase[BASE_T_OFFSET].add(offsets.get(readNum));
                        break;
                    default:
                        break; // TODO: INDEL AWARENESS
                }
            }
        }
        return new Pair(readsByBase, offsetsByBase);
    }

    public static Pair<List<SAMRecord>, List<Integer>> thresholdReadsByQuality(List<SAMRecord> reads, List<Integer> offsets, byte qThresh) {
        List<SAMRecord> threshReads;
        List<Integer> threshOffsets;
        if (reads == null) {
            threshReads = null;
            threshOffsets = null;
        } else if (qThresh <= 0) {
            threshReads = reads;
            threshOffsets = offsets;
        } else {
            threshReads = new ArrayList();
            threshOffsets = new ArrayList();

            for (int readNo = 0; readNo < reads.size(); readNo++) {
                if (reads.get(readNo).getBaseQualities()[offsets.get(readNo)] >= qThresh) {
                    threshReads.add(reads.get(readNo));
                    threshOffsets.add(offsets.get(readNo));
                } // else do nothing
            }
        }

        return new Pair<List<SAMRecord>,List<Integer>>(threshReads, threshOffsets);
    }

    public static Pair<List<SAMRecord>,List<Integer>> thresholdReadsByMappingQuality( List<SAMRecord> reads, List<Integer> offsets, int mapQual ) {
        List<SAMRecord> goodMapReads;
        List<Integer> goodMapOffsets;
        if ( reads == null ) {
            goodMapReads = null;
            goodMapOffsets = null;
        } else if ( mapQual < 0 ) {
            goodMapReads = reads;
            goodMapOffsets = offsets;
        } else {
            goodMapReads = new ArrayList<SAMRecord>();
            goodMapOffsets = new ArrayList<Integer>();

            for ( int readNo = 0; readNo < reads.size(); readNo ++ ) {
                if ( reads.get(readNo).getMappingQuality() > mapQual ) {
                    goodMapReads.add(reads.get(readNo));
                    goodMapOffsets.add(offsets.get(readNo));
                }
            }
        }

        return new Pair<List<SAMRecord>,List<Integer>>(goodMapReads,goodMapOffsets);
    }

    public static Pair<List<SAMRecord>,List<Integer>> thresholdReadsByQuality(Pair<List<SAMRecord>,List<Integer>> readPair, byte qThresh) {
        return thresholdReadsByQuality(readPair.getFirst(),readPair.getSecond(),qThresh);
    }

    public static int getBaseOffset(char base) {
        switch (base) {
            case 'A':
            case 'a':
                return getBaseAOffset();
            case 'C':
            case 'c':
                return getBaseCOffset();
            case 'G':
            case 'g':
                return getBaseGOffset();
            case 'T':
            case 't':
                return getBaseTOffset();
            default:
                return -1;
        }
        //TODO: indel offsets
    }

    public static int getBaseAOffset() {
        return BASE_A_OFFSET;
    }

    public static int getBaseCOffset() {
        return BASE_C_OFFSET;
    }

    public static int getBaseGOffset() {
        return BASE_G_OFFSET;
    }

    public static int getBaseTOffset() {
        return BASE_T_OFFSET;
    }

    public static List<Byte> getReadBaseQualities(List<SAMRecord> reads, List<Integer> offsets) {
        List<Byte> qualities = new ArrayList<Byte>(reads.size());
        for (int readNo = 0; readNo < reads.size(); readNo++) {
            qualities.add(reads.get(readNo).getBaseQualities()[offsets.get(readNo)]);
        }

        return qualities;
    }

    public static double calculateLogLikelihoodOfSample(ReadOffsetQuad readsBySupport, double alleleFreq) {
        List<Byte> qSNP = getReadBaseQualities(readsBySupport.getFirstReads(), readsBySupport.getFirstOffsets());
        List<Byte> qRef = getReadBaseQualities(readsBySupport.getSecondReads(), readsBySupport.getSecondOffsets());
        Pair<Double,Double> logsumSNP = qListToSumLogProbabilities(true, qSNP, 1/alleleFreq);
        Pair<Double,Double> logsumRef = qListToSumLogProbabilities(false, qRef, 1/alleleFreq);
        return 0.0 - logsumSNP.getFirst() - logsumRef.getFirst() + logsumSNP.getSecond() + logsumRef.getSecond();
    }

    public static double calculateLogLikelihoodOfSample(Pair<Pair<List<SAMRecord>,List<SAMRecord>>,Pair<List<Integer>,List<Integer>>> snpReadsRefReads, int nIndivids) {
        List<Byte> qListSnps = getReadBaseQualities(snpReadsRefReads.getFirst().getFirst(),snpReadsRefReads.getSecond().getFirst());
        List<Byte> qListRefs = getReadBaseQualities(snpReadsRefReads.getFirst().getSecond(),snpReadsRefReads.getSecond().getSecond());
        Pair<Double,Double> logsumSNP = qListToSumLogProbabilities(true,qListSnps, 2.0*nIndivids);
        Pair<Double,Double> logsumRef = qListToSumLogProbabilities(false,qListRefs, 2.0*nIndivids);

        return 0.0 - logsumSNP.first - logsumRef.first + logsumSNP.second + logsumRef.second;
    }


    public static Pair<Double,Double> qListToSumLogProbabilities(boolean listRepresentsSNPObservations, List<Byte> qList, double denom)
    {
        double logProbObserveXAndSNPTrue = 0; // note "error" for SNP is observing a ref
        double logProbObserveXAndRefTrue = 0;// and "error" for ref is observing a SNP

        for (byte qual : qList) {
            double p_err = QualityUtils.qualToErrorProb(qual);
            if (listRepresentsSNPObservations) {
                logProbObserveXAndSNPTrue += Math.log10((1 - p_err) / denom +((denom - 1)*p_err) / denom);
                logProbObserveXAndRefTrue += Math.log10(p_err);
            } else {
                logProbObserveXAndSNPTrue += Math.log10((denom - 1) * (1 - p_err)/denom + p_err/denom);
                logProbObserveXAndRefTrue+= Math.log10(1 -p_err);
            }
        }

        return new Pair<Double,Double>(logProbObserveXAndSNPTrue,logProbObserveXAndRefTrue);
    }

    public static ReadOffsetQuad coinTossPartition(List<SAMRecord> reads, List<Integer> offsets, double alleleFreq) {
        ArrayList<SAMRecord> snpReads = new ArrayList<SAMRecord>();
        ArrayList<Integer> snpOffsets = new ArrayList<Integer>();
        ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();
        ArrayList<Integer> refOffsets = new ArrayList<Integer>();

        for( int i = 0; i < reads.size(); i ++ ) {
            if ( Math.random() < alleleFreq ) {
                snpReads.add(reads.get(i));
                snpOffsets.add(offsets.get(i));
            } else {
                refReads.add(reads.get(i));
                refOffsets.add(offsets.get(i));
            }
        }

        return new ReadOffsetQuad(snpReads,snpOffsets,refReads,refOffsets);
    }

}

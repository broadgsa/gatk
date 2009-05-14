package org.broadinstitute.sting.utils.duplicates;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

public class DupUtils {
    public static boolean usableDuplicate( SAMRecord read1, SAMRecord read2 ) {
            return read1 != read2 && read1.getReadLength() == read2.getReadLength();
    }

    public static Pair<SAMRecord, SAMRecord> combinedReadPair( List<SAMRecord> duplicateReads ) {
        if ( duplicateReads.size() < 4 )
            return null;

        SAMRecord c1 = combine2Duplicates(duplicateReads.get(0),duplicateReads.get(1));
        SAMRecord c2 = combine2Duplicates(duplicateReads.get(2),duplicateReads.get(3));
        return new Pair<SAMRecord, SAMRecord>(c1, c2);
    }

    public static SAMRecord sample3rdRead( List<SAMRecord> duplicateReads, SAMRecord read1, SAMRecord read2 ) {
         if ( duplicateReads.size() <= 2 ) {
             // no third unique read is available
             return null;
         } else {
             for ( SAMRecord read3 : duplicateReads ) {
                 if ( usableDuplicate(read1, read3) && usableDuplicate(read2, read3) )
                     return read3;
             }

             return null;
         }
    }

    public static SAMRecord tmpCopyRead(SAMRecord read) {
        SAMRecord copy = new SAMRecord(read.getHeader());
        copy.setReadName(read.getReadName());
        //copy.setReadString(final String value) {
        copy.setReadBases(read.getReadBases());
        copy.setBaseQualities(read.getBaseQualities());
        copy.setReferenceName(read.getReferenceName());
        copy.setReferenceIndex(read.getReferenceIndex());
        copy.setMateReferenceName(read.getMateReferenceName());
        copy.setMateReferenceIndex(read.getMateReferenceIndex());
        copy.setAlignmentStart(read.getAlignmentStart());
                //copy.setAlignmentEnd(read.getAlignmentEnd());
        copy.setMateAlignmentStart(read.getMateAlignmentStart());
        copy.setInferredInsertSize(read.getInferredInsertSize());
        copy.setMappingQuality(read.getMappingQuality());
        copy.setCigar(read.getCigar());
        copy.setFlags(read.getFlags());

        return copy;
    }

    public static SAMRecord combine2Duplicates(SAMRecord read1, SAMRecord read2) {
        byte[] read1Bases = read1.getReadBases();
        byte[] read1Quals = read1.getBaseQualities();
        byte[] read2Bases = read2.getReadBases();
        byte[] read2Quals = read2.getBaseQualities();

        byte[] bases = new byte[read1Bases.length];
        byte[] quals = new byte[read1Bases.length];

        SAMRecord c = tmpCopyRead(read1);
        for ( int i = 0; i < read1Bases.length; i++) {
            byte base1 = read1Bases[i];
            byte base2 = read2Bases[i];
            byte qual1 = read1Quals[i];
            byte qual2 = read2Quals[i];

            Pair<Byte, Integer> combined = combine2BasesAndQuals(base1, base2, qual1, qual2);

            bases[i] = combined.getFirst();
            quals[i] = QualityUtils.boundQual(combined.getSecond());

//            if ( DEBUG )
//                logger.debug(String.format("Combining %s (Q%2d) with %s (Q%2d) -> %s (Q%2d)%s%n",
//                    (char)base1, qual1, (char)base2, qual2, (char)bases[i], quals[i],
//                    base1 == base2 ? "" : " [MISMATCH]"));
        }
        c.setReadBases(bases);
        c.setBaseQualities(quals);

        return c;
    }

    public static Pair<Byte, Integer> combine2BasesAndQuals(byte base1, byte base2, int qual1, int qual2) {
        byte cbase;
        int cqual;
                
        if ( base1 == base2 ) {
            // agreement
            cbase = base1;
            cqual = qual1 + qual2;
        } else {
            // disagreement
            cbase = qual1 > qual2 ? base1 : base2;
            //cqual = Math.max(qual1, qual2);
            cqual = Math.max(qual1, qual2) - Math.min(qual1, qual2);

        }

        return new Pair<Byte, Integer>(cbase, cqual);
    }

    public static SAMRecord combineDuplicates(List<SAMRecord> duplicates, int maxQScore) {
        if ( duplicates.size() == 0 )
            return null;

        // make the combined read by copying the first read and setting the
        // bases and quals to new arrays
        SAMRecord comb = tmpCopyRead(duplicates.get(0));
        comb.setDuplicateReadFlag(false);
        int readLen = comb.getReadBases().length;
        byte[] bases = new byte[readLen];
        byte[] quals = new byte[readLen];

        for ( int i = 0; i < readLen; i++ ) {
            //System.out.printf("I is %d%n", i);
            //for ( SAMRecord read : duplicates ) {
            //    System.out.printf("dup base %c %d%n", (char)read.getReadBases()[i], read.getBaseQualities()[i]);
            //}
            Pair<Byte, Byte> baseAndQual = combineBaseProbs(duplicates, i, maxQScore);
            // baseAndQual = combineBasesByConsensus(duplicates, i);
            // baseAndQual = combineDupBasesAtI(duplicates, i);
            bases[i] = baseAndQual.getFirst();
            quals[i] = baseAndQual.getSecond();            
        }


        comb.setBaseQualities(quals);
        comb.setReadBases(bases);

        return comb;
    }

    private static Pair<Byte, Byte> baseProbs2BaseAndQual(double[] probs, int maxQScore) {
        char bestBase = 0;
        double bestProb = Double.NEGATIVE_INFINITY;
        double sumProbs = 0;
        for ( int i = 0; i < 4; i++ ) {
            sumProbs += Math.pow(10, probs[i]);
            //System.out.printf("Bestprob is %f > %f%n", bestProb, probs[i]);
            if ( probs[i] > bestProb ) {
                bestBase = BaseUtils.baseIndexToSimpleBase(i);
                bestProb = probs[i];
            }
        }
        Arrays.sort(probs);
        double normalizedP = Math.pow(10, bestProb) / sumProbs;
        double normalizedQ = 1 - normalizedP;
        double eps = Math.pow(10, -maxQScore/10.0);
        byte qual = QualityUtils.probToQual(normalizedP, eps);
        if ( false ) {
            System.out.printf("Best base is %s %.8f%n", bestBase, bestProb);
            System.out.printf("2nd  base is %.8f%n", probs[1]);
            System.out.printf("normalized P %.8f%n", normalizedP);
            System.out.printf("normalized Q %.8f%n", 1 - normalizedP);
            System.out.printf("max Q        %2d%n", maxQScore);
            System.out.printf("eps          %.8f%n", eps);
            System.out.printf("encoded    Q %2d%n", qual);
        }

        return new Pair<Byte, Byte>((byte)bestBase, qual);
    }

    private static void print4BaseQuals(String header, double[] probs) {
        System.out.printf("%s log10(P(b)) is ", header);
        for ( int i = 0; i < 4; i++ ) {
            System.out.printf("%c=%+.8f ", BaseUtils.baseIndexToSimpleBase(i), probs[i]);
        }
        System.out.printf("%n");
    }

    private static Pair<Byte, Byte> combineBaseProbs(List<SAMRecord> duplicates, int readOffset, int maxQScore) {
        List<Integer> offsets = BasicPileup.constantOffset(duplicates, readOffset);
        ArrayList<Byte> bases = BasicPileup.basePileup(duplicates, offsets);
        ArrayList<Byte> quals = BasicPileup.qualPileup(duplicates, offsets);
        final boolean debug = false;

        // calculate base probs
        double[] qualSums = {0.0, 0.0, 0.0, 0.0};
        if ( debug ) print4BaseQuals("start", qualSums);
        for ( int i = 0; i < bases.size(); i++ ) {
            char base = (char)(byte)bases.get(i);
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            byte qual = quals.get(i);
            double pqual = QualityUtils.qualToProb(qual);
            for ( int j = 0; j < 4; j++) {
                qualSums[j] += Math.log10(j == baseIndex ?  pqual : (1 - pqual)/3);
            }
            if ( debug ) print4BaseQuals(String.format("%c Q%2d", base, qual), qualSums);
        }
        if ( debug ) print4BaseQuals("final", qualSums);

        Pair<Byte, Byte> combined = baseProbs2BaseAndQual(qualSums, maxQScore);
        if ( debug ) 
            System.out.printf("%s %s => %c Q%s%n",
                    BasicPileup.basePileupAsString(duplicates, offsets),
                    BasicPileup.qualPileupAsString(duplicates, offsets),
                    (char)(byte)combined.getFirst(), combined.getSecond());
        return combined;
    }


    private static Pair<Byte, Byte> combineDupBasesAtI(List<SAMRecord> duplicates, int readOffset) {
        //System.out.printf("combineDupBasesAtI() dups=%s, readOffset=%d%n", duplicates.size(), readOffset);
        List<Integer> offsets = BasicPileup.constantOffset(duplicates, readOffset);
        ArrayList<Byte> bases = BasicPileup.basePileup(duplicates, offsets);
        ArrayList<Byte> quals = BasicPileup.qualPileup(duplicates, offsets);

        // start the recursion
        byte combinedBase = bases.get(0);
        int combinedQual  = quals.get(0);
        for ( int i = 1; i < bases.size(); i++ ) {
            // recursively combine evidence
            byte nextBase = bases.get(i);
            byte nextQual = quals.get(i);
            Pair<Byte, Integer> combinedI = combine2BasesAndQuals(combinedBase, nextBase, combinedQual, nextQual);

            if ( false )
                System.out.printf("Combining at %d: %s (Q%2d) with %s (Q%2d) -> %s (Q%2d)%s%n",
                        i,
                        (char)combinedBase, combinedQual,
                        (char)nextBase, nextQual,
                        (char)((byte)combinedI.getFirst()), combinedI.getSecond(),
                        combinedBase == nextBase ? "" : " [MISMATCH]");
            combinedBase = combinedI.getFirst();
            combinedQual = combinedI.getSecond();
        }

        if ( false )
            System.out.printf("%s %s => %c Q%s => Q%d%n",
                    BasicPileup.basePileupAsString(duplicates, offsets),
                    BasicPileup.qualPileupAsString(duplicates, offsets),
                    (char)combinedBase, combinedQual, QualityUtils.boundQual(combinedQual));
        return new Pair<Byte, Byte>(combinedBase, QualityUtils.boundQual(combinedQual));
    }

    private static Pair<Byte, Byte> combineBasesByConsensus(List<SAMRecord> duplicates, int readOffset) {
        //System.out.printf("combineDupBasesAtI() dups=%s, readOffset=%d%n", duplicates.size(), readOffset);
        List<Integer> offsets = BasicPileup.constantOffset(duplicates, readOffset);
        String bases = BasicPileup.basePileupAsString(duplicates, offsets);
        ArrayList<Byte> quals = BasicPileup.qualPileup(duplicates, offsets);
        byte combinedBase  = BasicPileup.consensusBase(bases);
        int baseMatches = Utils.countOccurances((char)combinedBase, bases);
        double maxP = QualityUtils.qualToProb(Utils.listMaxByte(quals));
        double mismatchRate = (double)baseMatches / bases.length();
        byte combinedQual = QualityUtils.probToQual(Math.min(maxP, mismatchRate), 0.0);

        if ( false )
            System.out.printf("%s %s => %c Q%s => Q%d%n",
                    BasicPileup.basePileupAsString(duplicates, offsets),
                    BasicPileup.qualPileupAsString(duplicates, offsets),
                    (char)combinedBase, combinedQual, QualityUtils.boundQual(combinedQual));
        
        return new Pair<Byte, Byte>(combinedBase, QualityUtils.boundQual(combinedQual));
    }
}
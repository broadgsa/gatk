package org.broadinstitute.sting.utils;

import net.sf.samtools.*;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 8:54:05 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class BasicPileup implements Pileup {

    public static final char DELETION_CHAR = 'D';

    protected boolean includeDeletions = false;

    public void setIncludeDeletionsInPileupString(boolean value) {
	includeDeletions = value;
    }

    public String getPileupString()
    {
        return String.format("%s: %s %s %s", getLocation(), getRef(), getBases(), getQuals());
    }

    public static List<Integer> constantOffset( List<SAMRecord> reads, int i ) {
        List<Integer> l = new ArrayList<Integer>(reads.size());
        for ( SAMRecord read : reads ) {
            l.add(i);
        }
        return l;
    }

    public static String basePileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
	return basePileupAsString( reads, offsets, false );
    }

    public static String basePileupAsString( List<SAMRecord> reads, List<Integer> offsets, boolean includeDeletions ) {
        StringBuilder bases = new StringBuilder();
        for ( byte base : basePileup(reads, offsets, includeDeletions)) {
            bases.append((char)base);
        }
        return bases.toString();
    }

    public static String baseWithStrandPileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
	return baseWithStrandPileupAsString( reads, offsets, false );
    }

    public static String baseWithStrandPileupAsString( List<SAMRecord> reads, List<Integer> offsets, boolean includeDeletions ) {
        StringBuilder bases = new StringBuilder();

        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

	    char base;
	    if ( offset == -1 ) {
		if ( includeDeletions )
		    base = DELETION_CHAR;
		else
		    continue;
	    } else {
		base = (char) read.getReadBases()[offset];
	    }

            base = Character.toUpperCase(base);
            if (read.getReadNegativeStrandFlag()) {
                base = Character.toLowerCase(base);
            }

            bases.append(base);
        }

        return bases.toString();
    }

    public static ArrayList<Byte> basePileup( List<SAMRecord> reads, List<Integer> offsets ) {
	return basePileup( reads, offsets, false );
    }

    public static ArrayList<Byte> basePileup( List<SAMRecord> reads, List<Integer> offsets, boolean includeDeletions ) {
        ArrayList<Byte> bases = new ArrayList<Byte>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
	    if ( offset == -1 ) {
		if ( includeDeletions )
		    bases.add((byte)DELETION_CHAR);
	    } else {
		bases.add(read.getReadBases()[offset]);
	    }
        }
        return bases;
    }

    public static ArrayList<Byte> qualPileup( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> quals = new ArrayList<Byte>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
	    // skip deletion sites
	    if ( offset == -1 )
		continue;
            byte qual = (byte)read.getBaseQualities()[offset];
            quals.add(qual);
        }
        return quals;
    }

    public static ArrayList<Byte> mappingQualPileup( List<SAMRecord> reads) {
        ArrayList<Byte> quals = new ArrayList<Byte>(reads.size());
        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            byte qual = (byte)read.getMappingQuality();
            quals.add(qual);
        }
        return quals;
    }

    public static String mappingQualPileupAsString( List<SAMRecord> reads) {
        return quals2String(mappingQualPileup(reads));
    }

    public static String quals2String( List<Byte> quals ) {
        StringBuilder qualStr = new StringBuilder();
        for ( int qual : quals ) {
            qual = Math.min(qual, 63);              // todo: fixme, this isn't a good idea
            char qualChar = (char) (33 + qual);     // todo: warning, this is illegal for qual > 63
            qualStr.append(qualChar);
        }

        return qualStr.toString();
    }

    public static String qualPileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        return quals2String(qualPileup(reads, offsets));
    }

    public static ArrayList<Byte> secondaryBasePileup( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> bases2 = new ArrayList<Byte>(reads.size());
        boolean hasAtLeastOneSQField = false;

        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            byte[] compressedQuals = (byte[]) read.getAttribute("SQ");
            byte base2;
            
            if (offset != -1 && compressedQuals != null && compressedQuals.length == read.getReadLength()) {
                base2 = (byte) BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(compressedQuals[offset]));
                hasAtLeastOneSQField = true;
            } else {
                base2 = (byte) '.';
            }
            bases2.add(base2);
        }
        return (hasAtLeastOneSQField ? bases2 : null);
    }

    public static String secondaryBasePileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        StringBuilder bases2 = new StringBuilder();
        ArrayList<Byte> sbases = secondaryBasePileup(reads, offsets);

        if (sbases == null) { return null; }

        ArrayList<Byte> pbases = basePileup(reads, offsets);

        Random generator = new Random();
        
        for (int pileupIndex = 0; pileupIndex < sbases.size(); pileupIndex++) {
            byte pbase = pbases.get(pileupIndex);
            byte sbase = sbases.get(pileupIndex);

            while (sbase == pbase) {
                sbase = (byte) BaseUtils.baseIndexToSimpleBase(generator.nextInt(4));
            }

            bases2.append((char) sbase);
        }

        /*
        for (byte base2 : secondaryBasePileup(reads, offsets)) {
            bases2.append((char) base2);
        }
        */

        return bases2.toString();
    }

    public static ArrayList<Byte> secondaryQualPileup( List<SAMRecord> reads, List<Integer> offsets ) {
        ArrayList<Byte> quals2 = new ArrayList<Byte>(reads.size());
        boolean hasAtLeastOneSQField = false;

        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            byte[] compressedQuals = (byte[]) read.getAttribute("SQ");
            byte qual2;
            if (offset != -1 && compressedQuals != null) {
                qual2 = QualityUtils.probToQual(QualityUtils.compressedQualityToProb(compressedQuals[offset]));
                hasAtLeastOneSQField = true;
            } else {
                qual2 = 0;
            }
            quals2.add(qual2);
        }
        return (hasAtLeastOneSQField ? quals2 : null);
    }

    public static String secondaryQualPileupAsString( List<SAMRecord> reads, List<Integer> offsets) {
        StringBuilder quals2 = new StringBuilder();
        ArrayList<Byte> sqquals = secondaryQualPileup(reads, offsets);

        if (sqquals == null) {
            return null;
        } else {
            for (byte qual2 : secondaryQualPileup(reads, offsets)) {
                quals2.append(qual2);
            }
            return quals2.toString();
        }
    }

    public static double[][] probDistPileup( List<SAMRecord> reads, List<Integer> offsets ) {
        double[][] dist = new double[reads.size()][4];

        for (int readIndex = 0; readIndex < dist.length; readIndex++) {
            SAMRecord read = reads.get(readIndex);

            String bases = read.getReadString();
            int offset = offsets.get(readIndex);
	    if ( offset == -1 )
		continue;

            int bestBaseIndex = BaseUtils.simpleBaseToBaseIndex(bases.charAt(offset));

            if (bestBaseIndex >= 0 && bestBaseIndex < 4) {
                dist[readIndex][bestBaseIndex] = QualityUtils.qualToProb(read.getBaseQualities()[offset]);

                byte[] sqs = (byte[]) read.getAttribute("SQ");
                if (sqs != null && QualityUtils.compressedQualityToBaseIndex(sqs[offset]) != bestBaseIndex) {
                    double epsilon = 1e-4;

                    int secondBestBaseIndex = QualityUtils.compressedQualityToBaseIndex(sqs[offset]);
                    //dist[readIndex][secondBestBaseIndex] = (1.0 - dist[readIndex][bestBaseIndex] - 2.0*epsilon);
                    dist[readIndex][secondBestBaseIndex] = 0.8*(1.0 - dist[readIndex][bestBaseIndex]);

                    for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                        if (baseIndex != bestBaseIndex && baseIndex != secondBestBaseIndex) {
                            //dist[readIndex][baseIndex] = epsilon;
                            dist[readIndex][baseIndex] = 0.1*(1.0 - dist[readIndex][bestBaseIndex]);
                        }
                    }
                } else {
                    for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                        if (baseIndex != bestBaseIndex) {
                            dist[readIndex][baseIndex] = (1.0 - dist[readIndex][bestBaseIndex])/3.0;
                        }
                    }
                }
            } else {
                for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                    dist[readIndex][baseIndex] = 0.25;
                }
            }
        }

        return dist;
    }
    
    public static String probDistPileupAsString( List<SAMRecord> reads, List<Integer> offsets ) {
        double[][] dist = probDistPileup(reads, offsets);

        String distString = String.format("     %c      %c      %c      %c\n", 'A', 'C', 'G', 'T');
        for (int readIndex = 0; readIndex < dist.length; readIndex++) {
            distString += "[ ";
            for (int baseIndex = 0; baseIndex < 4; baseIndex++) {
                distString += String.format("%4.4f ", dist[readIndex][baseIndex]);
            }
            distString += "]\n";
        }

        return distString;
    }

    public static String[] indelPileup( List<SAMRecord> reads, List<Integer> offsets ) 
    {
        String[] indels = new String[reads.size()];

        for (int i = 0; i < reads.size(); i++)
        {
            SAMRecord read = reads.get(i);
            Cigar cigar    = read.getCigar();
	    int offset     = offsets.get(i);

            String cigar_string = read.getCigarString();	 
			if (! (cigar_string.contains("I") || cigar_string.contains("D"))) { indels[i] = "null"; continue; }

			//System.out.printf("%s:%d %s %s %s ", read.getReferenceName(), read.getAlignmentStart(), read.getReadName(), read.getReadString(), cigar_string);
            int k = 0;
            for (int j = 0; j < cigar.numCigarElements(); j++)
            {
                CigarOperator operator = cigar.getCigarElement(j).getOperator();
                int           length   = cigar.getCigarElement(j).getLength();
                if (operator == CigarOperator.M) 
                { 
                    k += length; 
                }
                else if ((k == offset+1) && (operator == CigarOperator.I))
                {
                    // this insertion is associated with this offset (kinda ;) ).
  					indels[i] = read.getReadString().substring(k, k+length);	 
					//System.out.printf("(I,%d,%d)", k, offset);
					break;
                }
                else if ((k != offset+1) && (operator == CigarOperator.I)) 
                { 
					//System.out.printf("(i,%d,%d)", k, offset);
                    k += length; 
                }
                else if ((k == offset) && (operator == CigarOperator.D))
                {
					// this deletion is associated with this offset.
  					indels[i] = length + "D";
					//System.out.printf("(D,%d,%d)", k, offset);
   					break;
                }
                else if (k >= offset) 
                {
                    // no indel here.
                    indels[i] = "null";
					//System.out.printf("(N,%d,%d)", k, offset);
                    break;
                }
            }
			if (indels[i] == null) { indels[i] = "null"; }
			//System.out.printf("\n");
        }

        return indels;
    }

    public static String pileupDiff(final Pileup a, final Pileup b)
    {
        return pileupDiff(a,b,true);
    }

    private static String maybeSorted( final String x, boolean sortMe )
    {
        if ( sortMe ) {
            byte[] bytes = x.getBytes();
            Arrays.sort(bytes);
            return new String(bytes);
        }
        else
            return x;
    }

    public static String pileupDiff(final Pileup a, final Pileup b, boolean orderDependent)
    {
        if ( a.size() != b.size() )
            return "Sizes not equal";
        if ( a.getLocation().compareTo(b.getLocation()) != 0 )
            return "Locations not equal";

        String aBases = maybeSorted(a.getBases(), ! orderDependent );
        String bBases = maybeSorted(b.getBases(), ! orderDependent );
        if ( ! aBases.toUpperCase().equals(bBases.toUpperCase()) )
            return "Bases not equal";

        String aQuals = maybeSorted(a.getQuals(), ! orderDependent );
        String bQuals = maybeSorted(b.getQuals(), ! orderDependent );
        if ( ! aQuals.equals(bQuals) )
            return "Quals not equal";

        return null;
    }

    public static class BaseCounts {
        int a, c, t, g;

        public BaseCounts(int a, int c, int t, int g) {
            this.a = a;
            this.c = c;
            this.t = t;
            this.g = g;
        }
    }

    public static int countBase(final char base, final String bases) {
        return Utils.countOccurrences(base, bases);
    }

    public static BaseCounts countBases(final String bases) {
        String canon = bases.toUpperCase();
        return new BaseCounts(Utils.countOccurrences('A', canon),
                Utils.countOccurrences('C', canon),
                Utils.countOccurrences('T', canon),
                Utils.countOccurrences('G', canon));
    }

    public static byte consensusBase(String bases) {
        BaseCounts counts = countBases( bases );
        int ACount = counts.a;
        int CCount = counts.c;
        int TCount = counts.t;
        int GCount = counts.g;

        int m = Math.max(ACount, Math.max(CCount, Math.max(TCount, GCount)));
        if ( ACount == m ) return 'A';
        if ( CCount == m ) return 'C';
        if ( TCount == m ) return 'T';
        if ( GCount == m ) return 'G';
        return 0;
    }
}



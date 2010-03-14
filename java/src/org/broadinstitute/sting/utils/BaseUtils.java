package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMRecord;

import java.util.Random;

/**
 * BaseUtils contains some basic utilities for manipulating nucleotides.
 */
public class BaseUtils {
    public final static char[] BASES = { 'A', 'C', 'G', 'T' };
    public final static char[] EXTENDED_BASES = { 'A', 'C', 'G', 'T', 'N', 'D' };
    // todo -- fix me (enums?)
    public static final byte DELETION_INDEX = 4;
    public static final byte NO_CALL_INDEX = 5; // (this is 'N')

    public static int gIndex = BaseUtils.simpleBaseToBaseIndex('G');
    public static int cIndex = BaseUtils.simpleBaseToBaseIndex('C');
    public static int aIndex = BaseUtils.simpleBaseToBaseIndex('A');
    public static int tIndex = BaseUtils.simpleBaseToBaseIndex('T');


    /// In genetics, a transition is a mutation changing a purine to another purine nucleotide (A <-> G) or
    // a pyrimidine to another pyrimidine nucleotide (C <-> T).
    // Approximately two out of every three single nucleotide polymorphisms (SNPs) are transitions.
    public enum BaseSubstitutionType {
        TRANSITION,         // A <-> G or C <-> T
        TRANSVERSION
    }

    /**
     * Returns the base substitution type of the 2 state SNP
     * @param base1
     * @param base2
     * @return
     */
    public static BaseSubstitutionType SNPSubstitutionType( byte base1, byte base2 ) {
        BaseSubstitutionType t = isTransition(base1, base2) ? BaseSubstitutionType.TRANSITION : BaseSubstitutionType.TRANSVERSION;
        //System.out.printf("SNPSubstitutionType( char %c, char %c ) => %s%n", base1, base2, t);
        return t;
    }

    public static boolean isTransition( byte  base1, byte  base2 ) {
        int b1 = simpleBaseToBaseIndex(base1);
        int b2 = simpleBaseToBaseIndex(base2);
        return b1 == 0 && b2 == 2 || b1 == 2 && b2 == 0 ||
               b1 == 1 && b2 == 3 || b1 == 3 && b2 == 1;
    }

    public static boolean isTransversion( byte  base1, byte  base2 ) {
        return ! isTransition(base1, base2);
    }

    /** Private constructor.  No instantiating this class! */
    private BaseUtils() {}

    static public boolean basesAreEqual(byte base1, byte base2) {
        return simpleBaseToBaseIndex((char)base1) == simpleBaseToBaseIndex((char)base2);
    }


    /**
     * Converts a IUPAC nucleotide code to a pair of bases
     *
     * @param code
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    static public char[] iupacToBases(char code) {
        char[] bases = new char[2];
        switch (code) {
            case '*':               // the wildcard character counts as an A
            case 'A':
            case 'a':
                bases[0] = bases[1] = 'A';
                break;
            case 'C':
            case 'c':
                bases[0] = bases[1] = 'C';
                break;
            case 'G':
            case 'g':
                bases[0] = bases[1] = 'G';
                break;
            case 'T':
            case 't':
                bases[0] = bases[1] = 'T';
                break;
            case 'R':
            case 'r':
                bases[0] = 'A';
                bases[1] = 'G';
                break;
            case 'Y':
            case 'y':
                bases[0] = 'C';
                bases[1] = 'T';
                break;
            case 'S':
            case 's':
                bases[0] = 'G';
                bases[1] = 'C';
                break;
            case 'W':
            case 'w':
                bases[0] = 'A';
                bases[1] = 'T';
                break;
            case 'K':
            case 'k':
                bases[0] = 'G';
                bases[1] = 'T';
                break;
            case 'M':
            case 'm':
                bases[0] = 'A';
                bases[1] = 'C';
                break;
            default:
                bases[0] = bases[1] = 'N';
        }
        return bases;
    }
    /**
     * Converts a simple base to a base index
     *
     * @param base  [AaCcGgTt]
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    static public int simpleBaseToBaseIndex(char base) {
        switch (base) {
            case '*':               // the wildcard character counts as an A
            case 'A':
            case 'a': return 0;

            case 'C':
            case 'c': return 1;

            case 'G':
            case 'g': return 2;

            case 'T':
            case 't': return 3;

            default: return -1;
        }
    }

    static public int extendedBaseToBaseIndex(char base) {
        switch (base) {
            case 'd':
            case 'D': return DELETION_INDEX;
            case 'n':
            case 'N': return NO_CALL_INDEX;

            default: return simpleBaseToBaseIndex(base);
        }
    }

    static public int simpleBaseToBaseIndex(byte base) {
        return simpleBaseToBaseIndex((char)base);
    }

    static public boolean isRegularBase(char base) {
        return simpleBaseToBaseIndex(base) != -1;
    }

    static public boolean isRegularBase(byte base) {
        return isRegularBase((char)base);
    }

    static public boolean isNBase(byte base) {
        return base == 'N';
    }


    /**
     * Converts a base index to a simple base
     *
     * @param baseIndex  0, 1, 2, 3
     * @return A, C, G, T, or '.' if the index can't be understood
     */
    static public char baseIndexToSimpleBase(int baseIndex) {
        switch (baseIndex) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return '.';
        }
    }

    /**
     * Converts a base index to a base index representing its cross-talk partner
     *
     * @param baseIndex  0, 1, 2, 3
     * @return 1, 0, 3, 2, or -1 if the index can't be understood
     */
    static public int crossTalkPartnerIndex(int baseIndex) {
        switch (baseIndex) {
            case 0: return 1; // A -> C
            case 1: return 0; // C -> A
            case 2: return 3; // G -> T
            case 3: return 2; // T -> G
            default: return -1;
        }
    }

    /**
     * Converts a base to the base representing its cross-talk partner
     *
     * @param base  [AaCcGgTt]
     * @return C, A, T, G, or '.' if the base can't be understood
     */
    static public char crossTalkPartnerBase(char base) {
        return baseIndexToSimpleBase(crossTalkPartnerIndex(simpleBaseToBaseIndex(base)));
    }

    /**
     * Return the complement of a base index.
     *
     * @param baseIndex  the base index (0:A, 1:C, 2:G, 3:T)
     * @return the complementary base index
     */
    static public byte complementIndex(int baseIndex) {
        switch (baseIndex) {
            case 0: return 3; // a -> t
            case 1: return 2; // c -> g
            case 2: return 1; // g -> c
            case 3: return 0; // t -> a
            default: return -1; // wtf?
        }
    }


    public static byte getSecondBase(final SAMRecord read, int offset) {
        byte base2 = '.'; // todo -- what should the default char really be?

        if (read.getAttribute("SQ") != null) {
            byte[] compressedQuals = (byte[]) read.getAttribute("SQ");

            if (offset != -1 && compressedQuals != null && compressedQuals.length == read.getReadLength()) {
                base2 = (byte) BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex(compressedQuals[offset]));
            }
        }
        else if (read.getAttribute("E2") != null) {
            String secondaries = (String) read.getAttribute("E2");
            if (offset != -1 && secondaries != null && secondaries.length() == read.getReadLength()) {
                base2 = (byte)secondaries.charAt(offset);
            }
        }
        else {
            base2 = 'N';
        }

        return base2;
    }

    /**
     * Perform a transition (A <-> G or C <-> T) on the base, or the specified base if it can't be done (i.e. an ambiguous base).
     *
     * @param base the base [AaCcGgTt]
     * @return the transition of the base, or the input base if it's not one of the understood ones
     */
    static public char transition(char base) {
        switch (base) {
            case 'A':
            case 'a': return 'G';
            case 'C':
            case 'c': return 'T';
            case 'G':
            case 'g': return 'A';
            case 'T':
            case 't': return 'C';
            default: return base;
        }
    }

    /**
     * Perform a transversion (A <-> C or G <-> T) on the base, or the specified base if it can't be done (i.e. an ambiguous base).
     *
     * @param base the base [AaCcGgTt]
     * @return the transversion of the base, or the input base if it's not one of the understood ones
     */
    static public char transversion(char base) {
        switch (base) {
            case 'A':
            case 'a': return 'C';
            case 'C':
            case 'c': return 'A';
            case 'G':
            case 'g': return 'T';
            case 'T':
            case 't': return 'G';
            default: return base;
        }
    }

   /**
     * Return the complement (A <-> T or C <-> G) of a base, or the specified base if it can't be complemented (i.e. an ambiguous base).
     *
     * @param base the base [AaCcGgTt]
     * @return the complementary base, or the input base if it's not one of the understood ones
     */
    static public char simpleComplement(char base) {
        switch (base) {
            case 'A':
            case 'a': return 'T';
            case 'C':
            case 'c': return 'G';
            case 'G':
            case 'g': return 'C';
            case 'T':
            case 't': return 'A';
            default: return base;
        }
    }

    /**
     * Reverse complement a byte array of bases (that is, chars casted to bytes, *not* base indices in byte form)
     *
     * @param bases the byte array of bases
     * @return the reverse complement of the base byte array
     */
    static public byte[] simpleReverseComplement(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = (byte)simpleComplement((char) bases[bases.length - 1 - i]);
        }

        return rcbases;
    }

    /**
     * Complement a byte array of bases (that is, chars casted to bytes, *not* base indices in byte form)
     *
     * @param bases  the byte array of bases
     * @return the complement of the base byte array
     */
    static public byte[] simpleComplement(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = (byte)simpleComplement((char) bases[i]);
        }

        return rcbases;
    }

    /**
     * Reverse complement a char array of bases
     *
     * @param bases the char array of bases
     * @return the reverse complement of the char byte array
     */
    static public char[] simpleReverseComplement(char[] bases) {
        char[] rcbases = new char[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = simpleComplement(bases[bases.length - 1 - i]);
        }

        return rcbases;
    }

    /**
     * Complement a char array of bases
     *
     * @param bases  the char array of bases
     * @return the complement of the base char array
     */
    static public char[] simpleComplement(char[] bases) {
        char[] rcbases = new char[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = simpleComplement(bases[i]);
        }

        return rcbases;
    }

    /**
     * Reverse complement a String of bases.  Preserves ambiguous bases.
     *
     * @param bases  the String of bases
     * @return the reverse complement of the String
     */
    static public String simpleReverseComplement(String bases) {
        return new String(simpleReverseComplement(bases.getBytes()));
    }


    /**
     * Complement a String of bases.  Preserves ambiguous bases.
     *
     * @param bases  the String of bases
     * @return the complement of the String
     */
    static public String simpleComplement(String bases) {
        return new String(simpleComplement(bases.getBytes()));
    }

    /**
     * Reverse a byte array of bases
     * 
     * @param bases  the byte array of bases
     * @return the reverse of the base byte array
     */
    static public byte[] reverse(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - i - 1];
        }

        return rcbases;
    }

    /**
     * Reverse an int array of bases
     *
     * @param bases  the int array of bases
     * @return the reverse of the base int array
     */
    static public int[] reverse(int[] bases) {
        int[] rcbases = new int[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - i - 1];
        }

        return rcbases;
    }

    /**
     * Reverse (NOT reverse-complement!!) a string
     *
     * @param bases  input string
     * @return the reversed string
     */
    static public String reverse(String bases) {
        return new String( reverse( bases.getBytes() )) ;
    }

    /**
     * For the most frequent base in the sequence, return the percentage of the read it constitutes.
     *
     * @param sequence  the read sequence
     * @return  the percentage of the read that's made up of the most frequent base
     */
    static public double mostFrequentBaseFraction(byte[] sequence) {
        int[] baseCounts = new int[4];

        for ( byte base : sequence ) {
            int baseIndex = simpleBaseToBaseIndex((char) base);

            if (baseIndex >= 0) {
                baseCounts[baseIndex]++;
            }
        }

        int mostFrequentBaseIndex = 0;
        for (int baseIndex = 1; baseIndex < 4; baseIndex++) {
            if (baseCounts[baseIndex] > baseCounts[mostFrequentBaseIndex]) {
                mostFrequentBaseIndex = baseIndex;
            }
        }

        return ((double) baseCounts[mostFrequentBaseIndex])/((double) sequence.length);
    }

    /**
     * Return a random base index (A=0, C=1, G=2, T=3).
     *
     * @return a random base index (A=0, C=1, G=2, T=3)
     */
    static public int getRandomBaseIndex() {
        return getRandomBaseIndex(-1);
    }

    /**
     * Return a random base index, excluding some base index.
     *
     * @param excludeBaseIndex the base index to exclude
     * @return a random base index, excluding the one specified (A=0, C=1, G=2, T=3)
     */
    static public int getRandomBaseIndex(int excludeBaseIndex) {
        int randomBaseIndex = excludeBaseIndex;

        Random generator = new Random();
        while (randomBaseIndex == excludeBaseIndex) {
            randomBaseIndex = generator.nextInt(4);
        }

        return randomBaseIndex;
    }

    /**
     * Return a random base (A, C, G, T).
     *
     * @return a random base (A, C, G, T)
     */
    static public char getRandomBase() {
        return getRandomBase('.');
    }

    /**
     * Return a random base, excluding some base.
     *
     * @param excludeBase the base to exclude
     * @return a random base, excluding the one specified (A, C, G, T)
     */
    static public char getRandomBase(char excludeBase) {
        return BaseUtils.baseIndexToSimpleBase(getRandomBaseIndex(BaseUtils.simpleBaseToBaseIndex(excludeBase)));
    }
    
    
    /** Computes the smallest period >= minPeriod for the specified string. The period is defined as such p, 
     * that for all  i = 0... seq.length-1,  seq[ i % p ] = seq[i] (or equivalently seq[i] = seq[i+p] for i=0...seq.length-1-p).
     *  The sequence does <i>not</i> have to contain whole number of periods. For instance, "ACACACAC" has a period 
     *  of 2 (it has a period of 4 as well), and so does
     * "ACACA"; similarly, smallest periods of "CTCCTC", "CTCCT", and "CTCC" are all equal to 3. The "trivial" period is 
     * the length of the string itself, and it will always be returned if no smaller period can be found in the specified period range
     * or if specified minPeriod is greater than the sequence length.
     *   
     * @param seq
     * @return
     */
    public static int sequencePeriod(byte[] seq, int minPeriod) {
    	int period = ( minPeriod > seq.length ? seq.length : minPeriod );
    	// we assume that bases [0,period-1] repeat themselves and check this assumption
    	// until we find correct period
    	
    	for ( int pos = period ; pos < seq.length ; pos++ ) {
    		
    		int offset = pos % period; // we are currenlty 'offset' bases into the putative repeat of period 'period'
    		                                                // if our current hypothesis holds, base[pos] must be the same as base[offset]
    		
    		if ( Character.toUpperCase( seq[pos] ) !=
    				Character.toUpperCase( seq[offset] )
    			) {
    			
    			// period we have been trying so far does not work.
    			// two possibilities:
    			// A) offset = 0, i.e. current position pos must be start of the next repeat, but it is not;
    			//      in this case only bases from start up to the current one, inclusive, may form a repeat, if at all;
    		   //       so period is at least pos+1 (remember, pos is 0-based), then on the next loop re-entrance 
    			//      pos will be autoincremented and we will be checking next base
    			// B) offset != 0, i.e. the current base breaks the repeat, but maybe it starts a new one?
    			//     hence we should first check if it matches the first base of the sequence, and to do that
    			//     we set period to pos  (thus trying the hypothesis that bases from start up to the current one,
    			//     non-inclusive are repeated hereafter), and decrement pos (this will re-test current base against the first base
    			// on the next loop re-entrance after pos is autoincremented)
    			if ( offset == 0 ) period = pos+1;
    			else period = pos-- ;
    		
    		} 
    	}
    	return period;
    }

    public static byte[] charSeq2byteSeq(char[] seqIn) {
        byte[] seqOut = new byte[seqIn.length];
        for ( int i = 0; i < seqIn.length; i++ ) {
            seqOut[i] = (byte)seqIn[i];
        }
        return seqOut;
    }

}

/* code snippet for testing sequencePeriod():
 * 
 *     	String str = "CCTTG";
    	int p = 0;
    	System.out.print("Periods of " + str +" are:");
    	while ( p < str.length() ) {
    		p = sequencePeriod(str, p+1);
        	System.out.print(" "+p);
    	}
    	System.out.println(); System.exit(1);
*/

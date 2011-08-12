package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;


/**
 * BaseUtils contains some basic utilities for manipulating nucleotides.
 */
public class BaseUtils {
    public final static byte A = (byte)'A';
    public final static byte C = (byte)'C';
    public final static byte G = (byte)'G';
    public final static byte T = (byte)'T';

    public final static byte N = (byte)'N';
    public final static byte D = (byte)'D';

    //
    // todo -- we need a generalized base abstraction using the Base enum.
    //
    public final static byte[] BASES = { 'A', 'C', 'G', 'T' };
    public final static byte[] EXTENDED_BASES = { 'A', 'C', 'G', 'T', 'N', 'D' };

    public enum Base {
        A ( 'A', 0 ),
        C ( 'C', 1 ),
        G ( 'G', 2 ),
        T ( 'T', 3 );

        byte b;
        int index;
        private Base(char base, int index) {
            this.b = (byte)base;
            this.index = index;
        }

        public byte getBase() { return b; }
        public char getBaseAsChar() { return (char)b; }
        public int getIndex() { return index; }

        public boolean sameBase(byte o) { return b == o; }
        public boolean sameBase(char o) { return b == (byte)o; }
        public boolean sameBase(int i)  { return index == i; }
    }


    // todo -- fix me (enums?)
    public static final byte DELETION_INDEX = 4;
    public static final byte NO_CALL_INDEX = 5; // (this is 'N')

    public static int gIndex = BaseUtils.simpleBaseToBaseIndex((byte)'G');
    public static int cIndex = BaseUtils.simpleBaseToBaseIndex((byte)'C');
    public static int aIndex = BaseUtils.simpleBaseToBaseIndex((byte)'A');
    public static int tIndex = BaseUtils.simpleBaseToBaseIndex((byte)'T');


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

    public static boolean isTransition( byte base1, byte  base2 ) {
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
        return simpleBaseToBaseIndex(base1) == simpleBaseToBaseIndex(base2);
    }

    static public boolean extendedBasesAreEqual(byte base1, byte base2) {
        return extendedBaseToBaseIndex(base1) == extendedBaseToBaseIndex(base2);
    }


    /**
     * Converts a IUPAC nucleotide code to a pair of bases
     *
     * @param code
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    @Deprecated
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
    static public int simpleBaseToBaseIndex(byte base) {
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


    /**
     * Converts a simple base to a base index
     *
     * @param base  [AaCcGgTt]
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    @Deprecated
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

    static public int extendedBaseToBaseIndex(byte base) {
        switch (base) {
            case 'd':
            case 'D': return DELETION_INDEX;
            case 'n':
            case 'N': return NO_CALL_INDEX;

            default: return simpleBaseToBaseIndex(base);
        }
    }

    @Deprecated
    static public boolean isRegularBase(char base) {
        return simpleBaseToBaseIndex(base) != -1;
    }

    static public boolean isRegularBase(byte base) {
        return simpleBaseToBaseIndex(base) != -1;
    }

    @Deprecated
    static public boolean isNBase(char base) {
        return isNBase((byte)base);
    }

    static public boolean isNBase(byte base) {
        return base == 'N' || base == 'n';
    }

    /**
     * Converts a base index to a simple base
     *
     * @param baseIndex  0, 1, 2, 3
     * @return A, C, G, T, or '.' if the index can't be understood
     */
    static public byte baseIndexToSimpleBase(int baseIndex) {
        switch (baseIndex) {
            case 0: return 'A';
            case 1: return 'C';
            case 2: return 'G';
            case 3: return 'T';
            default: return '.';
        }
    }

    @Deprecated
    static public char baseIndexToSimpleBaseAsChar(int baseIndex) {
        return (char)baseIndexToSimpleBase(baseIndex);
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
    @Deprecated
    static public char crossTalkPartnerBase(char base) {
        return (char)baseIndexToSimpleBase(crossTalkPartnerIndex(simpleBaseToBaseIndex(base)));
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

   /**
     * Return the complement (A <-> T or C <-> G) of a base, or the specified base if it can't be complemented (i.e. an ambiguous base).
     *
     * @param base the base [AaCcGgTt]
     * @return the complementary base, or the input base if it's not one of the understood ones
     */
    static public byte simpleComplement(byte base) {
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

    @Deprecated
    static public char simpleComplement(char base) {
        return (char)simpleComplement((byte)base);
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
            rcbases[i] = simpleComplement(bases[bases.length - 1 - i]);
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
            rcbases[i] = simpleComplement(bases[i]);
        }

        return rcbases;
    }

    /**
     * Reverse complement a char array of bases
     *
     * @param bases the char array of bases
     * @return the reverse complement of the char byte array
     */
    @Deprecated
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
    @Deprecated
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
    @Deprecated
    static public String simpleReverseComplement(String bases) {
        return new String(simpleReverseComplement(bases.getBytes()));
    }


    /**
     * Complement a String of bases.  Preserves ambiguous bases.
     *
     * @param bases  the String of bases
     * @return the complement of the String
     */
    @Deprecated
    static public String simpleComplement(String bases) {
        return new String(simpleComplement(bases.getBytes()));
    }

    /**
     * Returns the index of the most common base in the basecounts array. To be used with
     * pileup.getBaseCounts.
     *
     * @param baseCounts counts of a,c,g,t in order.
     * @return the index of the most common base
     */
    static public int mostFrequentBaseIndex(int[] baseCounts) {
        int mostFrequentBaseIndex = 0;
        for (int baseIndex = 1; baseIndex < 4; baseIndex++) {
            if (baseCounts[baseIndex] > baseCounts[mostFrequentBaseIndex]) {
                mostFrequentBaseIndex = baseIndex;
            }
        }
        return mostFrequentBaseIndex;
    }

    static public int mostFrequentBaseIndexNotRef(int[] baseCounts, int refBaseIndex) {
        int tmp = baseCounts[refBaseIndex];
        baseCounts[refBaseIndex] = -1;
        int result = mostFrequentBaseIndex(baseCounts);
        baseCounts[refBaseIndex] = tmp;
        return result;
    }

    static public int mostFrequentBaseIndexNotRef(int[] baseCounts, byte refSimpleBase) {
        return mostFrequentBaseIndexNotRef(baseCounts, simpleBaseToBaseIndex(refSimpleBase));
    }

    /**
     * Returns the most common base in the basecounts array. To be used with pileup.getBaseCounts.
     *
     * @param  baseCounts counts of a,c,g,t in order.
     * @return the most common base
     */
    static public byte mostFrequentSimpleBase(int[] baseCounts) {
        return baseIndexToSimpleBase(mostFrequentBaseIndex(baseCounts));
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
            int baseIndex = simpleBaseToBaseIndex(base);

            if (baseIndex >= 0) {
                baseCounts[baseIndex]++;
            }
        }

        int mostFrequentBaseIndex = mostFrequentBaseIndex(baseCounts);

        return ((double) baseCounts[mostFrequentBaseIndex])/((double) sequence.length);
    }

    // --------------------------------------------------------------------------------
    //
    // random bases
    //
    // --------------------------------------------------------------------------------

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

        while (randomBaseIndex == excludeBaseIndex) {
            randomBaseIndex = GenomeAnalysisEngine.getRandomGenerator().nextInt(4);
        }

        return randomBaseIndex;
    }

    /**
     * Return a random base (A, C, G, T).
     *
     * @return a random base (A, C, G, T)
     */
    @Deprecated
    static public byte getRandomBase() {
        return getRandomBase('.');
    }

    /**
     * Return a random base, excluding some base.
     *
     * @param excludeBase the base to exclude
     * @return a random base, excluding the one specified (A, C, G, T)
     */
    @Deprecated
    static public byte getRandomBase(char excludeBase) {
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

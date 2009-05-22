package org.broadinstitute.sting.utils;

/**
 * BaseUtils contains some basic utilities for manipulating nucleotides.
 *
 * @author Kiran Garimella
 */
public class BaseUtils {

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
    public static BaseSubstitutionType SNPSubstitutionType( char base1, char base2 ) {
        BaseSubstitutionType t = isTransition(base1, base2) ? BaseSubstitutionType.TRANSITION : BaseSubstitutionType.TRANSVERSION;
        //System.out.printf("SNPSubstitutionType( char %c, char %c ) => %s%n", base1, base2, t);
        return t;
    }

    public static boolean isTransition( char base1, char base2 ) {
        int b1 = simpleBaseToBaseIndex(base1);
        int b2 = simpleBaseToBaseIndex(base2);
        return b1 == 0 && b2 == 2 || b1 == 2 && b2 == 0 ||
               b1 == 1 && b2 == 3 || b1 == 3 && b2 == 1;
    }

    public static boolean isTransversion( char base1, char base2 ) {
        return ! isTransition(base1, base2);
    }

    /** Private constructor.  No instantiating this class! */
    private BaseUtils() {}

    static public boolean basesAreEqual(byte base1, byte base2) {
        return simpleBaseToBaseIndex((char)base1) == simpleBaseToBaseIndex((char)base2);
    }


    /**
     * Converts a simple base to a base index
     *
     * @param base  [AaCcGgTt]
     * @return 0, 1, 2, 3, or -1 if the base can't be understood
     */
    static public int simpleBaseToBaseIndex(char base) {
        switch (base) {
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

    /**
     * Return the complement of a base.
     *
     * @param base  the base [AaCcGgTt]
     * @return the complementary base
     */
    static public byte simpleComplement(char base) {
        switch (base) {
            case 'A':
            case 'a': return 'T';
            case 'C':
            case 'c': return 'G';
            case 'G':
            case 'g': return 'C';
            case 'T':
            case 't': return 'A';
            default: return '.';
        }
    }

    /**
     * Reverse complement a byte array of bases (that is, chars casted to bytes, *not* base indices in byte form)
     *
     * @param bases  the byte array of bases
     * @return the reverse complement of the base byte array
     */
    static public byte[] simpleReverseComplement(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = simpleComplement((char) bases[bases.length - 1]);
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
        char[] rcbases = new char[bases.length()];

        for (int i = 0; i < bases.length(); i++) {
            char base = bases.charAt(bases.length() - 1);
            char rcbase = (base == 'N' || base == '.') ? base : (char) simpleComplement(base);

            rcbases[i] = rcbase;
        }

        return new String(rcbases);
    }

    /**
     * Reverse a byte array of bases
     * @param bases  the byte array of bases
     * @return the reverse of the base byte array
     */
    static public byte[] reverse(byte[] bases) {
        byte[] rcbases = new byte[bases.length];

        for (int i = 0; i < bases.length; i++) {
            rcbases[i] = bases[bases.length - 1];
        }

        return rcbases;
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
}

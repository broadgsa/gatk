package org.broadinstitute.sting.gatk.contexts.variantcontext;

import org.broadinstitute.sting.utils.BaseUtils;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.util.Collection;

/**
 * Immutable representation of an allele
 *
 * Types of alleles:
 *
 * Ref: a t C g a // C is the reference base
 *
 *    : a t G g a // C base is a G in some individuals
 *
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 *    : a t CAg a // A base is inserted w.r.t. the reference sequence
 *
 * In these cases, where are the alleles?
 *
 * SNP polymorphism of C/G  -> { C , G } -> C is the reference allele
 * 1 base deletion of C     -> { C , - } -> C is the reference allele
 * 1 base insertion of A    -> { - ; A } -> Null is the reference allele
 *
 * Suppose I see a the following in the population:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t G g a // C base is a G in some individuals
 *    : a t - g a // C base is deleted w.r.t. the reference
 *
 * How do I represent this?  There are three segregating alleles:
 *
 *  { C , G , - }
 *
 * Now suppose I have this more complex example:
 *
 * Ref: a t C g a // C is the reference base
 *    : a t - g a
 *    : a t - - a
 *    : a t CAg a
 *
 * There are actually four segregating alleles:
 *
 *   { C g , - g, - -, and CAg } over bases 2-4
 *
 * However, the molecular equivalence explicitly listed above is usually discarded, so the actual
 * segregating alleles are:
 *
 *   { C g, g, -, C a g }
 *
 * Critically, it should be possible to apply an allele to a reference sequence to create the
 * correct haplotype sequence:
 *
 * Allele + reference => haplotype
 *
 * For convenience, we are going to create Alleles where the GenomeLoc of the allele is stored outside of the
 * Allele object itself.  So there's an idea of an A/C polymorphism independent of it's surrounding context.
 *
 * Given list of alleles it's possible to determine the "type" of the variation
 *
 *      A / C @ loc => SNP with
 *      - / A => INDEL
 *
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion.
 *
 * Alelle also supports is concept of a NO_CALL allele.  This Allele represents a haplotype that couldn't be
 * determined. This is usually represented by a '.' allele.
 *
 * Note that Alleles store all bases as bytes, in **UPPER CASE**.  So 'atc' == 'ATC' from the perspective of an
 * Allele.

 * @author ebanks, depristo
 */
public class Allele implements Comparable<Allele> {
    private static final byte[] EMPTY_ALLELE_BASES = new byte[0];

    private boolean isRef = false;
    private boolean isNull = false;
    private boolean isNoCall = false;

    private byte[] bases = null;

    public final static String NULL_ALLELE_STRING = "-";
    public final static String NO_CALL_STRING = ".";
    /** A generic static NO_CALL allele for use */
    public final static Allele NO_CALL = new Allele(NO_CALL_STRING);

    /**
     * Create a new Allele that includes bases and if tagged as the reference allele if isRef == true.  If bases
     * == '-', a Null allele is created.  If bases ==  '.', a no call Allele is created.
     *
     * @param bases the DNA sequence of this variation, '-', of '.'
     * @param isRef should we make this a reference allele?
     * @throws IllegalArgumentException if bases contains illegal characters or is otherwise malformated
     */
    public Allele(byte[] bases, boolean isRef) {
        if ( bases == null )
            throw new IllegalArgumentException("Constructor: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");

        // standardize our representation of null allele and bases
        if ( wouldBeNullAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            isNull = true;
        } else if ( wouldBeNoCallAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            isNoCall = true;
            if ( isRef ) throw new IllegalArgumentException("Cannot tag a NoCall allele as the reference allele");
        }
        else
            bases = new String(bases).toUpperCase().getBytes(); // todo -- slow performance

        this.isRef = isRef;
        this.bases = bases;

        if ( ! acceptableAlleleBases(bases) )
            throw new IllegalArgumentException("Unexpected base in allele bases " + new String(bases));
    }

    /**
     * Do the bases represent the null allele?
     */
    public static boolean wouldBeNullAllele(byte[] bases) {
        return (bases.length == 1 && bases[0] == '-') || bases.length == 0;
    }

    /** Do the bases represent the NO_CALL allele? */
    public static boolean wouldBeNoCallAllele(byte[] bases) {
        return bases.length == 1 && bases[0] == '.';
    }

    /** Do the bases represent the null allele? */
    public static boolean acceptableAlleleBases(String bases) {
        return acceptableAlleleBases(bases.getBytes());
    }

    /** Can we create an allele from bases, including NO_CALL and Null alleles? */
    public static boolean acceptableAlleleBases(byte[] bases) {
        if ( wouldBeNullAllele(bases) || wouldBeNoCallAllele(bases) )
            return true;

        for ( byte b : bases ) {
            if ( ! BaseUtils.isRegularBase(b) ) {
                return false;
            }
        }

        return true;
    }

    /**
     * @see Allele(byte[], boolean)
     *
     * @param bases
     * @param isRef
     */
    public Allele(String bases, boolean isRef) {
        this(bases.getBytes(), isRef);
    }

    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases
     */
    public Allele(String bases) { this(bases, false); }

    /**
     * Creates a non-Ref allele.  @see Allele(byte[], boolean) for full information
     *
     * @param bases
     */
    public Allele(byte[] bases) { this(bases, false); }

    // ---------------------------------------------------------------------------------------------------------
    //
    // accessor routines
    //
    // ---------------------------------------------------------------------------------------------------------

    /** Returns true if this is the null allele */
    public boolean isNull()             { return isNull; }
    /** Returns true if this is not the null allele */
    public boolean isNonNull()          { return ! isNull(); }

    /** Returns true if this is the NO_CALL allele */
    public boolean isNoCall()           { return isNoCall; }
    /** Returns true if this is the not the NO_CALL allele */
    public boolean isCalled()           { return ! isNoCall(); }

    /** Returns true if this Allele is the reference allele */
    public boolean isReference()        { return isRef; }
    /** Returns true if this Allele is not the reference allele */
    public boolean isNonReference()     { return ! isReference(); }

    /** Returns a nice string representation of this object */
    public String toString() {
        return (isNull() ? "-" : ( isNoCall() ? "." : new String(getBases()))) + (isReference() ? "*" : "");
    }

    /**
     * Return the DNA bases segregating in this allele.  Note this isn't reference polarized,
     * so the Null allele is represented by a vector of length 0
     *
     * @return the segregating bases
     */
    public byte[] getBases() { return bases; }

    /**
     * @param other  the other allele
     *
     * @return true if these alleles are equal
     */
    public boolean equals(Allele other) {
        return equals(other, false);
    }

    /**
     * Returns true if this and other are equal.  If ignoreRefState is true, then doesn't require both alleles has the
     * same ref tag
     *
     * @param other
     * @param ignoreRefState
     * @return
     */
    public boolean equals(Allele other, boolean ignoreRefState) {
        return (isRef == other.isRef || ignoreRefState) && isNull == other.isNull && isNoCall == other.isNoCall && this.basesMatch(other.getBases());
    }

    /**
     * Returns true if this Alelle contains the same bases as test, regardless of its reference status.  Also handles
     * Null and NO_CALL alleles
     *
     * @param test
     * @return
     */
    public boolean basesMatch(byte[] test) { return bases == test || Arrays.equals(bases, test); }

    /**
     * Returns true if this Alelle contains the same bases as test, regardless of its reference status.  Also handles
     * Null and NO_CALL alleles
     *
     * @param test
     * @return
     */
    public boolean basesMatch(String test) { return basesMatch(test.toUpperCase().getBytes()); }

    /**
     * Returns true if this Alelle contains the same bases as test, regardless of its reference status.  Also handles
     * Null and NO_CALL alleles
     *
     * @param test
     * @return
     */
    public boolean basesMatch(Allele test) { return basesMatch(test.getBases()); }

    /**
     * Returns the length of this allele.  Null and NO_CALL alleles have 0 length.
     * @return
     */
    public int length() {
        return bases.length;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // useful static functions
    //
    // ---------------------------------------------------------------------------------------------------------

    public static Allele getMatchingAllele(Collection<Allele> allAlleles, String alleleBases) {
        return getMatchingAllele(allAlleles, alleleBases.getBytes());
    }

    public static Allele getMatchingAllele(Collection<Allele> allAlleles, byte[] alleleBases) {
        for ( Allele a : allAlleles ) {
            if ( a.basesMatch(alleleBases) ) {
                return a;
            }
        }

        if ( wouldBeNoCallAllele(alleleBases) )
            return NO_CALL;
        else
            return null;    // couldn't find anything
    }

    public static List<Allele> resolveAlleles(List<Allele> possibleAlleles, List<String> alleleStrings) {
        List<Allele> myAlleles = new ArrayList<Allele>(alleleStrings.size());

        for ( String alleleString : alleleStrings ) {
            Allele allele = getMatchingAllele(possibleAlleles, alleleString);

            if ( allele == null ) {
                if ( Allele.wouldBeNoCallAllele(alleleString.getBytes()) ) {
                    allele = new Allele(alleleString);
                } else {
                    throw new IllegalArgumentException("Allele " + alleleString + " not present in the list of alleles " + possibleAlleles);
                }
            }

            myAlleles.add(allele);
        }

        return myAlleles;
    }

    public int compareTo(Allele other) {
        if ( isReference() && other.isNonReference() )
            return -1;
        else if ( isNonReference() && other.isReference() ) 
            return 1;
        else
            return new String(getBases()).compareTo(new String(other.getBases())); // todo -- potential performance issue
    }
}

package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.utils.BaseUtils;

import java.util.Arrays;

/**
 * @author ebanks, depristo
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
 * 1 base insertion of A    -> { - ; A } -> NULL is the reference allele
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
 * If you know where allele is the reference, you can determine whether the variant is an insertion or deletion
 */
public class Allele {
    private static final byte[] EMPTY_ALLELE_BASES = new byte[0];
//    private static final byte[] NULL_ALLELE_BASES = new byte[0];
//    private static final byte[] NO_CALL_ALLELE_BASES = ".".getBytes();

    private boolean isRef = false;
    private boolean isNull = false;
    private boolean isNoCall = false;

    private byte[] bases = null;

    public final static Allele NO_CALL = new Allele(".");

    public Allele(byte[] bases, boolean isRef) {
        if ( bases == null )
            throw new IllegalArgumentException("Constructor: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");

        // standardize our representation of null allele and bases
        if ( wouldBeNullAllele(bases) ) {
            bases = EMPTY_ALLELE_BASES;
            isNull = true;
        }
        if ( wouldBeNoCallAllele(bases) ) {
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

    public final static boolean wouldBeNullAllele(byte[] bases) {
        return (bases.length == 1 && bases[0] == '-') || bases.length == 0;
    }

    public final static boolean wouldBeNoCallAllele(byte[] bases) {
        return bases.length == 1 && bases[0] == '.';
    }


    public final static boolean acceptableAlleleBases(String bases) {
        return acceptableAlleleBases(bases.getBytes());
    }
    
    public final static boolean acceptableAlleleBases(byte[] bases) {
        if ( (bases.length == 1 && bases[0] == '-') || bases.length == 0)
            return true;

        for ( byte b : bases ) {
            if ( ! BaseUtils.isRegularBase(b) ) {
                return false;
            }
        }

        return true;
    }

    /** null allele creation method */
    public Allele(boolean isRef) {
        this("", isRef);
    }

    public Allele(String bases, boolean isRef) {
        this(bases.getBytes(), isRef);
    }

    public Allele()             { this(false); }
    public Allele(String bases) { this(bases, false); }
    public Allele(byte[] bases) { this(bases, false); }

    //
    //
    // accessor routines
    //
    //
    public boolean isNull()             { return isNull; }
    public boolean isNonNull()          { return ! isNull(); }

    public boolean isNoCall()           { return isNoCall; }
    public boolean isCalled()           { return ! isNoCall(); }

    public boolean isReference()        { return isRef; }
    public boolean isNonReference()     { return ! isReference(); }

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
        return isRef == other.isRef && isNull == other.isNull && isNoCall == other.isNoCall && this.basesMatch(other.getBases());
    }

    // todo -- notice case insensitivity
    public boolean basesMatch(byte[] test) { return bases == test || Arrays.equals(bases, test); }
    public boolean basesMatch(String test) { return basesMatch(test.toUpperCase().getBytes()); }
    public boolean basesMatch(Allele test) { return basesMatch(test.getBases()); }

    public int length() {
        return bases.length;
    }
}

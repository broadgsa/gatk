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
    private boolean isRef = false;
    private byte[] bases = null;

    public Allele(byte[] bases, boolean isRef) {
        bases = new String(bases).toUpperCase().getBytes(); // todo -- slow performance
        this.isRef = isRef;

        if ( bases == null )
            throw new IllegalArgumentException("Constructor: the Allele base string cannot be null; use new Allele() or new Allele(\"\") to create a Null allele");

        this.bases = bases;
        for ( byte b : bases ) {
            if ( ! BaseUtils.isRegularBase(b) ) {
                throw new IllegalArgumentException("Unexpected base in allele bases " + new String(bases));
            }
        }
    }

    /** null allele creation method */
    public Allele(boolean isRef) {
        this("", isRef);
    }

    public Allele(String bases, boolean isRef) {
        this(bases.getBytes(), isRef);
    }

    //
    //
    // accessor routines
    //
    //
    public boolean isNullAllele()       { return length() == 0; }
    public boolean isNonNullAllele()    { return ! isNullAllele(); }

    public boolean isReference()        { return isRef; }
    public boolean isNonReference()     { return ! isReference(); }


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
        return Arrays.equals(bases, other.getBases());
    }

    public int length() {
        return bases.length;
    }
}

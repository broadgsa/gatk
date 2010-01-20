package org.broadinstitute.sting.gatk.refdata;


/**
 * @author ebanks
 *         <p/>
 *         Class Allele
 *         <p/>
 *         This class emcompasses all the basic information about an allele
 */
public class Allele {

    private AlleleType type;

    private String bases;


    // the types of variants we currently allow
    public enum AlleleType {
        REFERENCE, SNP, INSERTION, DELETION, INVERSION, UNKNOWN_POINT_ALLELE
    }

    public Allele(AlleleType type, String bases) {
        this.type = type;
        if ( bases == null )
            throw new IllegalArgumentException("Constructor: the Allele base string cannot be null");
        if ( type == AlleleType.DELETION && bases.length() > 0 )
            throw new IllegalArgumentException("Constructor: deletions cannot have observed bases");
        if ( (type == AlleleType.REFERENCE || type == AlleleType.SNP || type == AlleleType.UNKNOWN_POINT_ALLELE) && bases.length() > 1 )
            throw new IllegalArgumentException("Constructor: point alleles cannot have more than one observed base");
        this.bases = bases.toUpperCase();
    }

    /**
     * convenience method for switching over the allele type
     *
     * @return the AlleleType of this allele
     **/
    public AlleleType getType() { return type; }

    /**
     * convenience method for SNPs
     *
     * @return true if this is a SNP, false otherwise
     */
    public boolean isSNP() { return type == AlleleType.SNP; }

    /**
     * convenience method for variants
     *
     * @return true if this is a variant allele, false if it's reference
     */
    public boolean isVariant() { return type != AlleleType.REFERENCE; }


    /**
     * convenience method for indels
     *
     * @return true if this is an indel, false otherwise
     */
    public boolean isIndel() { return type == AlleleType.INSERTION || type == AlleleType.DELETION; }


    /**
     * For deletions, this method returns an empty String.
     * For everything else, observed bases for the allele are returned.
     *
     * @return the bases, as a string
     */
    public String getBases() { return bases; }

    /**
     * @param other  the other allele
     *
     * @return true if these alleles are equal
     */
    public boolean equals(Allele other) {
        return type == other.getType() && bases.equals(other.getBases());
    }

}
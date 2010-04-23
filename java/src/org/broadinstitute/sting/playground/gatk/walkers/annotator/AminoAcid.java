package org.broadinstitute.sting.playground.gatk.walkers.annotator;

/**
 * Represents a single amino acid.
 */
public class AminoAcid {
    private String name;
    private String code;
    private String letter;


    /**
     * Constructor.
     *
     * @param letter The 1 letter code. (eg. I). This is '*' for the stop codon.
     * @param name The full name of the AA (eg. Isoleucine).
     * @param code The 3 letter code. (eg. Ile).
     */
    public AminoAcid( String letter, String name, String code) {
        this.name = name;
        this.code = code;
        this.letter = letter;
    }

    /** Equality based on the amino acid code. */
    public boolean equals(Object o) {
        if (this == o) { return true; }
        if (o == null || !(o instanceof AminoAcid)) { return false; }

        final AminoAcid aminoAcid = (AminoAcid) o;
        return !(getCode() != null ? !getCode().equals(aminoAcid.getCode()) : aminoAcid.getCode() != null);
    }

    /** Hashes the three letter code. */
    public int hashCode() {
        return (getCode() != null ? getCode().hashCode() : 0);
    }

    /**
     * Returns the full name of this AA.
     */
    public String getName() { return name; }

    /**
     * Returns the 1 letter code for this AA.
     */
    public String getLetter() { return letter; }

    /**
     * Returns the 3 letter code.
     */
    public String getCode() { return code; }


    /** Returns true if the amino acid is really just a stop codon. */
    public boolean isStop() {
        return "*".equals(getLetter());
    }
}

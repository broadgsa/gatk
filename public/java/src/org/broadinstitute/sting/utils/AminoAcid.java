/*
 * Copyright (c) 2010 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils;

/**
 * Represents a single amino acid.
 */
public class AminoAcid {
    private String name;
    private String threeLetterCode;
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
        this.threeLetterCode = code;
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
    public String getCode() { return threeLetterCode; }


    /** Returns true if the amino acid is really just a stop codon. */
    public boolean isStop() {
        return "*".equals(getLetter());
    }

    /** Returns true if the amino acid is really just a stop codon. */
    public boolean isUnknown() {
        return "X".equals(getLetter());
    }

    public String toString() {
        return name;
    }
}

/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package org.broadinstitute.sting.gatk.refdata;

import edu.mit.broad.picard.variation.KnownVariant;
import edu.mit.broad.picard.variation.VariantType;
import edu.mit.broad.picard.variation.DbSnpFileReader;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Collections;
import java.util.Arrays;

/**
 * @author alecw@broadinstitute.org
 */
public class rodPicardDbSNP implements VariationRod {

    private final KnownVariant knownVariant;
    private final MyGenomeLoc loc;
    private final List<String> alleles;

    public rodPicardDbSNP(final KnownVariant knownVariant) {
        this.knownVariant = knownVariant;
        alleles = Collections.unmodifiableList(Arrays.asList(knownVariant.getObservedAlleles()));
        loc = new MyGenomeLoc(knownVariant);
    }

    /**
     * are we bi-allelic?
     */
    public boolean isBiallelic() {
        return alleles.size() == 2;
    }

    /**
     * get the frequency of this variant, if we're a variant.  If we're reference this method
     * should return 0.  If we can't provide an alternate allele frequency, this should also
     * return 0.
     * <p/>
     * WARNING: This method is only valid for biAllelic data, the contract is to check isBiallelic()
     * before calling this method
     *
     * @return double the minor allele frequency
     */
    public double getNonRefAlleleFrequency() {
        return knownVariant.getMinorAlleleFrequency();
    }

    /**
     * A convenience method, for switching over the variation type
     *
     * @return the VARIANT_TYPE of the current variant
     */
    public VARIANT_TYPE getType() {
        switch (knownVariant.getType()) {
            case SNP:
                return VARIANT_TYPE.SNP;
            case insertion:
            case deletion:
                return VARIANT_TYPE.INDEL;
        }
        return null;
    }

    /**
     * are we a SNP? If not we're a Indel/deletion or the reference.  This method must be call before you use
     * the convenience methods getAlternativeBaseForSNP or getReferenceForSNP, to ensure that you're working with a SNP
     *
     * @return true if we're a SNP
     */
    public boolean isSNP() {
        return knownVariant.getType() == VariantType.SNP;
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isInsertion() {
        return knownVariant.getType() == VariantType.insertion;
    }

    /**
     * are we an deletion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isDeletion() {
        return knownVariant.getType() == VariantType.deletion;
    }

    /**
     * are we a variant that represents the reference allele?
     *
     * @return false if we're a variant(indel, delete, SNP, etc), true if we're hom ref
     */
    public boolean isReference() {
        return false;  // snp locations are never "reference", there's always a variant
    }

    /**
     * are we an insertion or a deletion? yes, then return true.  No? false.
     *
     * @return true if we're an insertion or deletion
     */
    public boolean isIndel() {
        return getType() == VARIANT_TYPE.INDEL;
    }

    public String getName() {
        return "PicarddbSNP";
    }

    public boolean parseLine(final Object header, final String[] parts) throws IOException {
        throw new UnsupportedOperationException("Not a text format so this should never be called.");
    }

    public String toSimpleString() {
        return getName() + ":" + knownVariant.getObservedAllelesString();
    }

    public String repl() {
        throw new UnsupportedOperationException("Not a text format so this should never be called.");
    }

    /**
     * Used by the ROD system to determine how to split input lines
     *
     * @return Regex string delimiter separating fields
     */
    public String delimiterRegex() {
        throw new UnsupportedOperationException("Not a text format so this should never be called.");
    }

    /**
     * get the location of this Variant
     *
     * @return a GenomeLoc
     */
    public GenomeLoc getLocation() {
        return loc;
    }

    public int compareTo(final ReferenceOrderedDatum that) {
        return getLocation().compareTo(that.getLocation());
    }

    /**
     * Backdoor hook to read header, meta-data, etc. associated with the file.  Will be
     * called by the ROD system before streaming starts
     *
     * @param source source data file on disk from which this rod stream will be pulled
     * @return a header object that will be passed to parseLine command
     */
    public Object initialize(final File source) throws FileNotFoundException {
        throw new UnsupportedOperationException("Not a text format so this should never be called.");
    }

    /**
     * get the reference base(s) for this Variant
     *
     * @return the reference base or bases, as a string
     */
    public String getReference() {
        // Not known without looking in reference fasta.
        return null;
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the postive number space log based error estimate
     */
    public double getNegLog10PError() {
        // Copied from rodDbSNP
        return 4; // -log10(0.0001)
    }

    /**
     * gets the alternate alleles.  This method should return all the alleles present at the location,
     * NOT including the reference base.  This is returned as a string list with no guarantee ordering
     * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
     * frequency).
     *
     * @return an alternate allele list
     */
    public List<String> getAlternateAlleleList() {
        throw new UnsupportedOperationException("Don't know which is the reference allele");
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles. If the reference base is not an allele in this varation
     * it will not be in the list (i.e. there is no guarantee that the reference base is in the list).
     *
     * @return an alternate allele list
     */
    public List<String> getAlleleList() {
        return alleles;
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    public char getAlternativeBaseForSNP() {
        throw new UnsupportedOperationException("Don't know which is the reference allele");
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    public char getReferenceForSNP() {
        throw new UnsupportedOperationException("Don't know which is the reference allele");
    }

    public static Iterator<rodPicardDbSNP> createIterator(final String trackName, final File f) throws IOException {
        return new rodPicardDbSNPIterator(f);
    }

    private static class MyGenomeLoc extends GenomeLoc {
        private MyGenomeLoc(final KnownVariant knownVariant) {
            // GenomeLoc is one-based closed interval.
            super(knownVariant.getReferenceSequence(), knownVariant.getSequenceIndex(), knownVariant.getStartPos(),
                    knownVariant.getEndPos());
        }
    }

    private static class rodPicardDbSNPIterator implements Iterator<rodPicardDbSNP> {
        private final DbSnpFileReader reader;

        private rodPicardDbSNPIterator(final File f) {
            reader = new DbSnpFileReader(f);
        }

        public boolean hasNext() {
            return reader.hasNext();
        }

        public rodPicardDbSNP next() {
            return new rodPicardDbSNP(reader.next());
        }

        public void remove() {
            reader.remove();
        }

    }

}

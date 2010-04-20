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

package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.glf.GLFReader;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;
import org.broadinstitute.sting.utils.genotype.glf.GLFSingleCall;
import org.broadinstitute.sting.utils.genotype.glf.GLFVariableLengthCall;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;


/**
 * @author aaron
 *         <p/>
 *         Class RodGLF
 *         <p/>
 *         the rod class for GLF data.
 */
public class RodGLF implements Iterator<RodGLF>, ReferenceOrderedDatum {
    public GLFReader mReader;
    private final String mName;
    private GenomeLoc mLoc;
    public GLFRecord mRecord;

    public RodGLF(String name) {
        mName = name;
    }

    /**
     * get the name
     *
     * @return the name
     */
    public String getName() {
        return mName;
    }

    /**
     * Backdoor hook to read header, meta-data, etc. associated with the file.  Will be
     * called by the ROD system before streaming starts
     *
     * @param source source data file on disk from which this rod stream will be pulled
     *
     * @return a header object that will be passed to parseLine command
     */
    public Object initialize(File source) throws FileNotFoundException {
        mReader = new GLFReader(source);
        return null;
    }

    public String toSimpleString() {
        return toString();
    }

    /** @return a string representation of the ROD GLF object */
    public String toString() {
        return String.format("%s\t%d\t%s\t%d\t%d\t%4.4f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
                             mLoc.getContig(),
                             mLoc.getStart(),
                             mRecord.getRefBase(),
                             mRecord.getReadDepth(),
                             mRecord.getRmsMapQ(),
                             getBestGenotypeValue(1),
                             ((GLFSingleCall) mRecord).getLikelihoods()[0],
                             ((GLFSingleCall) mRecord).getLikelihoods()[1],
                             ((GLFSingleCall) mRecord).getLikelihoods()[2],
                             ((GLFSingleCall) mRecord).getLikelihoods()[3],
                             ((GLFSingleCall) mRecord).getLikelihoods()[4],
                             ((GLFSingleCall) mRecord).getLikelihoods()[5],
                             ((GLFSingleCall) mRecord).getLikelihoods()[6],
                             ((GLFSingleCall) mRecord).getLikelihoods()[7],
                             ((GLFSingleCall) mRecord).getLikelihoods()[8],
                             ((GLFSingleCall) mRecord).getLikelihoods()[9]


        );
    }

    public String repl() {
        return this.toString();
    }

    /**
     * Used by the ROD system to determine how to split input lines
     *
     * @return Regex string delimiter separating fields
     */
    public String delimiterRegex() {
        return "";
    }

    /**
     * return a genome loc representing the current location
     *
     * @return the geonome loc
     */
    public GenomeLoc getLocation() {
        return mLoc;
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    public String getReference() {
        return mRecord.getRefBase().toString();
    }

    /** are we bi-allelic? */
    public boolean isBiallelic() {
        return true;
    }

    /**
     * Returns true if all observed alleles are reference alleles. All is<Variant> methods (where Variant=SNP,Insertion, etc) should
     * return false at such site to ensure consistency. This method is included for use with genotyping calls (isGenotype()==true), it makes
     * no sense for, e.g. dbSNP and should return false for the latter.
     *
     * @return
     */
    public boolean isReference() {
        return (!isSNP());
    }

    /**
     * are we an insertion or a deletion? yes, then return true.  No? Well, false it is.
     *
     * @return true if we're an insertion or deletion
     */
    public boolean isIndel() {
        return (isDeletion() || isInsertion());
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of
     *
     * @return a char, representing the alternate base
     */
    public char getAlternativeBaseForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        List<String> alleles = this.getAlternateAlleleList();
        if (alleles.size() != 1) throw new StingException("We're not biAllelic()");
        return Utils.stringToChar(alleles.get(0));
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    public char getReferenceForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        return Utils.stringToChar(getReference());
    }

    /**
     * Is this variant a SNP?
     *
     * @return true or false
     */
    public boolean isSNP() {
        return ((mRecord.getRecordType() == GLFRecord.RECORD_TYPE.SINGLE) &&
                (!getBestGenotype(1).toString().equals(refString(mRecord.getRefBase().toChar()))));
    }

    /**
     * return a string representing the reference
     *
     * @param ref the reference character
     *
     * @return a string for the homozygous ref in a diploid
     */
    private static String refString(char ref) {
        return new String(new char[]{ref, ref});
    }

    /**
     * Get the nth best genotype (one based), i.e. to get the best genotype pass in 1,
     * the second best 2, etdc.
     *
     * @param nthBest the nth best genotype to get (1 based, NOT ZERO BASED)
     *
     * @return a GENOTYPE object representing the nth best genotype
     */
    public LikelihoodObject.GENOTYPE getBestGenotype(int nthBest) {
        Integer[] sorted = MathUtils.sortPermutation(((GLFSingleCall) mRecord).getLikelihoods());
        return LikelihoodObject.GENOTYPE.values()[sorted[nthBest - 1]];
    }

    /**
     * Get the nth best genotype value (one based), i.e. to get the best genotype pass in 1,
     * the second best 2, etdc.
     *
     * @param nthBest the nth best genotype value to get
     *
     * @return a GENOTYPE object representing the nth best genotype
     */
    public double getBestGenotypeValue(int nthBest) {
        Integer[] sorted = MathUtils.sortPermutation(((GLFSingleCall) mRecord).getLikelihoods());
        return (((GLFSingleCall) mRecord).getLikelihoods())[sorted[nthBest - 1]];
    }

    /**
     * Is this variant an insertion? The contract requires isIndel() to return true
     * if this method returns true.
     *
     * @return true or false
     */
    public boolean isInsertion() {
        return ((mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE) &&
                ((GLFVariableLengthCall) mRecord).getIndelLen1() > 0);
    }

    /**
     * Is this variant a deletion? The contract requires isIndel() to return true
     * if isDeletion() returns true.
     *
     * @return true or false
     */
    public boolean isDeletion() {
        return ((mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE) &&
                ((GLFVariableLengthCall) mRecord).getIndelLen1() < 0);
    }

    /**
     * Returns minor allele frequency.
     *
     * @return
     */
    public double getNonRefAlleleFrequency() {
        return 0;
    }

    /**
     * Returns phred-mapped confidence in variation event (e.g. MAQ's SNP confidence, or AlleleCaller's best vs. ref).
     *
     * @return
     */
    public double getNegLog10PError() {
        String ref = new String() + mRecord.getRefBase() + mRecord.getRefBase();
        int index = 0;
        for (LikelihoodObject.GENOTYPE g : LikelihoodObject.GENOTYPE.values()) {
            if (g.toString().equals(ref)) break;
            index++;
        }
        return Math.abs(getBestGenotypeValue(1) - ((GLFSingleCall) mRecord).getLikelihoods()[index]) / GLFRecord.LIKELIHOOD_SCALE_FACTOR;
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
        LikelihoodObject.GENOTYPE genotype = getBestGenotype(1);
        List<String> ret = new ArrayList<String>();
        for (char c : genotype.toString().toCharArray()) {
            if (!String.valueOf(c).equals(this.getReference())) ret.add(String.valueOf(c));
        }
        return ret;
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    public List<String> getAlleleList() {
        LikelihoodObject.GENOTYPE genotype = getBestGenotype(1);
        List<String> list = new ArrayList<String>();
        if (genotype.toString().contains(this.getReference())) list.add(this.getReference());
        for (char c : genotype.toString().toCharArray())
            if (c != Utils.stringToChar(getReference()))
                list.add(String.valueOf(c));
        return list;
    }

    public int length() {
        return 1;
    }

    public int compareTo(ReferenceOrderedDatum that) {
        return this.mLoc.compareTo(that.getLocation());
    }

    /**
     * the parse line, which is not used by the GLF rod
     *
     * @param header the header to pass in
     * @param parts  the string object
     *
     * @return false, alwayss
     * @throws java.io.IOException
     */
    public boolean parseLine(Object header, String[] parts) throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean hasNext() {
        return (mReader.hasNext());
    }

    /**
     * @return the next element in the iteration.
     * @throws NoSuchElementException - iterator has no more elements.
     */
    @Override
    public RodGLF next() {
        if (!this.hasNext()) throw new NoSuchElementException("RodGLF next called on iterator with no more elements");
        mRecord = mReader.next();
        mLoc = GenomeLocParser.createGenomeLoc(mRecord.getContig(), mRecord.getPosition(), mRecord.getPosition());
        return this;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("GLF Rods don't support the remove() function");
    }

    public static RodGLF createIterator(String name, File file) {
        RodGLF glf = new RodGLF(name);
        try {
            glf.initialize(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find file " + file);
        }
        return glf;
    }

}


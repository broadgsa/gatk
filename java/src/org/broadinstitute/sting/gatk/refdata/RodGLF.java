package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.glf.GLFReader;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;
import org.broadinstitute.sting.utils.genotype.glf.SinglePointCall;
import org.broadinstitute.sting.utils.genotype.glf.VariableLengthCall;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class RodGLF
 *         <p/>
 *         the rod class for GLF data.
 */
public class RodGLF implements ReferenceOrderedDatum, AllelicVariant, Iterator<RodGLF> {
    static int count = 0;
    public GLFReader mReader;
    private final String mName;
    private GenomeLoc mLoc;
    public GLFRecord mRecord;

    public RodGLF(String name) {
        mName = name;
    }

    /**
     * get the name
     * @return the name
     */
    @Override
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
    @Override
    public Object initialize(File source) throws FileNotFoundException {
        mReader = new GLFReader(source);
        return null;
    }

    @Override
    public String toSimpleString() {
        return toString();
    }

    /**
     * @return a string representation of the ROD GLF object
     */
    public String toString() {
        return String.format("%s\t%d\t%s\t%d\t%d\t%4.4f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
                mLoc.getContig(),
                mLoc.getStart(),
                mRecord.getRefBase(),
                mRecord.getReadDepth(),
                mRecord.getRmsMapQ(),
                getBestGenotypeValue(1),
                ((SinglePointCall)mRecord).getLikelihoods()[0],
                ((SinglePointCall)mRecord).getLikelihoods()[1],
                ((SinglePointCall)mRecord).getLikelihoods()[2],
                ((SinglePointCall)mRecord).getLikelihoods()[3],
                ((SinglePointCall)mRecord).getLikelihoods()[4],
                ((SinglePointCall)mRecord).getLikelihoods()[5],
                ((SinglePointCall)mRecord).getLikelihoods()[6],
                ((SinglePointCall)mRecord).getLikelihoods()[7],
                ((SinglePointCall)mRecord).getLikelihoods()[8],
                ((SinglePointCall)mRecord).getLikelihoods()[9]


        );
    }

    @Override
    public String repl() {
        return this.toString();
    }

    /**
     * Used by the ROD system to determine how to split input lines
     *
     * @return Regex string delimiter separating fields
     */
    @Override
    public String delimiterRegex() {
        return "";
    }

    /**
     * return a genome loc representing the current location
     * @return the geonome loc
     */
    @Override
    public GenomeLoc getLocation() {
        return mLoc;
    }

    /**
     * Returns bases in the reference allele as a String. String can be empty (as in insertion into
     * the reference), can contain a single character (as in SNP or one-base deletion), or multiple characters
     * (for longer indels).
     *
     * @return reference allele, forward strand
     */
    @Override
    public String getRefBasesFWD() {
        if (mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE) return "";
        return String.valueOf(mRecord.getRefBase());
    }

    /**
     * Returns reference (major) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * This doesn't make much sense to me, what if the best genotype is hom non-ref?
     *
     * @return reference base on the forward strand
     */
    @Override
    public char getRefSnpFWD() throws IllegalStateException {
        if (!isSNP()) throw new IllegalStateException("Current GLF Record is not a SNP");
        return mRecord.getRefBase().toChar();
    }

    /**
     * Returns bases in the alternative allele as a String. String can be empty (as in deletion from
     * the reference), can contain a single character (as in SNP or one-base insertion), or multiple characters
     * (for longer indels).
     *
     * @return alternative allele, forward strand
     */
    @Override
    public String getAltBasesFWD() {
        return getBestGenotype(2).toString();
    }

    /**
     * Returns alternative (minor) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return alternative allele base on the forward starnd
     */
    @Override
    public char getAltSnpFWD() throws IllegalStateException {
        if (!isSNP()) {
            throw new IllegalStateException("Not a SNP");
        }
        String str = getBestGenotype(1).toString();
        if (String.valueOf(str.charAt(0)).equals(mRecord.getRefBase().toString())) {
            return str.charAt(1);
        }
        return str.charAt(0);
    }

    /**
     * Returns true if all observed alleles are reference alleles. All is<Variant> methods (where Variant=SNP,Insertion, etc) should
     * return false at such site to ensure consistency. This method is included for use with genotyping calls (isGenotype()==true), it makes
     * no sense for, e.g. dbSNP and should return false for the latter.
     *
     * @return
     */
    @Override
    public boolean isReference() {
        return (!isSNP());
    }

    /**
     * Is this variant a SNP?
     *
     * @return true or false
     */
    @Override
    public boolean isSNP() {
        return ((mRecord.getRecordType() == GLFRecord.RECORD_TYPE.SINGLE) &&
                (!getBestGenotype(1).toString().equals(refString(mRecord.getRefBase().toChar()))));
    }

    /**
     * return a string representing the reference
     * @param ref the reference character
     * @return a string for the homozygous ref in a diploid
     */
    private static String refString(char ref) {
        return new String(new char[]{ref, ref});
    }

    /**
     * Get the nth best genotype (one based), i.e. to get the best genotype pass in 1,
     * the second best 2, etdc.
     *
     * @param nthBest the nth best genotype to get
     *
     * @return a GENOTYPE object representing the nth best genotype
     */
    public LikelihoodObject.GENOTYPE getBestGenotype(int nthBest) {
        Integer[] sorted = Utils.SortPermutation(((SinglePointCall) mRecord).getLikelihoods());
        return LikelihoodObject.GENOTYPE.values()[sorted[nthBest-1]];
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
        Integer[] sorted = Utils.SortPermutation(((SinglePointCall) mRecord).getLikelihoods());
        return (((SinglePointCall) mRecord).getLikelihoods())[sorted[nthBest-1]];
    }

    /**
     * Is this variant an insertion? The contract requires isIndel() to return true
     * if this method returns true.
     *
     * @return true or false
     */
    @Override
    public boolean isInsertion() {
        return ((mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE) &&
                ((VariableLengthCall) mRecord).getIndelLen1() > 0);
    }

    /**
     * Is this variant a deletion? The contract requires isIndel() to return true
     * if isDeletion() returns true.
     *
     * @return true or false
     */
    @Override
    public boolean isDeletion() {
        return ((mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE) &&
                ((VariableLengthCall) mRecord).getIndelLen1() < 0);
    }

    /**
     * Is this variant an insertion or a deletion? The contract requires
     * this to be true if either isInsertion() or isDeletion() returns true. However,
     * this method is currently allowed to return true even if neither of isInsertion()
     * and isDeletion() does.
     *
     * @return
     */
    @Override
    public boolean isIndel() {
        return (mRecord.getRecordType() == GLFRecord.RECORD_TYPE.VARIABLE);
    }

    /**
     * Returns minor allele frequency.
     *
     * @return
     */
    @Override
    public double getMAF() {
        return 0;
    }

    /**
     * Returns heterozygosity, a more accessible general feature of a variant.
     *
     * @return
     */
    @Override
    public double getHeterozygosity() {
        return 0;
    }

    /**
     * Is this variant an actual genotype (such as individual call from sequencing, HapMap chip etc), or
     * population allelic variant (call from pooled sequencing, dbSNP site etc). Only if variant is a genotype, there
     * is a meaningful question of, e.g., whether it is a het, or homogeneous non-ref.
     *
     * @return true if this variant is an actual genotype.
     */
    @Override
    public boolean isGenotype() {
        return true;
    }

    /**
     * Returns phred-mapped confidence in variation event (e.g. MAQ's SNP confidence, or AlleleCaller's best vs. ref).
     *
     * @return
     */
    @Override
    public double getVariationConfidence() {
        String ref = new String() + mRecord.getRefBase() + mRecord.getRefBase();
        int index = 0;
        for (LikelihoodObject.GENOTYPE g: LikelihoodObject.GENOTYPE.values()) {
            if (g.toString().equals(ref)) break;
            index++;
        }
        return Math.abs(getBestGenotypeValue(1) - ((SinglePointCall)mRecord).getLikelihoods()[index]) / GLFRecord.LIKELIHOOD_SCALE_FACTOR;
    }

    /**
     * Returns phred-mapped confidence in called genotype (e.g. MAQ's consensus confidence, or AlleleCaller's
     * best vs next-best.
     *
     * @return
     */
    @Override
    public double getConsensusConfidence() {
        return Math.abs(getBestGenotypeValue(1) - getBestGenotypeValue(2)) / GLFRecord.LIKELIHOOD_SCALE_FACTOR;
    }

    /**
     * Returns actual observed genotype. Allowed to return more than two alleles (@see #getPloidy()). If this variant
     * is not a genotype, should throw an IllegalStateException.
     *
     * @return
     */
    @Override
    public List<String> getGenotype() throws IllegalStateException {
        List<String> ret = new ArrayList<String>();
        ret.add(this.getBestGenotype(1).toString());
        return ret;
    }

    /**
     * Return actual number of observed alleles (chromosomes) in the genotype. If this variant is not a genotype,
     * should throw IllegalStateException.
     *
     * @return
     */
    @Override
    public int getPloidy() throws IllegalStateException {
        return 2; // we're so haploid it hurts
    }

    /**
     * Returns true if the site has at most two known or observed alleles (ref and non-ref), or false if there are > 2 allelic variants known or observed. When
     * the implementing class is a genotype, alleles should be always counted including the reference allele whether it was observed in the particular
     * individual or not: i.e. if the reference is 'C', then both 'CA' and 'AA' genotypes must be reported as bi-allelic, while 'AT' is <i>not</i> bi-allelic (since there are
     * two allelic variants, 'A' and 'T' <i>in addition</i> to the (known) reference variant 'C').
     *
     * @return
     */
    @Override
    public boolean isBiallelic() {
        return false;
    }

    @Override
    public int compareTo(ReferenceOrderedDatum that) {
        return this.mLoc.compareTo(that.getLocation());
    }

    /**
     * the parse line, which is not used by the GLF rod
     * @param header the header to pass in
     * @param parts the string object
     * @return false, alwayss
     * @throws java.io.IOException
     */
    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        return false;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public boolean hasNext() {
        return (mReader.hasNext());
    }

    @Override
    public RodGLF next() {
        mRecord = mReader.next();
        mLoc = GenomeLocParser.createGenomeLoc(mReader.getReferenceName(), mReader.getCurrentLocation(), mReader.getCurrentLocation());
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


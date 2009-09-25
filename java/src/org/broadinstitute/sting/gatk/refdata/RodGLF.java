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
public class RodGLF implements VariationRod, Iterator<RodGLF> {
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
     *
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

    /** @return a string representation of the ROD GLF object */
    public String toString() {
        return String.format("%s\t%d\t%s\t%d\t%d\t%4.4f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",
                             mLoc.getContig(),
                             mLoc.getStart(),
                             mRecord.getRefBase(),
                             mRecord.getReadDepth(),
                             mRecord.getRmsMapQ(),
                             getBestGenotypeValue(1),
                             ((SinglePointCall) mRecord).getLikelihoods()[0],
                             ((SinglePointCall) mRecord).getLikelihoods()[1],
                             ((SinglePointCall) mRecord).getLikelihoods()[2],
                             ((SinglePointCall) mRecord).getLikelihoods()[3],
                             ((SinglePointCall) mRecord).getLikelihoods()[4],
                             ((SinglePointCall) mRecord).getLikelihoods()[5],
                             ((SinglePointCall) mRecord).getLikelihoods()[6],
                             ((SinglePointCall) mRecord).getLikelihoods()[7],
                             ((SinglePointCall) mRecord).getLikelihoods()[8],
                             ((SinglePointCall) mRecord).getLikelihoods()[9]


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
     *
     * @return the geonome loc
     */
    @Override
    public GenomeLoc getLocation() {
        return mLoc;
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
        return mRecord.getRefBase().toString();
    }

    /** are we bi-allelic? */
    @Override
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
    @Override
    public boolean isReference() {
        return (!isSNP());
    }

    /**
     * are we an insertion or a deletion? yes, then return true.  No? Well, false it is.
     *
     * @return true if we're an insertion or deletion
     */
    @Override
    public boolean isIndel() {
        return (isDeletion() || isInsertion());
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException in the case
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        if (getAlternateBases().charAt(0) == this.getReference().charAt(0))
            return getAlternateBases().charAt(1);
        return getAlternateBases().charAt(0);

    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        if (!this.isSNP()) throw new IllegalStateException("we're not a SNP");
        if (getAlternateBases().charAt(0) == this.getReference().charAt(0))
            return getAlternateBases().charAt(0);
        return getAlternateBases().charAt(1);

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
     * @param nthBest the nth best genotype to get
     *
     * @return a GENOTYPE object representing the nth best genotype
     */
    public LikelihoodObject.GENOTYPE getBestGenotype(int nthBest) {
        Integer[] sorted = Utils.SortPermutation(((SinglePointCall) mRecord).getLikelihoods());
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
        Integer[] sorted = Utils.SortPermutation(((SinglePointCall) mRecord).getLikelihoods());
        return (((SinglePointCall) mRecord).getLikelihoods())[sorted[nthBest - 1]];
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
     * get the base representation of this Variant
     *
     * @return a string, of ploidy
     */
    @Override
    public String getAlternateBases() {
        return this.getBestGenotype(0).toString();
    }

    /**
     * gets the alternate bases.  Use this method if teh allele count is greater then 2
     *
     * @return
     */
    @Override
    public List<String> getAlternateBaseList() {
        List<String> list = new ArrayList<String>();
        list.add(this.getAlternateBases());
        return list; 
    }

    /**
     * Returns minor allele frequency.
     *
     * @return
     */
    @Override
    public double getNonRefAlleleFrequency() {
        return 0;
    }

    /** @return the VARIANT_TYPE of the current variant */
    @Override
    public VARIANT_TYPE getType() {
        if (this.isSNP()) return VARIANT_TYPE.SNP;
        else if (this.isInsertion() || this.isDeletion()) return VARIANT_TYPE.INDEL;
        else return VARIANT_TYPE.REFERENCE;
    }


    /**
     * Returns phred-mapped confidence in variation event (e.g. MAQ's SNP confidence, or AlleleCaller's best vs. ref).
     *
     * @return
     */
    @Override
    public double getNegLog10PError() {
        String ref = new String() + mRecord.getRefBase() + mRecord.getRefBase();
        int index = 0;
        for (LikelihoodObject.GENOTYPE g : LikelihoodObject.GENOTYPE.values()) {
            if (g.toString().equals(ref)) break;
            index++;
        }
        return Math.abs(getBestGenotypeValue(1) - ((SinglePointCall) mRecord).getLikelihoods()[index]) / GLFRecord.LIKELIHOOD_SCALE_FACTOR;
    }

    public int length() {
        return 1;
    }

    @Override
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


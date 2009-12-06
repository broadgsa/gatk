package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class RodVCF
 *         <p/>
 *         An implementation of the ROD for VCF.
 */
public class RodVCF extends BasicReferenceOrderedDatum implements VariationRod, VariantBackedByGenotype, Iterator<RodVCF> {
    // our VCF related information
    private VCFReader mReader;
    public VCFRecord mCurrentRecord;

    public RodVCF(String name) {
        super(name);
    }

    private RodVCF(String name, VCFRecord currentRecord, VCFReader reader) {
        super(name);
        mCurrentRecord = currentRecord;
        mReader = reader;
    }

    public void assertNotNull() {
        if ( mCurrentRecord == null ) {
            throw new UnsupportedOperationException("The current Record is null");
        }
    }

    public void assertBiAllelic() {
        if ( !isBiallelic() )
            throw new StingException("This VCF rod is not bi-allelic.");
    }

    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        throw new UnsupportedOperationException("RodVCF does not support the parseLine method");
    }

    public Object initialize(final File source) throws FileNotFoundException {
        if ( mReader == null )
            mReader = new VCFReader(source);
        return mReader.getHeader();
    }

    @Override
    public String toString() {
        return (mCurrentRecord != null ? mCurrentRecord.toStringEncoding(mReader.getHeader()) : "");
    }

    public static RodVCF createIterator(String name, File file) {
        RodVCF vcf = new RodVCF(name);
        try {
            vcf.initialize(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find file " + file);
        }
        return vcf;
    }

    public boolean hasNonRefAlleleFrequency() {
        assertNotNull();
        return mCurrentRecord.getNonRefAlleleFrequency() > 0.0;
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    public double getNonRefAlleleFrequency() {
        assertNotNull();
        return mCurrentRecord.getNonRefAlleleFrequency();
    }

    public boolean hasStrandBias() {
        assertNotNull();
        return this.mCurrentRecord.getInfoValues().containsKey("SB");
    }

    /**
     * get the strand bias of this variant
     *
     * @return StrandBias with the stored slod
     */
    public double getStrandBias() {
        return hasStrandBias() ? Double.valueOf(this.mCurrentRecord.getInfoValues().get("SB")) : 0.0;
    }

    /** @return the VARIANT_TYPE of the current variant */
    public VARIANT_TYPE getType() {
        assertNotNull();
        return mCurrentRecord.getType();
    }

    /**
     * are we a SNP? If not we're a Indel/deletion
     *
     * @return true if we're a SNP
     */
    public boolean isSNP() {
        assertNotNull();
        assertBiAllelic();
        return mCurrentRecord.isSNP();
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isInsertion() {
        assertNotNull();
        assertBiAllelic();
        return mCurrentRecord.isInsertion();
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isDeletion() {
        assertNotNull();
        assertBiAllelic();
        return mCurrentRecord.isDeletion();
    }

    @Override
    public GenomeLoc getLocation() {
        assertNotNull();
        return mCurrentRecord.getLocation();
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    public String getReference() {
        assertNotNull();
        return mCurrentRecord.getReference();
    }

    /** are we bi-allelic? */
    public boolean isBiallelic() {
        assertNotNull();
        return mCurrentRecord.isBiallelic();
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    public double getNegLog10PError() {
        assertNotNull();
        return mCurrentRecord.getNegLog10PError();
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
        assertNotNull();
        return mCurrentRecord.getAlternateAlleleList();
    }

    /**
     * gets the alleles.  This method should return all the alleles present at the location,
     * including the reference base.  The first allele should always be the reference allele, followed
     * by an unordered list of alternate alleles.
     *
     * @return an alternate allele list
     */
    public List<String> getAlleleList() {
        assertNotNull();
        return mCurrentRecord.getAlleleList();
    }

    /**
     * are we truely a variant, given a reference
     *
     * @return false if we're a variant(indel, delete, SNP, etc), true if we're not
     */
    public boolean isReference() {
        assertNotNull();
        return mCurrentRecord.isReference();
    }

    /**
     * are we an insertion or a deletion? yes, then return true.  No? Well, false then.
     *
     * @return true if we're an insertion or deletion
     */
    public boolean isIndel() {
        assertNotNull();
        return mCurrentRecord.isIndel();
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    public char getAlternativeBaseForSNP() {
        assertNotNull();
        return mCurrentRecord.getAlternativeBaseForSNP();
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    public char getReferenceForSNP() {
        assertNotNull();
        return mCurrentRecord.getReferenceForSNP();
    }

    /**
     * get the genotype
     *
     * @return a map in lexigraphical order of the genotypes
     */
    public Genotype getCalledGenotype() {
        assertNotNull();
        return mCurrentRecord.getCalledGenotype();
    }

    /**
     * get the genotypes
     *
     * @return a list of the genotypes
     */
    public List<Genotype> getGenotypes() {
        assertNotNull();
        return mCurrentRecord.getGenotypes();
    }

    /**
     * get the genotypes
     *
     * @return a list of the genotypes
     */
    public List<VCFGenotypeRecord> getVCFGenotypeRecords() {
        assertNotNull();
        return mCurrentRecord.getVCFGenotypeRecords();
    }

    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    public boolean hasGenotype(DiploidGenotype x) {
        assertNotNull();
        return mCurrentRecord.hasGenotype(x);
    }

    public VCFHeader getHeader() {
        return mReader.getHeader();
    }

    public boolean hasNext() {
        return mReader.hasNext();
    }

    public RodVCF next() {
        mCurrentRecord = mReader.next();
        return new RodVCF(name, mCurrentRecord, mReader);
    }

    public void remove() {
        throw new UnsupportedOperationException("The remove operation is not supported for a VCF rod");
    }
}

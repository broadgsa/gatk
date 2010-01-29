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
import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class RodVCF
 *         <p/>
 *         An implementation of the ROD for VCF.
 */
public class RodVCF extends BasicReferenceOrderedDatum implements VariationRod, VariantBackedByGenotype, Iterator<RodVCF> {
    public VCFReader getReader() {
        return mReader;
    }

    // our VCF related information
    private VCFReader mReader;

    public VCFRecord getRecord() {
        return mCurrentRecord;
    }

    public VCFRecord mCurrentRecord;

    public RodVCF(String name) {
        super(name);
    }

    public RodVCF(String name, VCFRecord currentRecord, VCFReader reader) {
        super(name);
        mCurrentRecord = currentRecord;
        mReader = reader;
    }

    public void assertNotNull() {
        if ( mCurrentRecord == null ) {
            throw new UnsupportedOperationException("The current Record is null");
        }
    }

    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        throw new UnsupportedOperationException("RodVCF does not support the parseLine method");
    }

    public Object initialize(final File source) throws FileNotFoundException {
        if ( mReader == null ) {
            mReader = new VCFReader(source);
        }
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
        return this.mCurrentRecord.getInfoValues().containsKey(VCFRecord.STRAND_BIAS_KEY);
    }

    /**
     * get the strand bias of this variant
     *
     * @return StrandBias with the stored slod
     */
    public double getStrandBias() {
        return hasStrandBias() ? Double.valueOf(this.mCurrentRecord.getInfoValues().get(VCFRecord.STRAND_BIAS_KEY)) : 0.0;
    }

    /** @return the VARIANT_TYPE of the current variant */
    public VARIANT_TYPE getType() {
        assertNotNull();
        return mCurrentRecord.getType();
    }

    public String getID() {
        assertNotNull();
        return mCurrentRecord.getID();
    }

    /**
     * are we a SNP? If not we're a Indel/deletion
     *
     * @return true if we're a SNP
     */
    public boolean isSNP() {
        assertNotNull();
        return mCurrentRecord.isSNP();
    }

    /**
     * are we a novel site? Is there a DBSNP identifier
     * or a hapmap entry for the site?
     */

    public boolean isNovel() {
        assertNotNull();
        return mCurrentRecord.isNovel();
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isInsertion() {
        assertNotNull();
        return mCurrentRecord.isInsertion();
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    public boolean isDeletion() {
        assertNotNull();
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

    public double getQual() {
        assertNotNull();
        return mCurrentRecord.getQual();
    }

    public boolean hasAlternateAllele() {
        assertNotNull();
        return mCurrentRecord.hasAlternateAllele();
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

    public boolean hasGenotypeData() {
        assertNotNull();
        return mCurrentRecord.hasGenotypeData();
    }

    /**
     * get the genotype
     *
     * // todo -- WTF is this?  This is a deeply unsafe call
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
     * Returns the genotype associated with sample, or null if the genotype is missing
     *
     * @param sampleName the name of the sample genotype to fetch
     * @return
     */
    public Genotype getGenotype(final String sampleName) {
        return mCurrentRecord.getGenotype(sampleName);
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

    public String[] getSampleNames() {
        assertNotNull();
        return mCurrentRecord.getSampleNames();
    }

//    public Map<String, Genotype> getSampleGenotypes() {
//        String[] samples = getSampleNames();
//        List<Genotype> genotypes = getGenotypes();
//        HashMap<String, Genotype> map = new HashMap<String, Genotype>();
//
//        for ( int i = 0; i < samples.length; i++ ) {
//            map.put(samples[i], genotypes.get(i));
//        }
//
//        return map;
//    }

    public Map<String, String> getInfoValues() {
        assertNotNull();
        return mCurrentRecord.getInfoValues();
    }

    public String[] getFilteringCodes() {
        assertNotNull();
        return mCurrentRecord.getFilteringCodes();
    }

    public boolean isFiltered() {
        assertNotNull();
        return mCurrentRecord.isFiltered();
    }

//    public boolean hasFilteringCodes() {
//        assertNotNull();
//        return mCurrentRecord.hasFilteringCodes();
//    }

    public String getFilterString() {
        assertNotNull();
        return mCurrentRecord.getFilterString();
    }

    public VCFHeader getHeader() {
        return mReader.getHeader();
    }

    public boolean hasNext() {
        return mReader.hasNext();
    }

    /**
     * @return the next element in the iteration.
     * @throws NoSuchElementException - iterator has no more elements.
     */
    public RodVCF next() {
        if (!this.hasNext()) throw new NoSuchElementException("RodVCF next called on iterator with no more elements");

        // get the next record
        VCFRecord rec = mReader.next();

        // make sure the next VCF record isn't before the current record (we'll accept at the same location, the
        // spec doesn't indicate, and it seems like a valid use case)
        GenomeLoc curPosition = null;
        if (mCurrentRecord != null) curPosition = mCurrentRecord.getLocation();
        if (curPosition != null && rec != null && curPosition.compareTo(rec.getLocation()) > 0)
            throw new StingException("The next VCF record appears to be before the current (current location => " + curPosition.toString() +
                                     ", the next record position => " + rec.getLocation().toString() + " with line : " + rec.toStringEncoding(mReader.getHeader()) + "). " +
                                     "Check to make sure the input VCF file is correctly sorted.");

        // save off the previous record.  This is needed given how iterators are used in the ROD system;
        // we need to save off the last record
        mCurrentRecord = rec;
        return new RodVCF(name, rec, mReader);
    }

    public void remove() {
        throw new UnsupportedOperationException("The remove operation is not supported for a VCF rod");
    }
}

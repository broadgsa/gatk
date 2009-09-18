package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;


/**
 * @author aaron
 *         <p/>
 *         Class RodVCF
 *         <p/>
 *         An implementation of the ROD for VCF.
 */
public class RodVCF extends BasicReferenceOrderedDatum implements Variation, VariantBackedByGenotype, Iterator<RodVCF> {
    // our VCF related information
    private VCFReader mReader;
    public VCFRecord mCurrentRecord;

    public RodVCF(String name) {
        super(name);
    }

    public void assertNotNull() {
        if (mCurrentRecord == null) {
            throw new UnsupportedOperationException("The current Record is null");
        }
    }

    @Override
    public boolean parseLine(Object header, String[] parts) throws IOException {
        throw new UnsupportedOperationException("We don't support the parse line");
    }

    public Object initialize(final File source) throws FileNotFoundException {
        if (mReader == null) {
            mReader = new VCFReader(source);
        }
        return mReader.getHeader();
    }

    @Override
    public String toString() {
        if (this.mCurrentRecord != null)
            return this.mCurrentRecord.toString();
        else
            return "";
    }

    public Iterator<RodVCF> createIterator(String name, File file) {
        RodVCF vcf = new RodVCF(name);
        try {
            vcf.initialize(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find file " + file);
        }
        return vcf;
    }

    /**
     * get the frequency of this variant
     *
     * @return VariantFrequency with the stored frequency
     */
    @Override
    public double getNonRefAlleleFrequency() {
        assertNotNull();
        if (this.mCurrentRecord.getInfoValues().containsKey("AF")) {
            return Double.valueOf(this.mCurrentRecord.getInfoValues().get("AF"));
        } else {
            // this is the poor man's AF
            if (this.mCurrentRecord.getInfoValues().containsKey("AC") && this.mCurrentRecord.getInfoValues().containsKey("AN")) {
                String splt[] = mCurrentRecord.getInfoValues().get("AC").split(",");
                if (splt.length > 0) {
                    return (Double.valueOf(splt[0]) / Double.valueOf(mCurrentRecord.getInfoValues().get("AN")));
                }
            }
        }

        return 0.0;
    }

    /** @return the VARIANT_TYPE of the current variant */
    @Override
    public VARIANT_TYPE getType() {
        if (this.isSNP()) return VARIANT_TYPE.SNP;
        else if (this.isIndel()) return VARIANT_TYPE.INDEL;
        return VARIANT_TYPE.REFERENCE;
    }

    /**
     * are we a SNP? If not we're a Indel/deletion
     *
     * @return true if we're a SNP
     */
    @Override
    public boolean isSNP() {
        this.assertNotNull();
        if (!mCurrentRecord.hasAlternateAllele())
            return false;
        for (String alt : this.mCurrentRecord.getAlternateAlleles()) {
            if (alt.length() != 1)
                return false;
        }
        return true;
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    @Override
    public boolean isInsertion() {
         this.assertNotNull();
        if (!mCurrentRecord.hasAlternateAllele())
            return false;
        for (String alt : this.mCurrentRecord.getAlternateAlleles()) {
            if (alt.startsWith("I"))
                return true;
        }
        return false;
    }

    /**
     * are we an insertion?
     *
     * @return true if we are, false otherwise
     */
    @Override
    public boolean isDeletion() {
         this.assertNotNull();
        if (!mCurrentRecord.hasAlternateAllele())
            return false;
        for (String alt : this.mCurrentRecord.getAlternateAlleles()) {
            if (alt.startsWith("D"))
                return true;
        }
        return false;
    }

    @Override
    public GenomeLoc getLocation() {
        this.assertNotNull();
        return GenomeLocParser.createGenomeLoc(mCurrentRecord.getChromosome(), mCurrentRecord.getPosition());
    }

    /**
     * get the reference base(s) at this position
     *
     * @return the reference base or bases, as a string
     */
    @Override
    public String getReference() {
        return String.valueOf(mCurrentRecord.getReferenceBase());
    }

    /** are we bi-allelic? */
    @Override
    public boolean isBiallelic() {
        return (this.getAlternateBases().size() == 1);
    }

    /**
     * get the -1 * (log 10 of the error value)
     *
     * @return the log based error estimate
     */
    @Override
    public double getNegLog10PError() {
        // we're -10  log(error), we have to divide by 10
        return mCurrentRecord.getQual() / 10.0;
    }

    /**
     * are we truely a variant, given a reference
     *
     * @return false if we're a variant(indel, delete, SNP, etc), true if we're not
     */
    @Override
    public boolean isReference() {
        return (!mCurrentRecord.hasAlternateAllele());
    }

    /**
     * gets the alternate bases.  If this is homref, throws an UnsupportedOperationException
     *
     * @return
     */
    @Override
    public String getAlternateBase() {
        if (!this.isBiallelic()) throw new UnsupportedOperationException("We're not biallelic, so please call getAlternateBases instead");
        return this.mCurrentRecord.getAlternateAlleles().get(0);
    }

    /**
     * gets the alternate bases.  If this is homref, throws an UnsupportedOperationException
     *
     * @return
     */
    @Override
    public List<String> getAlternateBases() {
        return this.mCurrentRecord.getAlternateAlleles();
    }

    /**
     * are we an insertion or a deletion? yes, then return true.  No? Well, false then.
     *
     * @return true if we're an insertion or deletion
     */
    @Override
    public boolean isIndel() {
        return (!isSNP());
    }

    /**
     * gets the alternate base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     * of
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getAlternativeBaseForSNP() {
        if (!isSNP()) throw new IllegalStateException("we're not a SNP");
        return mCurrentRecord.getAlternateAlleles().get(0).charAt(0);
    }

    /**
     * gets the reference base is the case of a SNP.  Throws an IllegalStateException if we're not a SNP
     *
     * @return a char, representing the alternate base
     */
    @Override
    public char getReferenceForSNP() {
        if (!isSNP()) throw new IllegalStateException("we're not a SNP");
        return mCurrentRecord.getReferenceBase();
    }

    /**
     * get the genotype
     *
     * @return a map in lexigraphical order of the genotypes
     */
    @Override
    public Genotype getCalledGenotype() {
        throw new UnsupportedOperationException("We don't support this right now"); 
    }

    /**
     * get the genotypes
     *
     * @return a list of the genotypes
     */
    @Override
    public List<Genotype> getGenotypes() {
        List<Genotype> genotypes = new ArrayList<Genotype>();
        if (!this.mCurrentRecord.hasGenotypeData()) {
            return genotypes;
        }
        double refQual = (this.getNegLog10PError());
        // add the reference
        for (VCFGenotypeRecord rec : mCurrentRecord.getVCFGenotypeRecords()) {
            double qual = 0;
            if (rec.getAlleles().equals(this.getReference()))
                qual = refQual;
            else if (rec.getFields().containsKey("GQ"))
                qual = Double.valueOf(rec.getFields().get("GQ")) / 10.0;
            genotypes.add(new BasicGenotype(this.getLocation(), Utils.join("", rec.getAlleles()), this.getReference().charAt(0), qual));
        }
        return genotypes;
    }

    /**
     * a private helper method
     *
     * @return an array in lexigraphical order of the likelihoods
     */
    private Genotype getGenotype(DiploidGenotype x) {
        if (x.toString().equals(getReference()))
            return new BasicGenotype(this.getLocation(), getReference(), this.getReference().charAt(0), 0);
        for (VCFGenotypeRecord record : mCurrentRecord.getVCFGenotypeRecords()) {
            if (Utils.join("", record.getAlleles()).equals(x.toString())) {
                double qual = 0.0;
                if (record.getAlleles().equals(this.getReference()))
                    qual = this.getNegLog10PError();
                else if (record.getFields().containsKey("GQ"))
                    qual = Double.valueOf(record.getFields().get("GQ")) / 10.0;
                return new BasicGenotype(this.getLocation(), Utils.join("", record.getAlleles()), this.getReference().charAt(0), qual);
            }
        }
        return null;
    }

    /**
     * do we have the specified genotype?  not all backedByGenotypes
     * have all the genotype data.
     *
     * @param x the genotype
     *
     * @return true if available, false otherwise
     */
    @Override
    public boolean hasGenotype(DiploidGenotype x) {
        if (getGenotype(x) != null)
            return true;
        return false;
    }

    @Override
    public boolean hasNext() {
        return (mReader.hasNext());
    }

    @Override
    public RodVCF next() {
        mCurrentRecord = mReader.next();
        return this;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException("You cannot remove from a VCF rod");
    }
}

package org.broadinstitute.sting.oneoffprojects.refdata;

import org.broadinstitute.sting.gatk.refdata.BasicReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.VariationRod;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.List;

/*
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 29, 2010
 */
public class HapmapVCFROD extends BasicReferenceOrderedDatum implements VariationRod, VariantBackedByGenotype, Iterator<HapmapVCFROD> {
// This is a (hopefully temporary) wrapper class for certain VCF files that we want to protect from
// utilities that grab genotypes or sample names across all VCF files

    private RodVCF rod;

    public VCFReader getReader() {
        return rod.getReader();
    }

    public VCFRecord getRecord() {
        return rod.getRecord();
    }

    public HapmapVCFROD(String name) {
        super(name);
    }

    public HapmapVCFROD(String name, VCFRecord currentRecord, VCFReader reader) {
        super(name);
        rod = new RodVCF(name,currentRecord,reader);
    }

    public HapmapVCFROD(String name, RodVCF rod) {
        super(name);
        this.rod = rod;
    }

    public Object initialize(final File source) throws FileNotFoundException {
        rod = new RodVCF(name);
        rod.initialize(source);
        return rod.getHeader();
    }

    public boolean parseLine(Object obj, String[] args) {
        try {
            return rod.parseLine(obj,args);
        } catch (Exception e) {
            throw new UnsupportedOperationException("Parse line not supported",e);
        }
    }

    public double getNegLog10PError() {
        return rod.getNegLog10PError();
    }

    public List<Genotype> getGenotypes() {
        return null;
        //return rod.getGenotypes();
    }

    public String getReference() {
        return rod.getReference();
    }

    public String toString() {
        return rod.toString();
    }

    public List<String> getAlternateAlleleList() {
        return rod.getAlternateAlleleList();
    }

    public boolean isDeletion() {
        return rod.isDeletion();
    }

    public GenomeLoc getLocation() {
        return rod.getLocation();
    }

    public boolean isBiallelic() {
        return rod.isBiallelic();
    }

    public boolean isIndel() {
        return rod.isIndel();
    }

    public Variation.VARIANT_TYPE getType() {
        return rod.getType();
    }

    public boolean isSNP() {
        return rod.isSNP();
    }

    public boolean isReference() {
        return rod.isReference();
    }

    public double getNonRefAlleleFrequency() {
        return rod.getNonRefAlleleFrequency();
    }

    public char getAlternativeBaseForSNP() {
        return rod.getAlternativeBaseForSNP();
    }

    public boolean isInsertion() {
        return rod.isInsertion();
    }

    public List<String> getAlleleList() {
        return rod.getAlleleList();
    }

    public Genotype getCalledGenotype() {
        return null;
        //return rod.getCalledGenotype();
    }

    public char getReferenceForSNP() {
        return rod.getReferenceForSNP();
    }

    public boolean hasGenotype(DiploidGenotype g) {
        return false;
        //return rod.hasGenotype(g);
    }

    public VCFHeader getHeader() {
        return rod.getHeader();
    }

    public boolean hasNext() {
        return rod.hasNext();
    }

    public HapmapVCFROD next() {
        return new HapmapVCFROD(name,rod.next());
    }

    public void remove() {
        rod.remove();
    }

    public static HapmapVCFROD createIterator(String name, File file) {
        RodVCF vcf = new RodVCF(name);
        try {
            vcf.initialize(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to find file " + file);
        }
        return new HapmapVCFROD(name,vcf);
    }

}

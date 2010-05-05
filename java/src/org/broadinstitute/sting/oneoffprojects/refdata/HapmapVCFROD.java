package org.broadinstitute.sting.oneoffprojects.refdata;

import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.refdata.BasicReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;
import java.util.List;

/*
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Jan 29, 2010
 */
public class HapmapVCFROD extends BasicReferenceOrderedDatum implements Iterator<HapmapVCFROD> {
// This is a (hopefully temporary) wrapper class for certain VCF files that we want to protect from
// utilities that grab genotypes or sample names across all VCF files

    private VCFRecord record;
    private VCFReader reader;

    public VCFReader getReader() {
        return reader;
    }

    public VCFRecord getRecord() {
        return record;
    }

    public HapmapVCFROD(String name) {
        super(name);
    }

    public HapmapVCFROD(String name, VCFRecord currentRecord, VCFReader reader) {
        super(name);
        this.record = currentRecord;
        this.reader = reader;
    }

    public Object initialize(final File source) throws FileNotFoundException {
        reader = new VCFReader(source);
        if (reader.hasNext()) record = reader.next();
        return reader.getHeader();
    }

    public boolean parseLine(Object obj, String[] args) {
        return false;
    }

    public double getNegLog10PError() {
        return record.getNegLog10PError();
    }
    public String getReference() {
        return record.getReference();
    }

    public String toString() {
        return record.toString();
    }

    public List<String> getAlternateAlleleList() {
        return record.getAlternateAlleleList();
    }

    public boolean isDeletion() {
        return record.isDeletion();
    }

    public GenomeLoc getLocation() {
        return GenomeLocParser.createGenomeLoc(record.getChr(),record.getStart());
    }

    public boolean isBiallelic() {
        return record.isBiallelic();
    }

    public boolean isIndel() {
        return record.isIndel();
    }

    public boolean isSNP() {
        return record.isSNP();
    }

    public boolean isReference() {
        return record.isReference();
    }

    public double getNonRefAlleleFrequency() {
        return record.getNonRefAlleleFrequency();
    }

    public char getAlternativeBaseForSNP() {
        return record.getAlternativeBaseForSNP();
    }

    public boolean isInsertion() {
        return record.isInsertion();
    }

    public List<String> getAlleleList() {
        return record.getAlleleList();
    }

    public char getReferenceForSNP() {
        return record.getReferenceForSNP();
    }

    public boolean hasGenotype(DiploidGenotype g) {
        return false;
        //return rod.hasGenotype(g);
    }

    public VCFHeader getHeader() {
        return record.getHeader();
    }

    public boolean hasNext() {
        return reader.hasNext();
    }

    public HapmapVCFROD next() {
        return new HapmapVCFROD(name, record, reader);
    }

    public void remove() {
        throw new UnsupportedOperationException("Unable to remove");
    }

    public static HapmapVCFROD createIterator(String name, File file) {
        HapmapVCFROD vcf = new HapmapVCFROD(name);
        try {
            vcf.initialize(file);
        } catch (FileNotFoundException e) {
            throw new StingException("Unable to load file",e);
        }
        return vcf;
    }

}

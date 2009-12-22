package org.broadinstitute.sting.gatk.io.stubs;

import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;
import java.util.Set;

/**
 * Stub providing a passthrough for VCF files.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFGenotypeWriterStub extends GenotypeWriterStub<VCFGenotypeWriter> implements VCFGenotypeWriter {
    /**
     * Construct a new stub with the given engine and target file.
     * @param engine The engine, for extracting command-line arguments, etc.
     * @param genotypeFile Target file into which to write genotyping data.
     */
    public VCFGenotypeWriterStub(GenomeAnalysisEngine engine, File genotypeFile) {
        super(engine,genotypeFile);
    }

    /**
     * Gets the format of this stub.  We may want to discontinue use of this method and rely on instanceof comparisons.
     * @return VCF always.
     */
    public GenotypeWriterFactory.GENOTYPE_FORMAT getFormat() {
        return GenotypeWriterFactory.GENOTYPE_FORMAT.VCF;
    }

    /**
     * initialize this VCF header
     *
     * @param sampleNames  the sample names
     * @param headerInfo  the optional header fields
     */
    public void writeHeader(Set<String> sampleNames, Set<VCFHeaderLine> headerInfo) {
        outputTracker.getStorage(this).writeHeader(sampleNames,headerInfo);
    }

    /**
     * Add a given VCF record to the given output.
     * @param vcfRecord Record to add.
     */
    public void addRecord(VCFRecord vcfRecord) {
        outputTracker.getStorage(this).addRecord(vcfRecord);
    }
}

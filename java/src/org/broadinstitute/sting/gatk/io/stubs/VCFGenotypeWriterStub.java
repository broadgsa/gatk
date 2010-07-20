package org.broadinstitute.sting.gatk.io.stubs;

import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeWriter;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;
import java.io.PrintStream;

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
     * Construct a new stub with the given engine and target stream.
     * @param engine The engine, for extracting command-line arguments, etc.
     * @param genotypeStream Target stream into which to write genotyping data.
     */
    public VCFGenotypeWriterStub(GenomeAnalysisEngine engine, PrintStream genotypeStream) {
        super(engine,genotypeStream);
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
     * @param header  the header
     */
    public void writeHeader(VCFHeader header) {
        outputTracker.getStorage(this).writeHeader(header);
    }

    /**
     * Add a given VCF file to the writer.
     * @param file  file from which to add records
     */
    public void append(File file) {
        outputTracker.getStorage(this).append(file);
    }
}

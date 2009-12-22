package org.broadinstitute.sting.gatk.io.stubs;

import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.glf.GLFGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;

/**
 * Stub providing a passthrough for GLF files.
 *
 * @author mhanna
 * @version 0.1
 */
public class GLFGenotypeWriterStub extends GenotypeWriterStub<GLFGenotypeWriter> implements GLFGenotypeWriter {
    /**
     * Construct a new stub with the given engine and target file.
     * @param engine The engine, for extracting command-line arguments, etc.
     * @param genotypeFile Target file into which to write genotyping data.
     */
    public GLFGenotypeWriterStub(GenomeAnalysisEngine engine, File genotypeFile) {
        super(engine,genotypeFile);
    }

    /**
     * Gets the format of this stub.  We may want to discontinue use of this method and rely on instanceof comparisons.
     * @return GLF always.
     */
    public GenotypeWriterFactory.GENOTYPE_FORMAT getFormat() {
        return GenotypeWriterFactory.GENOTYPE_FORMAT.GLF;
    }

    /**
     * Write the GLF header to the target file.
     * @param headerText The header to write.
     */
    public void writeHeader(String headerText) {
        outputTracker.getStorage(this).writeHeader(headerText);
    }

    /**
     * add a GLF record to the output file
     *
     * @param contigName   the contig name
     * @param contigLength the contig length
     * @param rec          the GLF record to write.
     */
    public void addGLFRecord(String contigName, int contigLength, GLFRecord rec) {
        outputTracker.getStorage(this).addGLFRecord(contigName, contigLength, rec);
    }
    
}

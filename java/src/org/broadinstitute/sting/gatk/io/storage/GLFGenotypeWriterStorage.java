package org.broadinstitute.sting.gatk.io.storage;

import org.broadinstitute.sting.utils.genotype.glf.GLFReader;
import org.broadinstitute.sting.utils.genotype.glf.GLFRecord;
import org.broadinstitute.sting.utils.genotype.glf.GLFGenotypeWriter;
import org.broadinstitute.sting.gatk.io.stubs.GenotypeWriterStub;

import java.io.File;

/**
 * Provides temporary and permanent storage for genotypes in GLF format.
 *
 * @author mhanna
 * @version 0.1
 */
public class GLFGenotypeWriterStorage extends GenotypeWriterStorage<GLFGenotypeWriter> implements GLFGenotypeWriter {
    /**
     * Creates new (permanent) storage for GLF genotype writers.
     * @param stub Stub containing appropriate input parameters.
     */
    public GLFGenotypeWriterStorage(GenotypeWriterStub stub) {
        super(stub);
    }

    /**
     * Creates new (temporary) storage for GLF genotype writers.
     * @param stub Stub containing appropriate input parameters.
     * @param target Target file for output data.
     */
    public GLFGenotypeWriterStorage(GenotypeWriterStub stub, File target) {
        super(stub,target);
    }

    /**
     * Write the geli header to the target file.
     * @param headerText The header to write.
     */
    public void writeHeader(String headerText) {
        ((GLFGenotypeWriter)writer).writeHeader(headerText);
    }

    /**
     * add a GLF record to the output file
     *
     * @param contigName   the contig name
     * @param contigLength the contig length
     * @param rec          the GLF record to write.
     */
    public void addGLFRecord(String contigName, int contigLength, GLFRecord rec) {
        ((GLFGenotypeWriter)writer).addGLFRecord(contigName,contigLength,rec);
    }

    /**
     * Merges the stream backing up this temporary storage into the target.
     * @param target Target stream for the temporary storage.  May not be null.
     */
    @Override
    public void mergeInto(GLFGenotypeWriter target) {
        GLFReader reader = new GLFReader(file);
        while ( reader.hasNext() ) {
             GLFRecord rec = reader.next();
            target.addGLFRecord(rec.getContig(),(int)rec.getPosition(),rec);
        }
        reader.close();

        file.delete();        
    }    
}

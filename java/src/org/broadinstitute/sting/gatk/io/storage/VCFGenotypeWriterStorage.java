package org.broadinstitute.sting.gatk.io.storage;

import org.broad.tribble.vcf.VCFHeader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeWriter;
import org.broadinstitute.sting.gatk.io.stubs.GenotypeWriterStub;

import java.io.File;

/**
 * Provides temporary and permanent storage for genotypes in VCF format.
 *
 * @author mhanna
 * @version 0.1
 */
public class VCFGenotypeWriterStorage extends GenotypeWriterStorage<VCFGenotypeWriter> implements VCFGenotypeWriter {
    /**
     * Creates new (permanent) storage for VCF genotype writers.
     * @param stub Stub containing appropriate input parameters.
     */
    public VCFGenotypeWriterStorage(GenotypeWriterStub stub) {
        super(stub);
    }

    /**
     * Creates new (temporary) storage for VCF genotype writers.
     * @param stub Stub containing appropriate input parameters.
     * @param target Target file for output data.
     */
    public VCFGenotypeWriterStorage(GenotypeWriterStub stub,File target) {
        super(stub,target);
    }

    /**
     * initialize this VCF header
     *
     * @param header  the header
     */
    public void writeHeader(VCFHeader header) {
        ((VCFGenotypeWriter)writer).writeHeader(header);
    }

    /**
     * Add a given VCF file to the writer.
     * @param file  file from which to add records
     */
    public void append(File file) {
        ((VCFGenotypeWriter)writer).append(file);
    }

    /**
     * Merges the stream backing up this temporary storage into the target.
     * @param target Target stream for the temporary storage.  May not be null.
     */
    public void mergeInto(VCFGenotypeWriter target) {
        target.append(file);
        file.delete();
    }    
}

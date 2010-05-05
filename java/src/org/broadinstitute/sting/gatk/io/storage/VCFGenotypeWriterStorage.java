package org.broadinstitute.sting.gatk.io.storage;

import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeWriter;
import org.broadinstitute.sting.gatk.io.stubs.GenotypeWriterStub;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;

import java.io.File;
import java.util.Set;

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
     * @param sampleNames  the sample names
     * @param headerInfo  the optional header fields
     */
    public void writeHeader(Set<String> sampleNames, Set<VCFHeaderLine> headerInfo) {
        ((VCFGenotypeWriter)writer).writeHeader(sampleNames,headerInfo);    
    }

    /**
     * Add a given VCF record to the given output.
     * @param vcfRecord Record to add.
     */
    public void addRecord(VCFRecord vcfRecord) {
        ((VCFGenotypeWriter)writer).addRecord(vcfRecord);    
    }

    /**
     * set the validation stringency
     * @param value   validation stringency value
     */
    public void setValidationStringency(VALIDATION_STRINGENCY value) {
        ((VCFGenotypeWriter)writer).setValidationStringency(value);
    }

    /**
     * Merges the stream backing up this temporary storage into the target.
     * @param target Target stream for the temporary storage.  May not be null.
     */
    public void mergeInto(VCFGenotypeWriter target) {
        VCFReader reader = new VCFReader(file);
        while ( reader.hasNext() )
            target.addRecord(reader.next());
        reader.close();

        file.delete();        
    }    
}

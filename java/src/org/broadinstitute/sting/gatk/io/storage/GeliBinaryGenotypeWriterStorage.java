package org.broadinstitute.sting.gatk.io.storage;

import org.broadinstitute.sting.utils.genotype.geli.GeliGenotypeWriter;
import org.broadinstitute.sting.gatk.io.stubs.GenotypeWriterStub;

import java.io.File;

import edu.mit.broad.picard.genotype.geli.GeliFileReader;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileHeader;

/**
 * Provides temporary and permanent storage for genotypes in Geli binary format.
 *
 * @author mhanna
 * @version 0.1
 */
public class GeliBinaryGenotypeWriterStorage extends GenotypeWriterStorage<GeliGenotypeWriter> implements GeliGenotypeWriter {
    /**
     * Creates new (permanent) storage for geli binary genotype writers.
     * @param stub Stub containing appropriate input parameters.
     */
    public GeliBinaryGenotypeWriterStorage(GenotypeWriterStub stub) {
        super(stub);
    }

    /**
     * Creates new (temporary) storage for geli binary genotype writers.
     * @param stub Stub containing appropriate input parameters.
     * @param target Target file for output data.
     */
    public GeliBinaryGenotypeWriterStorage(GenotypeWriterStub stub, File target) {
        super(stub,target);
    }

    /**
     * Write the geli header to the target file.
     * @param header The header to write.
     */
    public void writeHeader(SAMFileHeader header) {
        ((GeliGenotypeWriter)writer).writeHeader(header);
    }

    /**
     * Writes the genotype likelihoods to the output.
     * @param gl genotype likelihoods to write.
     */
    public void addGenotypeLikelihoods(GenotypeLikelihoods gl) {
        ((GeliGenotypeWriter)writer).addGenotypeLikelihoods(gl);
    }


    /**
     * Merges the stream backing up this temporary storage into the target.
     * @param target Target stream for the temporary storage.  May not be null.
     */
    public void mergeInto(GeliGenotypeWriter target) {
        GeliFileReader reader = new GeliFileReader(file);
        while ( reader.hasNext() )
            target.addGenotypeLikelihoods(reader.next());
        reader.close();

        file.delete();        
    }
}

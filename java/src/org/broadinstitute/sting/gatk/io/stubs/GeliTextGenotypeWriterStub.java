package org.broadinstitute.sting.gatk.io.stubs;

import org.broadinstitute.sting.utils.genotype.geli.GeliGenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;

import java.io.File;

import net.sf.samtools.SAMFileHeader;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;

/**
 * Stub providing a passthrough for geli text files.
 *
 * @author mhanna
 * @version 0.1
 */
public class GeliTextGenotypeWriterStub extends GenotypeWriterStub<GeliGenotypeWriter> implements GeliGenotypeWriter {
    /**
     * Construct a new stub with the given engine and target file.
     * @param engine The engine, for extracting command-line arguments, etc.
     * @param genotypeFile Target file into which to write genotyping data.
     */
    public GeliTextGenotypeWriterStub(GenomeAnalysisEngine engine, File genotypeFile) {
        super(engine,genotypeFile);
    }

    /**
     * Gets the format of this stub.  We may want to discontinue use of this method and rely on instanceof comparisons.
     * @return GELI always.
     */
    public GenotypeWriterFactory.GENOTYPE_FORMAT getFormat() {
        return GenotypeWriterFactory.GENOTYPE_FORMAT.GELI;
    }

    /**
     * Write the geli header to the target file.
     * @param header The header to write.
     */
    public void writeHeader(SAMFileHeader header) {
        outputTracker.getStorage(this).writeHeader(header);
    }

    /**
     * Writes the genotype likelihoods to the output.
     * @param gl genotype likelihoods to write.
     */
    public void addGenotypeLikelihoods(GenotypeLikelihoods gl) {
        outputTracker.getStorage(this).addGenotypeLikelihoods(gl);
    }
}

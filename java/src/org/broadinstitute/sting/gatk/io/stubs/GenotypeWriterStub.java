/*
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.io.stubs;

import java.io.File;
import java.util.List;

import org.broadinstitute.sting.gatk.io.OutputTracker;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.VariationCall;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import net.sf.samtools.SAMFileHeader;

/**
 * A stub for routing and management of genotype reading and writing.
 *
 * @author ebanks
 * @version 0.1
 */
public abstract class GenotypeWriterStub<T extends GenotypeWriter> implements Stub<T>, GenotypeWriter {

    /**
     * Engine to use for collecting attributes for the output SAM file.
     */
    private final GenomeAnalysisEngine engine;

    /**
     * The file that this stub should write to.  Should be passed along to
     * whatever happens to create the StreamConnector.
     */
    private final File genotypeFile;

    /**
     * Connects this stub with an external stream capable of serving the
     * requests of the consumer of this stub.
     */
    protected OutputTracker outputTracker = null;

    /**
     * Create a new stub given the requested file.
     * @param engine        GATK engine.
     * @param genotypeFile  file to (ultimately) create.
     */
    public GenotypeWriterStub(GenomeAnalysisEngine engine,File genotypeFile) {
        this.engine = engine;
        this.genotypeFile = genotypeFile;
    }

    /**
     * Retrieves the file to (ultimately) be created.
     * @return The file.  Must not be null.
     */
    public File getFile() {
        return genotypeFile;
    }

    /**
     * Retrieves the header to use when creating the new file.
     * @return header to use when creating the new file.
     */
    public SAMFileHeader getSAMFileHeader() {
        return engine.getSAMFileHeader();
    }

    /**
     * Retrieves the format to use when creating the new file.
     * @return format to use when creating the new file.
     */
    public abstract GenotypeWriterFactory.GENOTYPE_FORMAT getFormat();

    /**
     * Registers the given streamConnector with this stub.
     * @param outputTracker The connector used to provide an appropriate stream.
     */
    public void register( OutputTracker outputTracker ) {
        this.outputTracker = outputTracker;
    }

    /**
     * @{inheritDoc}
     */
    public void addGenotypeCall(Genotype call) {
        outputTracker.getStorage(this).addGenotypeCall(call);
    }

    /**
     * @{inheritDoc}
     */
    public void addNoCall(int position) {
        outputTracker.getStorage(this).addNoCall(position);
    }

    /**
     * @{inheritDoc}
     */
    public void addMultiSampleCall(List<Genotype> genotypes, VariationCall variation) {
        outputTracker.getStorage(this).addMultiSampleCall(genotypes, variation);
    }

    /**
     * @{inheritDoc}
     */
    public boolean supportsMultiSample() {
        return outputTracker.getStorage(this).supportsMultiSample();
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        outputTracker.getStorage(this).close();
    }

}
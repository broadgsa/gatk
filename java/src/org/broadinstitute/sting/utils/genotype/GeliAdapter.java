package org.broadinstitute.sting.utils.genotype;

import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import edu.mit.broad.picard.genotype.geli.GeliFileWriter;
import net.sf.samtools.SAMFileHeader;

import java.io.File;


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

/**
 * @author aaron
 * @version 1.0
 *
 * Class GeliAdapter
 * Adapts the Geli file writer to the Genotype writer interface
 */
public class GeliAdapter implements GenotypeWriter {


    // the geli file writer we're adapting
    private final GeliFileWriter writer;

    /**
     * wrap a GeliFileWriter in the Genotype writer interface
     * @param writeTo where to write to
     * @param fileHeader the file header to write out
     */
    public GeliAdapter( File writeTo, final SAMFileHeader fileHeader ) {
        this.writer = GeliFileWriter.newInstanceForPresortedRecords(writeTo,fileHeader);
    }


    /**
     * add a single point genotype call to the
     *
     * @param position      the position on the contig
     * @param referenceBase the reference base
     * @param readDepth     the read depth at the specified position
     * @param likelihoods   the likelihoods of each of the possible alleles
     */
    @Override
    public void addGenotypeCall( int position, float rmsMapQuals, char referenceBase, int readDepth, LikelihoodObject likelihoods ) {
        writer.addGenotypeLikelihoods(likelihoods.convert(writer.getFileHeader(), 1, position, (byte)referenceBase));
    }

    /**
     * add a variable length call to the genotyper
     *
     * @param position      the position on the genome
     * @param rmsMapQuals   the root mean square of the mapping qualities
     * @param readDepth     the read depth
     * @param refBase       the reference base
     * @param firstHomZyg   the first homozygous indel
     * @param secondHomZyg  the second homozygous indel (if present, null if not)
     * @param hetLikelihood the heterozygous likelihood
     */
    @Override
    public void addVariableLengthCall( int position, float rmsMapQuals, int readDepth, char refBase, IndelLikelihood firstHomZyg, IndelLikelihood secondHomZyg, byte hetLikelihood ) {
        throw new UnsupportedOperationException("Geli format does not support variable length allele calls");
    }

    /**
     * add a no call to the genotype file, if supported.
     *
     * @param position
     * @param readDepth
     */
    @Override
    public void addNoCall( int position, int readDepth ) {
        throw new UnsupportedOperationException("Geli format does not support no-calls");
    }

    /** finish writing, closing any open files. */
    @Override
    public void close() {
        if (this.writer != null) {
            this.writer.close();
        }
    }
}

package org.broadinstitute.sting.utils.genotype.geli;

import edu.mit.broad.picard.genotype.geli.GeliFileWriter;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import net.sf.samtools.SAMFileHeader;
import org.broad.tribble.util.variantcontext.VariantContext;

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
 * @author aaron, ebanks
 * @version 1.0
 *          <p/>
 *          Class GeliAdapter
 *          Adapts the Geli file writer to the Genotype writer interface
 */
public class GeliAdapter implements GeliGenotypeWriter {

    // the file we're writing to
    private File writeTo = null;

    // the geli file writer we're adapting
    private GeliFileWriter writer = null;

    /**
     * wrap a GeliFileWriter in the Genotype writer interface
     *
     * @param writeTo    where to write to
     */
    public GeliAdapter(File writeTo) {
        this.writeTo = writeTo;
    }

    /**
     * wrap a GeliFileWriter in the Genotype writer interface
     *
     * @param fileHeader the file header to write out
     */
    public void writeHeader(final SAMFileHeader fileHeader) {
        this.writer = GeliFileWriter.newInstanceForPresortedRecords(writeTo, fileHeader);
    }

    public void addGenotypeLikelihoods(GenotypeLikelihoods gl) {
        if ( writer == null )
            throw new IllegalStateException("The Geli Header must be written before records can be added");

        writer.addGenotypeLikelihoods(gl);
    }

    /**
     * Add a genotype, given a variant context
     *
     * @param vc  the variant context representing the call to add
     * @param refBase not used by this writer
     */
    public void add(VariantContext vc, byte refBase) {
        throw new UnsupportedOperationException("We no longer support writing Geli");
    }

    /** finish writing, closing any open files. */
    public void close() {
        if (this.writer != null) {
            this.writer.close();
        }
    }
}

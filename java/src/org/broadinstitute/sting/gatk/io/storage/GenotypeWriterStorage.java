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

package org.broadinstitute.sting.gatk.io.storage;

import java.io.*;
import java.util.List;
import java.util.Set;
import java.util.HashSet;

import org.broadinstitute.sting.gatk.io.stubs.GenotypeWriterStub;
import org.broadinstitute.sting.utils.genotype.*;
import org.broadinstitute.sting.utils.genotype.glf.*;
import org.broadinstitute.sting.utils.genotype.geli.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.SampleUtils;
import edu.mit.broad.picard.genotype.geli.GeliFileReader;

/**
 * Provides temporary storage for GenotypeWriters.
 *
 * @author ebanks
 * @version 0.1
 */
public class GenotypeWriterStorage implements GenotypeWriter, Storage<GenotypeWriter> {
    private final File file;
    private final GenotypeWriter writer;

    public GenotypeWriterStorage( GenotypeWriterStub stub ) {
        this(stub, stub.getFile());   
    }

    public GenotypeWriterStorage( GenotypeWriterStub stub, File file ) {
        this.file = file;
        writer = GenotypeWriterFactory.create(stub.getFormat(), file);
        Set<String> samples = SampleUtils.getSAMFileSamples(stub.getSAMFileHeader());
        GenotypeWriterFactory.writeHeader(writer, stub.getSAMFileHeader(), samples, new HashSet<VCFHeaderLine>());
    }

    public void mergeInto( GenotypeWriter targetStream ) {

        // TODO -- This is ugly, but there is no GenotypeWriter interface since
        // TODO -- VCFReaders need to be separated out for compatability with Tribble
        // TODO -- and the adapters don't all implement a common interface.  Fix me.  Please.

        // VCF
        if ( targetStream instanceof VCFGenotypeWriterAdapter ) {
            VCFReader reader = new VCFReader(file);
            while ( reader.hasNext() )
                ((VCFGenotypeWriterAdapter)targetStream).addRecord(reader.next());
            reader.close();
        }

        // GELI TEXT
        else if ( targetStream instanceof GeliTextWriter ) {
            GeliFileReader reader = new GeliFileReader(file);
            while ( reader.hasNext() )
                ((GeliTextWriter)targetStream).addGenotypeLikelihoods(reader.next());
            reader.close();
        }

        // GELI BINARY
        else if ( targetStream instanceof GeliAdapter ) {
            GeliFileReader reader = new GeliFileReader(file);
            while ( reader.hasNext() )
                ((GeliAdapter)targetStream).addGenotypeLikelihoods(reader.next());
            reader.close();
        }

        // GLF
        else if ( targetStream instanceof GLFWriter ) {
            GLFReader reader = new GLFReader(file);
            while ( reader.hasNext() ) {
                 GLFRecord rec = reader.next();                    
                ((GLFWriter)targetStream).addGLFRecord(rec.getContig(),(int)rec.getPosition(),rec);
            }
            reader.close();
        }

        file.delete();
    }

    public void addGenotypeCall(Genotype call) {
        writer.addGenotypeCall(call);
    }

    public void addNoCall(int position) {
        writer.addNoCall(position);
    }

    public void addMultiSampleCall(List<Genotype> genotypes, VariationCall variation) {
        writer.addMultiSampleCall(genotypes, variation);
    }

    public boolean supportsMultiSample() {
        return writer.supportsMultiSample();
    }

    public void close() {
        writer.close();
    }

}
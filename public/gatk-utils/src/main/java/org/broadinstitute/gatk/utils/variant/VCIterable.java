/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.variant;

import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.FeatureCodecHeader;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.broadinstitute.gatk.utils.collections.Pair;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;

/*
* NOTE: Refactored out of GATKVCFUtils
*/
public class VCIterable<SOURCE> implements Iterable<VariantContext>, Iterator<VariantContext> {
    final SOURCE source;
    final FeatureCodec<VariantContext, SOURCE> codec;
    final VCFHeader header;

    VCIterable(final SOURCE source, final FeatureCodec<VariantContext, SOURCE> codec, final VCFHeader header) {
        this.source = source;
        this.codec = codec;
        this.header = header;
    }

    /**
     * Utility class to read all of the VC records from a file
     *
     * @param file
     * @param codec
     * @return
     * @throws java.io.IOException
     */
    public final static <SOURCE> Pair<VCFHeader, VCIterable<SOURCE>> readAllVCs( final File file, final FeatureCodec<VariantContext, SOURCE> codec) throws IOException {
        // read in the features
        SOURCE source = codec.makeSourceFromStream(new FileInputStream(file));
        FeatureCodecHeader header = codec.readHeader(source);
        final VCFHeader vcfHeader = (VCFHeader)header.getHeaderValue();
        return new Pair<>(vcfHeader, new VCIterable<>(source, codec, vcfHeader));
    }

    @Override
    public Iterator<VariantContext> iterator() {
        return this;
    }

    @Override
    public boolean hasNext() {
        return ! codec.isDone(source);
    }

    @Override
    public VariantContext next() {
        try {
            final VariantContext vc = codec.decode(source);
            return vc == null ? null : vc.fullyDecode(header, false);
        } catch ( IOException e ) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void remove() {
    }
}

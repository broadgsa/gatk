/*
 * Copyright (c) 2012, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.codecs.bcf2;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

/**
 * Lazy version of genotypes decoder for BCF2 genotypes
 *
 * @author Mark DePristo
 * @since 5/12
 */
class BCF2LazyGenotypesDecoder implements LazyGenotypesContext.LazyParser {
    final protected static Logger logger = Logger.getLogger(BCF2LazyGenotypesDecoder.class);

    // the essential information for us to use to decode the genotypes data
    // initialized when this lazy decoder is created, as we know all of this from the BCF2Codec
    // and its stored here again for code cleanliness
    private final BCF2Codec codec;
    private final ArrayList<Allele> siteAlleles;
    private final int nSamples;
    private final int nFields;

    BCF2LazyGenotypesDecoder(final BCF2Codec codec, final ArrayList<Allele> alleles, final int nSamples, final int nFields) {
        this.codec = codec;
        this.siteAlleles = alleles;
        this.nSamples = nSamples;
        this.nFields = nFields;
    }

    @Override
    public LazyGenotypesContext.LazyData parse(final Object data) {
        logger.info("Decoding BCF genotypes for " + nSamples + " samples with " + nFields + " fields each");

        // load our byte[] data into the decoder
        final BCF2Decoder decoder = new BCF2Decoder(((BCF2Codec.LazyData)data).bytes);

        // TODO -- fast path for sites only

        // go ahead and decode everyone
        final List<String> samples = new ArrayList<String>(codec.getHeader().getGenotypeSamples());

        if ( samples.size() != nSamples )
            throw new UserException.MalformedBCF2("GATK currently doesn't support reading BCF2 files with " +
                    "different numbers of samples per record.  Saw " + samples.size() +
                    " samples in header but have a record with " + nSamples + " samples");

        // create and initialize the genotypes array
        final ArrayList<GenotypeBuilder> builders = new ArrayList<GenotypeBuilder>(nSamples);
        for ( int i = 0; i < nSamples; i++ ) {
            builders.add(new GenotypeBuilder(samples.get(i)));
        }

        for ( int i = 0; i < nFields; i++ ) {
            // get the field name
            final int offset = (Integer) decoder.decodeTypedValue();
            final String field = codec.getDictionaryString(offset);

            // the type of each element
            final byte typeDescriptor = decoder.readTypeDescriptor();
            final BCF2GenotypeFieldDecoders.Decoder fieldDecoder = codec.getGenotypeFieldDecoder(field);
            try {
                fieldDecoder.decode(siteAlleles, field, decoder, typeDescriptor, builders);
            } catch ( ClassCastException e ) {
                throw new UserException.MalformedBCF2("BUG: expected encoding of field " + field
                        + " inconsistent with the value observed in the decoded value");
            }
        }

        final ArrayList<Genotype> genotypes = new ArrayList<Genotype>(nSamples);
        for ( final GenotypeBuilder gb : builders )
            genotypes.add(gb.make());

        return new LazyGenotypesContext.LazyData(genotypes, codec.getHeader().getSampleNamesInOrder(), codec.getHeader().getSampleNameToOffset());
    }
}

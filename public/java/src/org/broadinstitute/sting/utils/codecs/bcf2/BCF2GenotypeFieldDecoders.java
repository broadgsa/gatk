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

import com.google.java.contract.Ensures;
import com.google.java.contract.Requires;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.GenotypeBuilder;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * An efficient
 *
 * @author Your Name
 * @since Date created
 */
public class BCF2GenotypeFieldDecoders {
    // initialized once per writer to allow parallel writers to work
    private final HashMap<String, Decoder> genotypeFieldDecoder = new HashMap<String, Decoder>();
    private final Decoder defaultDecoder = new GenericDecoder();

    public BCF2GenotypeFieldDecoders(final VCFHeader header) {
        // TODO -- fill in appropriate decoders for each FORMAT field in the header

        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_KEY, new GTDecoder());
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_FILTER_KEY, new FLDecoder());
        genotypeFieldDecoder.put(VCFConstants.DEPTH_KEY, new DPDecoder());
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_ALLELE_DEPTHS, new ADDecoder());
        genotypeFieldDecoder.put(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY, new PLDecoder());
        genotypeFieldDecoder.put(VCFConstants.GENOTYPE_QUALITY_KEY, new GQDecoder());
    }

    // -----------------------------------------------------------------
    //
    // Genotype field decoder
    //
    // -----------------------------------------------------------------

    /**
     * Return decoder appropriate for field, or the generic decoder if no
     * specialized one is bound
     * @param field the GT field to decode
     * @return a non-null decoder
     */
    @Requires("field != null")
    @Ensures("result != null")
    public Decoder getDecoder(final String field) {
        final Decoder d = genotypeFieldDecoder.get(field);
        return d == null ? defaultDecoder : d;
    }

    /**
     * Decoder a field (implicit from creation) encoded as
     * typeDescriptor in the decoder object in the GenotypeBuilders
     * one for each sample in order.
     *
     * The way this works is that this decode method
     * iterates over the builders, decoding a genotype field
     * in BCF2 for each sample from decoder.
     *
     * This system allows us to easily use specialized
     * decoders for specific genotype field values. For example,
     * we use a special decoder to directly read the BCF2 data for
     * the PL field into a int[] rather than the generic List of Integer
     */
    public interface Decoder {
        public void decode(final List<Allele> siteAlleles,
                           final String field,
                           final BCF2Decoder decoder,
                           final byte typeDescriptor,
                           final List<GenotypeBuilder> gbs);
    }

    private class GTDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            for ( final GenotypeBuilder gb : gbs ) {
                // TODO -- fast path for size == 2 (diploid)
                final List<Integer> encoded = (List<Integer>)decoder.decodeTypedValue(typeDescriptor);
                if ( encoded == null )
                    // no called sample GT = .
                    gb.alleles(null);
                else {
                    // we have at least some alleles to decode
                    final List<Allele> gt = new ArrayList<Allele>(encoded.size());

                    for ( final Integer encode : encoded ) {
                        if ( encode == null ) {
                            // absent, as are all following by definition
                            break;
                        } else {
                            final int offset = encode >> 1;
                            if ( offset == 0 )
                                gt.add(Allele.NO_CALL);
                            else
                                gt.add(siteAlleles.get(offset - 1));
                        }
                    }

                    gb.alleles(gt);
                }
            }
        }
    }

    private class DPDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            for ( final GenotypeBuilder gb : gbs ) {
                final Object value = decoder.decodeTypedValue(typeDescriptor);
                if ( value != null )
                    gb.DP((Integer)value);
            }
        }
    }

    private class GQDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            for ( final GenotypeBuilder gb : gbs ) {
                final Object value = decoder.decodeTypedValue(typeDescriptor);
                if ( value != null )
                    gb.GQ((Integer)value);
            }
        }
    }

    private class ADDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            for ( final GenotypeBuilder gb : gbs ) {
                final int[] AD = decodeIntArray(decoder.decodeTypedValue(typeDescriptor));
                if ( AD != null )
                    gb.AD(AD);
            }
        }
    }

    private class PLDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            for ( final GenotypeBuilder gb : gbs ) {
                final int[] PL = decodeIntArray(decoder.decodeTypedValue(typeDescriptor));
                if ( PL != null )
                    gb.PL(PL);
            }
        }
    }

    private class FLDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            throw new ReviewedStingException("Genotype filter not implemented in BCF2 yet");
        }
    }

    private class GenericDecoder implements Decoder {
        @Override
        public void decode(final List<Allele> siteAlleles, final String field, final BCF2Decoder decoder, final byte typeDescriptor, final List<GenotypeBuilder> gbs) {
            for ( final GenotypeBuilder gb : gbs ) {
                Object value = decoder.decodeTypedValue(typeDescriptor);
                if ( value != null ) { // don't add missing values
                    if ( value instanceof List && ((List)value).size() == 1)
                        value = ((List)value).get(0);
                    gb.attribute(field, value);
                }
            }
        }
    }

    private static final int[] decodeIntArray(final Object value) {
        // todo -- decode directly into int[]
        final List<Integer> pls = (List<Integer>)value;
        if ( pls != null ) { // we have a PL field
            final int[] x = new int[pls.size()];
            for ( int j = 0; j < x.length; j++ )
                x[j] = pls.get(j);
            return x;
        } else
            return null;
    }
}

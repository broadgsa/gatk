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
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
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

//    @Override
//    public LazyGenotypesContext.LazyData parse(final Object data) {
//        logger.info("Decoding BCF genotypes for " + nSamples + " samples with " + nFields + " fields each");
//
//        // load our bytep[] data into the decoder
//        final BCF2Decoder decoder = new BCF2Decoder((byte[])data);
//
//        // go ahead and decode everyone
//        final List<String> samples = new ArrayList<String>(codec.getHeader().getGenotypeSamples());
//
//        if ( samples.size() != nSamples )
//            throw new UserException.MalformedBCF2("GATK currently doesn't support reading BCF2 files with " +
//                    "different numbers of samples per record.  Saw " + samples.size() +
//                    " samples in header but have a record with " + nSamples + " samples");
//
//        final Map<String, List<Object>> fieldValues = decodeGenotypeFieldValues(decoder, nFields, nSamples);
//        final ArrayList<Genotype> genotypes = new ArrayList<Genotype>(nSamples);
//        for ( int i = 0; i < nSamples; i++ ) {
//            // all of the information we need for each genotype, with default values
//            final String sampleName = samples.get(i);
//            List<Allele> alleles = null;
//            boolean isPhased = false;
//            double log10PError = VariantContext.NO_LOG10_PERROR;
//            Set<String> filters = null;
//            Map<String, Object> attributes = null;
//            double[] log10Likelihoods = null;
//
//            for ( final Map.Entry<String, List<Object>> entry : fieldValues.entrySet() ) {
//                final String field = entry.getKey();
//                Object value = entry.getValue().get(i);
//                try {
//                    if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
//                        alleles = decodeGenotypeAlleles(siteAlleles, (List<Integer>)value);
//                    } else if ( field.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
//                        if ( value != BCF2Type.INT8.getMissingJavaValue() )
//                            log10PError = ((Integer)value) / -10.0;
//                    } else if ( field.equals(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY) ) {
//                        final List<Integer> pls = (List<Integer>)value;
//                        if ( pls != null ) { // we have a PL field
//                            log10Likelihoods = new double[pls.size()];
//                            for ( int j = 0; j < log10Likelihoods.length; j++ ) {
//                                final double d = pls.get(j);
//                                log10Likelihoods[j] = d == -0.0 ? 0.0 : d / -10.0;
//                            }
//                        }
//                    } else if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
//                        throw new ReviewedStingException("Genotype filters not implemented in GATK BCF2");
//                        //filters = new HashSet<String>(values.get(i));
//                    } else { // add to attributes
//                        if ( value != null ) { // don't add missing values
//                            if ( attributes == null ) attributes = new HashMap<String, Object>(nFields);
//                            if ( value instanceof List && ((List)value).size() == 1)
//                                value = ((List)value).get(0);
//                            attributes.put(field, value);
//                        }
//                    }
//                } catch ( ClassCastException e ) {
//                    throw new UserException.MalformedBCF2("BUG: expected encoding of field " + field
//                            + " inconsistent with the value observed in the decoded value in the "
//                            + " BCF file.  Value was " + value);
//                }
//            }
//
//            if ( alleles == null ) throw new UserException.MalformedBCF2("BUG: no alleles found");
//
//            final Genotype g = new Genotype(sampleName, alleles, log10PError, filters, attributes, isPhased, log10Likelihoods);
//            genotypes.add(g);
//        }
//
//        return new LazyGenotypesContext.LazyData(genotypes, codec.getHeader().getSampleNamesInOrder(), codec.getHeader().getSampleNameToOffset());
//    }

    @Override
    public LazyGenotypesContext.LazyData parse(final Object data) {
        logger.info("Decoding BCF genotypes for " + nSamples + " samples with " + nFields + " fields each");

        // load our bytep[] data into the decoder
        final BCF2Decoder decoder = new BCF2Decoder(((BCF2Codec.LazyData)data).bytes);

        // go ahead and decode everyone
        final List<String> samples = new ArrayList<String>(codec.getHeader().getGenotypeSamples());

        if ( samples.size() != nSamples )
            throw new UserException.MalformedBCF2("GATK currently doesn't support reading BCF2 files with " +
                    "different numbers of samples per record.  Saw " + samples.size() +
                    " samples in header but have a record with " + nSamples + " samples");

        final Map<String, List<Object>> fieldValues = decodeGenotypeFieldValues(decoder, nFields, nSamples);
        final ArrayList<Genotype> genotypes = new ArrayList<Genotype>(nSamples);
        final GenotypeBuilder gb = new GenotypeBuilder();
        for ( int i = 0; i < nSamples; i++ ) {
            // all of the information we need for each genotype, with default values
            gb.reset();
            gb.name(samples.get(i));

            for ( final Map.Entry<String, List<Object>> entry : fieldValues.entrySet() ) {
                final String field = entry.getKey();
                Object value = entry.getValue().get(i);
                try {
                    if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
                        gb.alleles(decodeGenotypeAlleles(siteAlleles, (List<Integer>)value));
                    } else if ( field.equals(VCFConstants.DEPTH_KEY) ) {
                        if ( value != BCF2Type.INT8.getMissingJavaValue() )
                            gb.DP((Integer)value);
                    } else if ( field.equals(VCFConstants.GENOTYPE_QUALITY_KEY) ) {
                        if ( value != BCF2Type.INT8.getMissingJavaValue() )
                            gb.GQ((Integer)value);
                    } else if ( field.equals(VCFConstants.PHRED_GENOTYPE_LIKELIHOODS_KEY) ) {
                        final int[] PLs = decodeIntArray(value);
                        if ( PLs != null )
                            gb.PL(PLs);
                    } else if ( field.equals(VCFConstants.GENOTYPE_ALLELE_DEPTHS) ) {
                        final int[] AD = decodeIntArray(value);
                        if ( AD != null )
                            gb.AD(AD);
                    } else if ( field.equals(VCFConstants.GENOTYPE_FILTER_KEY) ) {
                        throw new ReviewedStingException("Genotype filters not implemented in GATK BCF2");
                        //filters = new HashSet<String>(values.get(i));
                    } else { // add to attributes
                        if ( value != null ) { // don't add missing values
                            if ( value instanceof List && ((List)value).size() == 1)
                                value = ((List)value).get(0);
                            gb.attribute(field, value);
                        }
                    }
                } catch ( ClassCastException e ) {
                    throw new UserException.MalformedBCF2("BUG: expected encoding of field " + field
                            + " inconsistent with the value observed in the decoded value in the "
                            + " BCF file.  Value was " + value);
                }
            }

            final Genotype g = gb.make();
            genotypes.add(g);
        }

        return new LazyGenotypesContext.LazyData(genotypes, codec.getHeader().getSampleNamesInOrder(), codec.getHeader().getSampleNameToOffset());
    }

    private final int[] decodeIntArray(final Object value) {
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

    private final List<Allele> decodeGenotypeAlleles(final ArrayList<Allele> siteAlleles, final List<Integer> encoded) {
        if ( encoded == null )
            // no called sample GT = .
            return Collections.emptyList();
        else {
            // we have at least some alleles to decode
            final List<Allele> gt = new ArrayList<Allele>(encoded.size());
            for ( final Integer encode : encoded ) {
                if ( encode == null ) // absent, as are all following by definition
                    return gt;
                else {
                    final int offset = encode >> 1;
                    if ( offset == 0 )
                        gt.add(Allele.NO_CALL);
                    else
                        gt.add(siteAlleles.get(offset - 1));
                }
            }

            return gt;
        }
    }

    private final Map<String, List<Object>> decodeGenotypeFieldValues(final BCF2Decoder decoder,
                                                                      final int nFields,
                                                                      final int nSamples) {
        assert (nFields > 0 && nSamples > 0) || (nFields == 0 && nSamples == 0);

        if ( nFields == 0 ) // fast path exit for sites only file
            return Collections.emptyMap();
        else {
            final Map<String, List<Object>> map = new LinkedHashMap<String, List<Object>>(nFields);

            for ( int i = 0; i < nFields; i++ ) {
                final int offset = (Integer) decoder.decodeTypedValue();
                final String field = codec.getDictionaryString(offset);

                // the type of each element
                final byte typeDescriptor = decoder.readTypeDescriptor();
                final List<Object> values = new ArrayList<Object>(nSamples);
                for ( int j = 0; j < nSamples; j++ )
                    values.add(decoder.decodeTypedValue(typeDescriptor));

                assert ! map.containsKey(field);

                map.put(field, values);
            }

            return map;
        }
    }
}

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

package org.broadinstitute.sting.utils.variantcontext.writer;

import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Encoder;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Type;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Utils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * [Short one sentence description of this walker]
 * <p/>
 * <p>
 * [Functionality of this walker]
 * </p>
 * <p/>
 * <h2>Input</h2>
 * <p>
 * [Input description]
 * </p>
 * <p/>
 * <h2>Output</h2>
 * <p>
 * [Output description]
 * </p>
 * <p/>
 * <h2>Examples</h2>
 * <pre>
 *    java
 *      -jar GenomeAnalysisTK.jar
 *      -T $WalkerName
 *  </pre>
 *
 * @author Your Name
 * @since Date created
 */
public abstract class BCF2FieldWriter {
    private final VCFHeader header;
    private final BCF2FieldEncoder fieldEncoder;

    protected BCF2FieldWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
        this.header = header;
        this.fieldEncoder = fieldEncoder;
    }

    protected VCFHeader getHeader() { return header; }
    protected BCF2FieldEncoder getFieldEncoder() {
        return fieldEncoder;
    }
    protected String getField() { return getFieldEncoder().getField(); }

    public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
        encoder.encodeTyped(fieldEncoder.getDictionaryOffset(), fieldEncoder.getDictionaryOffsetType());
    }

    public void done(final BCF2Encoder encoder, final VariantContext vc) throws IOException { } // TODO -- overload done so that we null out values and test for correctness

    @Override
    public String toString() {
        return "BCF2FieldWriter " + getClass().getSimpleName() + " with encoder " + getFieldEncoder();
    }

    // --------------------------------------------------------------------------------
    //
    // Sites writers
    //
    // --------------------------------------------------------------------------------

    public static abstract class SiteWriter extends BCF2FieldWriter {
        protected SiteWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        public abstract void site(final BCF2Encoder encoder, final VariantContext vc) throws IOException;
    }

    public static class GenericSiteWriter extends SiteWriter {
        public GenericSiteWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        @Override
        public void site(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            final Object rawValue = vc.getAttribute(getField(), null);
            final BCF2Type type = getFieldEncoder().getType(rawValue);
            if ( rawValue == null ) {
                // the value is missing, just write in null
                encoder.encodeType(0, type);
            } else {
                final int valueCount = getFieldEncoder().getBCFFieldCount(vc, rawValue);
                encoder.encodeType(valueCount, type);
                getFieldEncoder().encodeValue(encoder, rawValue, type);
            }
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Genotypes writers
    //
    // --------------------------------------------------------------------------------

    public static abstract class GenotypesWriter extends BCF2FieldWriter {
        int nValuesPerGenotype = -1;
        BCF2Type encodingType = null;

        protected GenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);

            if ( fieldEncoder.hasFixedCount() ) {
                nValuesPerGenotype = getFieldEncoder().getFixedCount();
            }
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            // writes the key information
            super.start(encoder, vc);

            // only update if we need to
            if ( ! getFieldEncoder().hasFixedCount() ) {
                if ( getFieldEncoder().hasContextDeterminedCount() )
                    // we are cheap -- just depends on genotype of allele counts
                    nValuesPerGenotype = getFieldEncoder().getContextDeterminedCount(vc);
                else
                    // we have to go fishing through the values themselves (expensive)
                    nValuesPerGenotype = computeMaxSizeOfGenotypeFieldFromValues(vc);
            }

            encoder.encodeType(nValuesPerGenotype, encodingType);
        }

        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            final Object fieldValue = g.getAttribute(getField(), null);
            getFieldEncoder().encodeValue(encoder, fieldValue, encodingType, nValuesPerGenotype);
        }

        public Object getGenotypeValue(final Genotype g) {
            return g.getAttribute(getField());
        }

        private final int computeMaxSizeOfGenotypeFieldFromValues(final VariantContext vc) {
            int size = -1;

            for ( final Genotype g : vc.getGenotypes() ) {
                final Object o = getGenotypeValue(g);
                size = Math.max(size, getFieldEncoder().getBCFFieldCount(vc, o));
            }

            return size;
        }
    }

    public static class FixedTypeGenotypesWriter extends GenotypesWriter {
        public FixedTypeGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);

            encodingType = getFieldEncoder().getFixedType();
        }
    }

    public static class IntegerTypeGenotypesWriter extends GenotypesWriter {
        public IntegerTypeGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            // the only value that is dynamic are integers
            final List<Integer> values = new ArrayList<Integer>(vc.getNSamples());
            for ( final Genotype g : vc.getGenotypes() ) {
                for ( final Object i : BCF2Utils.toList(g.getAttribute(getField(), null)) ) {
                    values.add((Integer)i); // we know they are all integers
                }
            }

            encodingType = BCF2Utils.determineIntegerType(values);
            super.start(encoder, vc);
        }
    }

    // TODO TODO TODO TODO TODO
    // TODO
    // TODO THIS ROUTINE NEEDS TO BE OPTIMIZED.  IT ACCOUNTS FOR A SIGNIFICANT AMOUNT OF THE
    // TODO RUNTIME FOR WRITING OUT BCF FILES WITH MANY GENOTYPES
    // TODO
    // TODO TODO TODO TODO TODO
    public static class IGFGenotypesWriter extends GenotypesWriter {
        final IntGenotypeFieldAccessors.Accessor ige;

        public IGFGenotypesWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder, final IntGenotypeFieldAccessors.Accessor ige) {
            super(header, fieldEncoder);
            this.ige = ige;

            if ( ! (fieldEncoder instanceof BCF2FieldEncoder.IntArray) )
                throw new ReviewedStingException("BUG: IntGenotypesWriter requires IntArray encoder for field " + getField());
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            encodingType = BCF2Type.INT8;
            for ( final Genotype g : vc.getGenotypes() ) {
                final int[] pls = ige.getValues(g);
                final BCF2Type plsType = getFieldEncoder().getType(pls);
                encodingType = BCF2Utils.maxIntegerType(encodingType, plsType);
                if ( encodingType == BCF2Type.INT32 )
                    break; // stop early
            }

            super.start(encoder, vc);
        }

        @Override
        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            getFieldEncoder().encodeValue(encoder, ige.getValues(g), encodingType, nValuesPerGenotype);
        }

        @Override
        public Object getGenotypeValue(final Genotype g) {
            return ige.getValues(g);
        }
    }

    // TODO TODO TODO TODO TODO
    // TODO
    // TODO we should really have a fast path for encoding diploid genotypes where
    // TODO we don't pay the overhead of creating the allele maps
    // TODO
    // TODO TODO TODO TODO TODO
    public static class GTWriter extends GenotypesWriter {
        Map<Allele, Integer> alleleMap = null;

        public GTWriter(final VCFHeader header, final BCF2FieldEncoder fieldEncoder) {
            super(header, fieldEncoder);
        }

        @Override
        public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            if ( vc.getNAlleles() > BCF2Utils.MAX_ALLELES_IN_GENOTYPES )
                throw new ReviewedStingException("Current BCF2 encoder cannot handle sites " +
                        "with > " + BCF2Utils.MAX_ALLELES_IN_GENOTYPES + " alleles, but you have "
                        + vc.getNAlleles() + " at " + vc.getChr() + ":" + vc.getStart());

            encodingType = BCF2Type.INT8;
            alleleMap = buildAlleleMap(vc);
            nValuesPerGenotype = vc.getMaxPloidy();
            super.start(encoder, vc);    //To change body of overridden methods use File | Settings | File Templates.
        }

        @Override
        public void addGenotype(final BCF2Encoder encoder, final VariantContext vc, final Genotype g) throws IOException {
            final List<Allele> alleles = g.getAlleles();
            final int samplePloidy = alleles.size();
            for ( int i = 0; i < nValuesPerGenotype; i++ ) {
                if ( i < samplePloidy ) {
                    // we encode the actual allele
                    final Allele a = alleles.get(i);
                    final int offset = alleleMap.get(a);
                    final int encoded = ((offset+1) << 1) | (g.isPhased() ? 0x01 : 0x00);
                    encoder.encodePrimitive(encoded, encodingType);
                } else {
                    // we need to pad with missing as we have ploidy < max for this sample
                    encoder.encodePrimitive(encodingType.getMissingBytes(), encodingType);
                }
            }
        }

        private final static Map<Allele, Integer> buildAlleleMap(final VariantContext vc) {
            final Map<Allele, Integer> alleleMap = new HashMap<Allele, Integer>(vc.getAlleles().size()+1);
            alleleMap.put(Allele.NO_CALL, -1); // convenience for lookup

            final List<Allele> alleles = vc.getAlleles();
            for ( int i = 0; i < alleles.size(); i++ ) {
                alleleMap.put(alleles.get(i), i);
            }

            return alleleMap;
        }
    }
}


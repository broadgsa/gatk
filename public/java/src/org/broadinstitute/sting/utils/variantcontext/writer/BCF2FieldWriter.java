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
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.io.IOException;

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
    private final BCF2FieldEncoder fieldEncoder;

    protected BCF2FieldWriter(final BCF2FieldEncoder fieldEncoder) {
        this.fieldEncoder = fieldEncoder;
    }

    protected BCF2FieldEncoder getFieldEncoder() {
        return fieldEncoder;
    }

    public void start(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
        encoder.encodeTyped(fieldEncoder.getDictionaryOffset(), fieldEncoder.getDictionaryOffsetType());
    }

    public void done(final BCF2Encoder encoder, final VariantContext vc) throws IOException { }

    @Override
    public String toString() {
        return "BCF2FieldWriter " + getClass().getSimpleName() + " with encoder " + getFieldEncoder();
    }

    public static abstract class SiteWriter extends BCF2FieldWriter {
        protected SiteWriter(final BCF2FieldEncoder fieldEncoder) {
            super(fieldEncoder);
        }

        public abstract void site(final BCF2Encoder encoder, final VariantContext vc) throws IOException;
    }

    public static class GenericSiteWriter extends SiteWriter {
        public GenericSiteWriter(final BCF2FieldEncoder fieldEncoder) {
            super(fieldEncoder);
        }

        @Override
        public void site(final BCF2Encoder encoder, final VariantContext vc) throws IOException {
            final Object rawValue = vc.getAttribute(getFieldEncoder().getField(), null);
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
}


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

import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.codecs.bcf2.BCF2Encoder;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.util.HashMap;
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
public class BCF2FieldWriterManager {
    final protected static Logger logger = Logger.getLogger(BCF2FieldWriterManager.class);
    final Map<String, BCF2FieldWriter.SiteWriter> siteWriters = new HashMap<String, BCF2FieldWriter.SiteWriter>();
    final Map<String, BCF2FieldWriter.GenotypesWriter> genotypesWriters = new HashMap<String, BCF2FieldWriter.GenotypesWriter>();
    final IntGenotypeFieldAccessors intGenotypeFieldAccessors = new IntGenotypeFieldAccessors();

    public BCF2FieldWriterManager() { }

    public void setup(final VCFHeader header, final BCF2Encoder encoder, final Map<String, Integer> dictionary) {
        for (final VCFHeaderLine line : header.getMetaData()) {
            if ( line instanceof VCFInfoHeaderLine ) {
                final String field = ((VCFInfoHeaderLine) line).getID();
                final BCF2FieldWriter.SiteWriter writer = createInfoWriter(header, (VCFInfoHeaderLine)line, encoder, dictionary);
                log(field, writer);
                siteWriters.put(field, writer);
            } else if ( line instanceof VCFFormatHeaderLine ) {
                final String field = ((VCFFormatHeaderLine) line).getID();
                final BCF2FieldWriter.GenotypesWriter writer = createGenotypesWriter(header, (VCFFormatHeaderLine)line, encoder, dictionary);
                log(field, writer);
                genotypesWriters.put(field, writer);
            }
        }
    }

    private final void log(final String field, final BCF2FieldWriter writer) {
        logger.info("Using writer " + writer);
    }

    // -----------------------------------------------------------------
    //
    // Master routine to look at the header, a specific line, and
    // build an appropriate SiteWriter for that header element
    //
    // -----------------------------------------------------------------

    private BCF2FieldWriter.SiteWriter createInfoWriter(final VCFHeader header,
                                                        final VCFInfoHeaderLine line,
                                                        final BCF2Encoder encoder,
                                                        final Map<String, Integer> dict) {
        return new BCF2FieldWriter.GenericSiteWriter(header, createFieldEncoder(line, encoder, dict, false));
    }

    private BCF2FieldEncoder createFieldEncoder(final VCFCompoundHeaderLine line,
                                                final BCF2Encoder encoder,
                                                final Map<String, Integer> dict,
                                                final boolean createGenotypesEncoders ) {

        if ( createGenotypesEncoders && intGenotypeFieldAccessors.getAccessor(line.getID()) != null ) {
            if ( line.getType() != VCFHeaderLineType.Integer )
                logger.warn("Warning: field " + line.getID() + " expected to encode an integer but saw " + line.getType() + " for record " + line);
            return new BCF2FieldEncoder.IntArray(line, encoder, dict);
        } else if ( createGenotypesEncoders && line.getID().equals(VCFConstants.GENOTYPE_KEY) ) {
            return new BCF2FieldEncoder.IntList(line, encoder, dict);
        } else {
            switch ( line.getType() ) {
                case Character:
                case String:
                    return new BCF2FieldEncoder.StringOrCharacter(line, encoder, dict);
                case Flag:
                    return new BCF2FieldEncoder.Flag(line, encoder, dict);
                case Float:
                    return new BCF2FieldEncoder.Float(line, encoder, dict);
                case Integer:
                    if ( line.getCountType() == VCFHeaderLineCount.INTEGER && line.getCount() == 1 )
                        return new BCF2FieldEncoder.AtomicInt(line, encoder, dict);
                    else
                        return new BCF2FieldEncoder.IntList(line, encoder, dict);
                default:
                    throw new ReviewedStingException("Unexpected type for field " + line.getID());
            }
        }
    }

    // -----------------------------------------------------------------
    //
    // Master routine to look at the header, a specific line, and
    // build an appropriate Genotypes for that header element
    //
    // -----------------------------------------------------------------

    private BCF2FieldWriter.GenotypesWriter createGenotypesWriter(final VCFHeader header,
                                                                  final VCFFormatHeaderLine line,
                                                                  final BCF2Encoder encoder,
                                                                  final Map<String, Integer> dict) {
        final String field = line.getID();
        final BCF2FieldEncoder fieldEncoder = createFieldEncoder(line, encoder, dict, true);

        if ( field.equals(VCFConstants.GENOTYPE_KEY) ) {
            return new BCF2FieldWriter.GTWriter(header, fieldEncoder);
        } else if ( intGenotypeFieldAccessors.getAccessor(field) != null ) {
            return new BCF2FieldWriter.IGFGenotypesWriter(header, fieldEncoder, intGenotypeFieldAccessors.getAccessor(field));
        } else if ( line.getType() == VCFHeaderLineType.Integer ) {
            return new BCF2FieldWriter.IntegerTypeGenotypesWriter(header, fieldEncoder);
        } else {
            return new BCF2FieldWriter.FixedTypeGenotypesWriter(header, fieldEncoder);
        }
    }

    // -----------------------------------------------------------------
    //
    // Accessors to get site / genotype writers
    //
    // -----------------------------------------------------------------

    public BCF2FieldWriter.SiteWriter getSiteFieldWriter(final String key) {
        return getWriter(key, siteWriters);
    }

    public BCF2FieldWriter.GenotypesWriter getGenotypeFieldWriter(final String key) {
        return getWriter(key, genotypesWriters);
    }

    public <T> T getWriter(final String key, final Map<String, T> map) {
        final T writer = map.get(key);
        if ( writer == null ) throw new ReviewedStingException("BUG: no writer found for " + key);
        return writer;
    }
}

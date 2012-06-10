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
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineCount;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
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

    public BCF2FieldWriterManager() { }

    public void setup(final VCFHeader header, final BCF2Encoder encoder, final Map<String, Integer> dictionary) {
        for (final VCFHeaderLine line : header.getMetaData()) {
            if ( line instanceof VCFInfoHeaderLine ) {
                final String field = ((VCFInfoHeaderLine) line).getID();
                final BCF2FieldWriter.SiteWriter writer = createInfoWriter((VCFInfoHeaderLine)line, encoder, dictionary);
                logger.info("Installing for field " + field + " field writer " + writer);
                siteWriters.put(field, writer);
            }
        }
    }

    private BCF2FieldWriter.SiteWriter createInfoWriter(final VCFInfoHeaderLine line, final BCF2Encoder encoder, final Map<String, Integer> dict) {
        BCF2FieldEncoder fieldEncoder = null;
        switch ( line.getType() ) {
            case Character:
            case String:
                fieldEncoder = new BCF2FieldEncoder.StringOrCharacter(line, encoder, dict);
                break;
            case Flag:
                fieldEncoder = new BCF2FieldEncoder.Flag(line, encoder, dict);
                break;
            case Float:
                fieldEncoder = new BCF2FieldEncoder.Float(line, encoder, dict);
                break;
            case Integer:
                if ( line.getCountType() == VCFHeaderLineCount.INTEGER && line.getCount() == 1 )
                    fieldEncoder = new BCF2FieldEncoder.AtomicInt(line, encoder, dict);
                else
                    fieldEncoder = new BCF2FieldEncoder.IntList(line, encoder, dict);
                break;
            default:
                throw new ReviewedStingException("Unexpected type for field " + line.getID());
        }

        return new BCF2FieldWriter.GenericSiteWriter(fieldEncoder);
    }

    public BCF2FieldWriter.SiteWriter getSiteFieldWriter(final String key) {
        final BCF2FieldWriter.SiteWriter writer = siteWriters.get(key);
        if ( writer == null ) throw new ReviewedStingException("BUG: no writer found for " + key);
        return writer;
    }
}

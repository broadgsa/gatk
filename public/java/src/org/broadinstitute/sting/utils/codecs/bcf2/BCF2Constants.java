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

import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;

import java.nio.charset.Charset;

public class BCF2Constants {
    public static final String VERSION_LINE_FORMAT = "fileformat=BCF2v%d.%d";
    public static final String VERSION_LINE = String.format(VCFHeader.METADATA_INDICATOR + VERSION_LINE_FORMAT, 0, 1);
    public static final String DICTIONARY_LINE_TAG = "dictionary";
    public static final String DICTIONARY_LINE_ENTRY_SEPARATOR = ",";

    public static final Charset BCF2_TEXT_CHARSET = Charset.forName("US-ASCII");  // TODO: enforce this!

    // Note that these values are prefixed by FFFFFF for convenience
    public static final int INT8_MISSING_VALUE  = 0xFFFFFF80;
    public static final int INT16_MISSING_VALUE = 0xFFFF8000;
    public static final int INT32_MISSING_VALUE = 0x80000000;
    public static final int FLOAT_MISSING_VALUE = 0x7F800001;
}

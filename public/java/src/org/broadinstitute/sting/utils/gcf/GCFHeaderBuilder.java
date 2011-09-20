/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.utils.gcf;

import org.broadinstitute.sting.utils.variantcontext.Allele;

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
public class GCFHeaderBuilder {
    Map<Allele, Integer> alleles = new HashMap<Allele, Integer>();
    Map<String, Integer> strings = new HashMap<String, Integer>();
    Map<String, Integer> samples = new HashMap<String, Integer>();

    public GCFHeader createHeader() {
        return new GCFHeader(alleles, strings, samples);
    }

    public int encodeString(final String chr)    { return encode(strings, chr); }
    public int encodeAllele(final Allele allele) { return encode(alleles, allele); }
    public int encodeSample(final String sampleName) { return encode(samples, sampleName); }

    private <T> int encode(Map<T, Integer> map, T key) {
        Integer v = map.get(key);
        if ( v == null ) {
            v = map.size();
            map.put(key, v);
        }
        return v;
    }
}

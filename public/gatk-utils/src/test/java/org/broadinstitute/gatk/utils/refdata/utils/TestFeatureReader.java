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

package org.broadinstitute.gatk.utils.refdata.utils;

import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.TribbleIndexedFeatureReader;

import java.io.IOException;

/**
 * Feature reader with additional test utilities. The iterators can be checked to see if they are closed.
 */
public class TestFeatureReader extends TribbleIndexedFeatureReader<Feature, Object> {
    public TestFeatureReader(String featurePath, FeatureCodec codec) throws IOException {
        super(featurePath, codec, true);
    }

    @Override
    @SuppressWarnings("unchecked")
    public CheckableCloseableTribbleIterator<Feature> iterator() throws IOException {
        return new CheckableCloseableTribbleIterator<Feature>(super.iterator());
    }

    @Override
    @SuppressWarnings("unchecked")
    public CheckableCloseableTribbleIterator<Feature> query(String chr, int start, int end) throws IOException {
        return new CheckableCloseableTribbleIterator<Feature>(super.query(chr, start, end));
    }
}

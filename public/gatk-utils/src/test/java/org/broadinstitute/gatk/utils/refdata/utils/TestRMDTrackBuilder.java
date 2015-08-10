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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.index.Index;
import org.broadinstitute.gatk.utils.refdata.tracks.FeatureManager;
import org.broadinstitute.gatk.utils.refdata.tracks.IndexDictionaryUtils;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrack;
import org.broadinstitute.gatk.utils.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import java.io.File;
import java.io.IOException;

/**
 * Extension of RMDTrackBuilder that creates TestFeatureReader's which in turn create CheckableCloseableTribbleIterator's.
 */
public class TestRMDTrackBuilder extends RMDTrackBuilder {
    private GenomeLocParser genomeLocParser;

    public TestRMDTrackBuilder(SAMSequenceDictionary dict, GenomeLocParser genomeLocParser) {
        // disable auto-index creation/locking in the RMDTrackBuilder for tests
        super(dict, genomeLocParser, null, true, null);
        this.genomeLocParser = genomeLocParser;
    }

    @Override
    public RMDTrack createInstanceOfTrack(RMDTriplet fileDescriptor) {
        String name = fileDescriptor.getName();
        File inputFile = new File(fileDescriptor.getFile());
        FeatureManager.FeatureDescriptor descriptor = getFeatureManager().getByTriplet(fileDescriptor);
        FeatureCodec codec = getFeatureManager().createCodec(descriptor, name, genomeLocParser, null);
        TestFeatureReader featureReader;
        Index index;
        try {
            // Create a feature reader that creates checkable tribble iterators.
            index = loadIndex(inputFile, codec);
            featureReader = new TestFeatureReader(inputFile.getAbsolutePath(), codec);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        SAMSequenceDictionary sequenceDictionary = IndexDictionaryUtils.getSequenceDictionaryFromProperties(index);
        return new RMDTrack(descriptor.getCodecClass(), name, inputFile, featureReader, sequenceDictionary, genomeLocParser, codec);
    }
}

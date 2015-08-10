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

package org.broadinstitute.gatk.engine.datasources.reads;

import org.broadinstitute.gatk.utils.sam.SAMReaderID;
import org.broadinstitute.gatk.engine.resourcemanagement.ThreadAllocation;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.GenomeLocParser;

import java.util.List;
import java.util.Collections;

/**
 * A mock locus shard, usable for infrastructure that requires a shard to behave properly.
 *
 * @author mhanna
 * @version 0.1
 */
public class MockLocusShard extends LocusShard {
    public MockLocusShard(final GenomeLocParser genomeLocParser,final List<GenomeLoc> intervals) {
        super(  genomeLocParser,
                new SAMDataSource(null, Collections.<SAMReaderID>emptyList(),new ThreadAllocation(),null,genomeLocParser),
                intervals,
                null);
    }
}

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

package org.broadinstitute.gatk.utils;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;

import java.io.File;

/**
 * Caliper microbenchmark of genome loc parser
 */
public class GenomeLocParserBenchmark extends SimpleBenchmark {
    private IndexedFastaSequenceFile seq;
    private final int ITERATIONS = 1000000;

    @Param({"NEW", "NONE"})
    GenomeLocParser.ValidationLevel validationLevel; // set automatically by framework

    @Param({"true", "false"})
    boolean useContigIndex; // set automatically by framework

    @Override protected void setUp() throws Exception {
        seq = new CachingIndexedFastaSequenceFile(new File("/Users/depristo/Desktop/broadLocal/localData/human_g1k_v37.fasta"));
    }
//
//    public void timeSequentialCreationFromGenomeLoc(int rep) {
//        final GenomeLocParser genomeLocParser = new GenomeLocParser(seq.getSequenceDictionary(), validationLevel);
//        GenomeLoc last = genomeLocParser.createGenomeLoc("1", 1, 1);
//        for ( int i = 0; i < rep; i++ ) {
//            for ( int j = 1; j < ITERATIONS; j++ ) {
//                if ( useContigIndex )
//                    last = genomeLocParser.createGenomeLoc(last.getContig(), last.getContigIndex(), last.getStart() + 1);
//                else
//                    last = genomeLocParser.createGenomeLoc(last.getContig(), last.getStart() + 1);
//            }
//        }
//    }
//
//    public void timeSequentialCreationFromGenomeLocOriginal(int rep) {
//        final GenomeLocParserOriginal genomeLocParser = new GenomeLocParserOriginal(seq.getSequenceDictionary());
//        GenomeLoc last = genomeLocParser.createGenomeLoc("1", 1, 1);
//        for ( int i = 0; i < rep; i++ ) {
//            for ( int j = 1; j < ITERATIONS; j++ ) {
//                if ( useContigIndex )
//                    last = genomeLocParser.createGenomeLoc(last.getContig(), last.getContigIndex(), last.getStart() + 1);
//                else
//                    last = genomeLocParser.createGenomeLoc(last.getContig(), last.getStart() + 1);
//            }
//        }
//    }

    public static void main(String[] args) {
        com.google.caliper.Runner.main(GenomeLocParserBenchmark.class, args);
    }
}

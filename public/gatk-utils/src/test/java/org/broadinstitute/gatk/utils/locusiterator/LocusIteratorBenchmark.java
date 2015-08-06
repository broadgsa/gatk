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

package org.broadinstitute.gatk.utils.locusiterator;

import com.google.caliper.Param;
import com.google.caliper.SimpleBenchmark;
import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.GenomeLocParser;
import org.broadinstitute.gatk.utils.QualityUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecord;

import java.util.LinkedList;
import java.util.List;

/**
 * Caliper microbenchmark of fragment pileup
 */
public class LocusIteratorBenchmark extends SimpleBenchmark {
    protected SAMFileHeader header;
    protected GenomeLocParser genomeLocParser;

    List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
    final int readLength = 101;
    final int nReads = 10000;
    final int locus = 1;

    @Param({"101M", "50M10I40M", "50M10D40M"})
    String cigar; // set automatically by framework

    @Override protected void setUp() {
        header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());

        for ( int j = 0; j < nReads; j++ ) {
            GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, readLength);
            read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
            final byte[] quals = new byte[readLength];
            for ( int i = 0; i < readLength; i++ )
                quals[i] = (byte)(i % QualityUtils.MAX_SAM_QUAL_SCORE);
            read.setBaseQualities(quals);
            read.setCigarString(cigar);
            reads.add(read);
        }
    }

//    public void timeOriginalLIBS(int rep) {
//        for ( int i = 0; i < rep; i++ ) {
//            final org.broadinstitute.gatk.utils.locusiterator.old.LocusIteratorByState libs =
//                    new org.broadinstitute.gatk.utils.locusiterator.old.LocusIteratorByState(
//                            new LocusIteratorByStateBaseTest.FakeCloseableIterator<SAMRecord>(reads.iterator()),
//                            LocusIteratorByStateBaseTest.createTestReadProperties(),
//                            genomeLocParser,
//                            LocusIteratorByState.sampleListForSAMWithoutReadGroups());
//
//            while ( libs.hasNext() ) {
//                AlignmentContext context = libs.next();
//            }
//        }
//    }
//
//    public void timeLegacyLIBS(int rep) {
//        for ( int i = 0; i < rep; i++ ) {
//            final org.broadinstitute.gatk.utils.locusiterator.legacy.LegacyLocusIteratorByState libs =
//                    new org.broadinstitute.gatk.utils.locusiterator.legacy.LegacyLocusIteratorByState(
//                            new LocusIteratorByStateBaseTest.FakeCloseableIterator<SAMRecord>(reads.iterator()),
//                            LocusIteratorByStateBaseTest.createTestReadProperties(),
//                            genomeLocParser,
//                            LocusIteratorByState.sampleListForSAMWithoutReadGroups());
//
//            while ( libs.hasNext() ) {
//                AlignmentContext context = libs.next();
//            }
//        }
//    }

    public void timeNewLIBS(int rep) {
        for ( int i = 0; i < rep; i++ ) {
            final org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState libs =
                    new org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState(
                            new LocusIteratorByStateBaseTest.FakeCloseableIterator<GATKSAMRecord>(reads.iterator()),
                            null, true, false,
                            genomeLocParser,
                            LocusIteratorByState.sampleListForSAMWithoutReadGroups());

            while ( libs.hasNext() ) {
                AlignmentContext context = libs.next();
            }
        }
    }

//    public void timeOriginalLIBSStateMachine(int rep) {
//        for ( int i = 0; i < rep; i++ ) {
//            for ( final SAMRecord read : reads ) {
//                final SAMRecordAlignmentState alignmentStateMachine = new SAMRecordAlignmentState(read);
//                while ( alignmentStateMachine.stepForwardOnGenome() != null ) {
//                    alignmentStateMachine.getGenomeOffset();
//                }
//            }
//        }
//    }

    public void timeAlignmentStateMachine(int rep) {
        for ( int i = 0; i < rep; i++ ) {
            for ( final GATKSAMRecord read : reads ) {
                final AlignmentStateMachine alignmentStateMachine = new AlignmentStateMachine(read);
                while ( alignmentStateMachine.stepForwardOnGenome() != null ) {
                    ;
                }
            }
        }
    }

    public static void main(String[] args) {
        com.google.caliper.Runner.main(LocusIteratorBenchmark.class, args);
    }
}

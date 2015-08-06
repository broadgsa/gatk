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

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecordIterator;
import org.apache.log4j.Logger;
import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.CommandLineProgram;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.sam.GATKSAMRecordIterator;
import org.broadinstitute.gatk.utils.*;
import org.broadinstitute.gatk.utils.fasta.CachingIndexedFastaSequenceFile;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Caliper microbenchmark of fragment pileup
 */
public class LIBSPerformance extends CommandLineProgram {
    private static Logger logger = Logger.getLogger(LIBSPerformance.class);

    @Input(fullName = "input_file", shortName = "I", doc = "SAM or BAM file(s)", required = true)
    public File samFile = null;

    @Input(fullName = "reference_sequence", shortName = "R", doc = "Reference sequence file", required = true)
    public File referenceFile = null;

    @Argument(fullName = "L", shortName = "L", doc = "Query location", required = false)
    public String location = null;

    @Argument(fullName = "dt", shortName = "dt", doc = "Enable downsampling", required = false)
    public boolean downsample = false;

    @Override
    public int execute() throws IOException {
        final IndexedFastaSequenceFile reference = new CachingIndexedFastaSequenceFile(referenceFile);
        final GenomeLocParser genomeLocParser = new GenomeLocParser(reference);

        final SAMFileReader reader = new SAMFileReader(samFile);

        SAMRecordIterator rawIterator;
        if ( location == null )
            rawIterator = reader.iterator();
        else {
            final GenomeLoc loc = genomeLocParser.parseGenomeLoc(location);
            rawIterator = reader.query(loc.getContig(), loc.getStart(), loc.getStop(), false);
        }

        final GATKSAMRecordIterator iterator = new GATKSAMRecordIterator(rawIterator);

        final Set<String> samples = new HashSet<String>();
        for ( final SAMReadGroupRecord rg : reader.getFileHeader().getReadGroups() )
            samples.add(rg.getSample());

        final LIBSDownsamplingInfo ds = new LIBSDownsamplingInfo(downsample, 250);

        final LocusIteratorByState libs =
                new LocusIteratorByState(
                        iterator,
                        ds,
                        true,
                        genomeLocParser,
                        samples,
                        false);

        final SimpleTimer timer = new SimpleTimer().start();
        int bp = 0;
        double lastElapsed = 0;
        while ( libs.hasNext() ) {
            AlignmentContext context = libs.next();
            bp++;
            if ( timer.getElapsedTime() - lastElapsed > 10 ) {
                logger.info(bp + " iterations at " + context.getLocation());
                lastElapsed = timer.getElapsedTime();
            }
        }
        logger.info(String.format("runtime in seconds: %.2f", timer.getElapsedTime()));

        return 0;
    }

//    private void syntheticTests() {
//        final int readLength = 101;
//        final int nReads = 10000;
//        final int locus = 1;
//
//        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
//        final GenomeLocParser genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
//
//        int nIterations = 0;
//        for ( final String cigar : Arrays.asList("101M", "50M10I40M", "50M10D40M") ) {
//            GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "read", 0, locus, readLength);
//            read.setReadBases(Utils.dupBytes((byte) 'A', readLength));
//            final byte[] quals = new byte[readLength];
//            for ( int i = 0; i < readLength; i++ )
//                quals[i] = (byte)(i % QualityUtils.MAX_SAM_QUAL_SCORE);
//            read.setBaseQualities(quals);
//            read.setCigarString(cigar);
//
//            for ( int j = 0; j < nReads; j++ ) {
//                for ( int i = 0; i < rep; i++ ) {
//                    switch ( op ) {
//                        case NEW_STATE:
//                        {
//                            final AlignmentStateMachine alignmentStateMachine = new AlignmentStateMachine(read);
//                            while ( alignmentStateMachine.stepForwardOnGenome() != null ) {
//                                nIterations++;
//                            }
//                        }
//                        break;
////                        case OLD_STATE:
////                        {
////                            final SAMRecordAlignmentState alignmentStateMachine = new SAMRecordAlignmentState(read);
////                            while ( alignmentStateMachine.stepForwardOnGenome() != null ) {
////                                alignmentStateMachine.getRead();
////                                nIterations++;
////                            }
////                        }
////                        break;
//                        case NEW_LIBS:
//                        {
//                            final List<GATKSAMRecord> reads = Collections.nCopies(30, read);
//                            final org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState libs =
//                                    new org.broadinstitute.gatk.utils.locusiterator.LocusIteratorByState(
//                                            new LocusIteratorByStateBaseTest.FakeCloseableIterator<GATKSAMRecord>(reads.iterator()),
//                                            LocusIteratorByStateBaseTest.createTestReadProperties(),
//                                            genomeLocParser,
//                                            LocusIteratorByState.sampleListForSAMWithoutReadGroups());
//
//                            while ( libs.hasNext() ) {
//                                AlignmentContext context = libs.next();
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//        System.out.printf("iterations %d%n", nIterations);
//    }

    /**
     * Required main method implementation.
     * @param argv Command-line argument text.
     * @throws Exception on error.
     */
    public static void main(String[] argv) throws Exception {
        int returnCode = 0;
        try {
            LIBSPerformance instance = new LIBSPerformance();
            start(instance, argv);
            returnCode = 0;
        } catch(Exception ex) {
            returnCode = 1;
            ex.printStackTrace();
            throw ex;
        } finally {
            System.exit(returnCode);
        }
    }

}

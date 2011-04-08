package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.picard.sam.SamPairUtil;
import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.File;
import java.util.*;

//import org.broadinstitute.sting.utils.SimpleTimer;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 *
 * @author depristo
 * @version 0.1
 */
public class ReducingSAMFileWriter implements SAMFileWriter {
    protected static final Logger logger = Logger.getLogger(ReducingSAMFileWriter.class);
    private static final boolean DEBUG = false;
    private static final boolean INVERT = false;
    private static final boolean PRINT_CONSENSUS_READS = false;

    /** The place where we ultimately write out our records */
    final SAMFileWriter finalDestination;
    Queue<SAMRecord> waitingReads = new LinkedList<SAMRecord>();
    LinkedList<SAMRecord> consensusReads = new LinkedList<SAMRecord>();
    final int IndelContextSize;
    final int SNPContextSize;
    Set<VariantContext> VCs = new HashSet<VariantContext>();

    public long getNCompressedReads() {
        return nCompressedReads;
    }

    long nCompressedReads = 0;

    /**
     *
     * @param header
     * @param outputFile
     * @param compressionLevel
     */
    public ReducingSAMFileWriter(final SAMFileHeader header,
                                 final File outputFile,
                                 final int compressionLevel,
                                 final int SNPContextSize,
                                 final int IndelContextSize) {
        this(new SAMFileWriterFactory().makeBAMWriter(header, true, outputFile, compressionLevel),
                SNPContextSize,
                IndelContextSize);
    }

    public ReducingSAMFileWriter(final SAMFileWriter finalDestination,
                                 final int SNPContextSize,
                                 final int IndelContextSize) {
        this.finalDestination = finalDestination;
        this.IndelContextSize = IndelContextSize;
        this.SNPContextSize = SNPContextSize;
    }

    public int getIndelContextSize() {
        return IndelContextSize;
    }

    public int getSNPContextSize() {
        return SNPContextSize;
    }

    /**
     * Retrieves the header to use when creating the new SAM file.
     * @return header to use when creating the new SAM file.
     */
    public SAMFileHeader getFileHeader() {
        return finalDestination.getFileHeader();
    }

    /**
     * @{inheritDoc}
     */
    public void addAlignment( SAMRecord newRead ) {
        if ( DEBUG ) logger.info("New read pos " + newRead.getAlignmentStart() + " OP = " + newRead.getAttribute("OP"));

        // if the new read is on a different contig, then we need to flush the queue and clear the map
        if ( waitingReads.size() > 0 && waitingReads.peek().getReferenceIndex() != newRead.getReferenceIndex()) {
            if ( DEBUG ) logger.warn("Flushing queue on move to new contig: " + newRead.getReferenceName());

            // TODO -- fixme
            while ( ! waitingReads.isEmpty() ) {
                // emit to disk
                finalDestination.addAlignment(waitingReads.remove());
            }
        }

        waitingReads.add(newRead);
        emitReadsIfPossible(newRead.getAlignmentStart());
    }

    public void addVariant(VariantContext vc) {
        VCs.add(vc);
    }

    private void emitReadsIfPossible(int alignmentStartOfLastRead) {
        //
        // 2 states:
        //   -- reads ending << vc.getLocation()
        //   -- read overlapping vc.getLocation() +/- context size
        //
        while ( ! waitingReads.isEmpty() ) { // there's something in the queue
            SAMRecord read = waitingReads.peek();
            if ( ! withinContextOfVariants(read, VCs) ) {
                addToConsensus(read);
                waitingReads.remove();
            } else if ( cannotOverlapFutureVC(read, alignmentStartOfLastRead) ) {
                emitConsensus();
                if ( ! INVERT ) finalDestination.addAlignment(read);
                waitingReads.remove();
            } else {
                return;
            }
        }
    }

    private void addToConsensus(SAMRecord read) {
        if ( ! read.getDuplicateReadFlag() && ! read.getNotPrimaryAlignmentFlag() && ! read.getReadUnmappedFlag() )
            consensusReads.add(read);
    }

    private void emitConsensus() {
        if ( ! consensusReads.isEmpty() ) {
            SAMRecord firstRead = consensusReads.peek();

            int start = firstRead.getAlignmentStart();
            int end = furtherestEnd(consensusReads);
            int len = end - start + 1;

            int[][] baseCounts = new int[len][4];
            for ( SAMRecord read : consensusReads ) {
                int readI = 0, refI = read.getAlignmentStart() - start;
                for ( CigarElement elt : read.getCigar().getCigarElements() ) {
                    int l = elt.getLength();
                    switch (elt.getOperator()) {
                        case N: // cannot handle these
                            break;
                        case H : case P : // ignore pads and hard clips
                            break;
                        case S : refI += l; // move the reference too, in addition to I
                        case I :
                            // todo - handle insertions?
                            readI += l;
                            break;
                        case D :
                            refI += l;
                            break;
                        case M :
                            while (readI < l) {
                                byte base = read.getReadBases()[readI++];
                                int baseI = BaseUtils.simpleBaseToBaseIndex(base);
                                if ( baseI >= 0 ) // no Ns
                                    baseCounts[refI][baseI]++;
                                refI++;
                            }
                            break;
                        default:
                            throw new ReviewedStingException("BUG: Unexpected CIGAR element " + elt + " in read " + read.getReadName());
                    }
                }
            }

            byte[] bases = new byte[len];
            byte[] quals = new byte[len];
            for ( int i = 0; i < len; i++) {
                final int maxI = maxBaseCount(baseCounts[i]);
                final int count = baseCounts[i][maxI];
                final byte base = count == 0 ? (byte)'N' : BaseUtils.baseIndexToSimpleBase(maxI);
                bases[i] = base;
                quals[i] = QualityUtils.boundQual(count, (byte)64);
            }

            SAMRecord consensus = new SAMRecord(firstRead.getHeader());
            consensus.setReferenceIndex(firstRead.getReferenceIndex());
            consensus.setReadName("Mark");
            consensus.setCigarString(String.format("%dM", len));
            consensus.setReadPairedFlag(false);
            consensus.setAlignmentStart(start);
            consensus.setReadBases(bases);
            consensus.setBaseQualities(quals);
            consensus.setMappingQuality(60);

            int nConsensus = consensusReads.size();
            nCompressedReads += nConsensus;
            logger.info(String.format("Compressing %5d reads into a single consensus", nConsensus));
            finalDestination.addAlignment(consensus);

            if ( INVERT && PRINT_CONSENSUS_READS )
                for ( SAMRecord read : consensusReads )
                    finalDestination.addAlignment(read);

            consensusReads.clear();
        }
    }

    private static int furtherestEnd(Collection<SAMRecord> reads) {
        int end = -1;
        for ( SAMRecord read : reads ) {
            end = Math.max(end, read.getAlignmentEnd());
        }
        return end;
    }

    private int maxBaseCount(int[] counts) {
        int maxI = 0;
        for ( int i = 0; i < counts.length; i++) {
            if ( counts[i] > counts[maxI] ) {
                maxI = i;
            }
        }
        return maxI;
    }

    private boolean cannotOverlapFutureVC(SAMRecord read, int alignmentStartOfLastRead) {
        return read.getAlignmentEnd() < alignmentStartOfLastRead - Math.max(SNPContextSize, IndelContextSize);
    }

    private boolean withinContextOfVariants(SAMRecord read, Collection<VariantContext> vcs) {
        for ( VariantContext vc : vcs ) {
            if ( withinContextOfVariant(read, vc) ) {
                return true;
            }
        }

        return false;
    }

    private boolean withinContextOfVariant(SAMRecord read, VariantContext vc) {
        if ( ! read.getReferenceName().equals(vc.getChr()) )
            return false;
        else if ( vc.isVariant() ) {
            int contextSize = vc.isSNP() ? SNPContextSize : IndelContextSize;
            int vcContextStart = vc.getStart() - contextSize;
            int vcContextEnd = vc.getEnd() + contextSize;
            boolean notInContext = read.getAlignmentEnd() < vcContextStart || read.getAlignmentStart() > vcContextEnd;
            return ! notInContext;
        } else {
            return false;
        }
    }

    /**
     * @{inheritDoc}
     */
    public void close() {
        // write out all of the remaining reads
        while ( ! waitingReads.isEmpty() ) { // there's something in the queue
            finalDestination.addAlignment(waitingReads.remove());
        }
        finalDestination.close();
    }
}

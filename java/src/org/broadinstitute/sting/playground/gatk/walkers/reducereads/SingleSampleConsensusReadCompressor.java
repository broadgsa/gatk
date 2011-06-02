package org.broadinstitute.sting.oneoffprojects.walkers.reducereads;

import net.sf.samtools.*;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.clipreads.ClippingOp;
import org.broadinstitute.sting.utils.clipreads.ClippingRepresentation;
import org.broadinstitute.sting.utils.clipreads.ReadClipper;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import org.broadinstitute.sting.utils.sam.ReadUtils;

import java.io.PrintStream;
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
public class SingleSampleConsensusReadCompressor implements ConsensusReadCompressor {
    protected static final Logger logger = Logger.getLogger(SingleSampleConsensusReadCompressor.class);
    private static final boolean DEBUG = false;
    private static final boolean INVERT = false;
    private static final boolean PRINT_CONSENSUS_READS = false;
    private static final int CYCLES_BEFORE_RETRY = 1000;
    private static final double MAX_FRACTION_DISAGREEING_BASES = 0.1;
    private static final ClippingRepresentation VARIABLE_READ_REPRESENTATION = ClippingRepresentation.SOFTCLIP_BASES;
    private static final double MIN_FRACT_BASES_FOR_VARIABLE_READ = 0.33;  // todo -- should be variable
    private static final int MIN_BASES_IN_VARIABLE_SPAN_TO_INCLUDE_READ = 10;
    protected final static String RG_POSTFIX = ".ReducedReads";
    public final static int REDUCED_READ_BASE_QUALITY = 30;

    // todo  -- should merge close together spans

    /** The place where we ultimately write out our records */
    Queue<SAMRecord> waitingReads = new LinkedList<SAMRecord>();

    final int readContextSize;
    final int maxReadsAtVariableSites;
    final int minBpForRunningConsensus;
    int retryTimer = 0;
    int consensusCounter = 0;

    final SAMReadGroupRecord reducedReadGroup;
    String contig = null;
    final GenomeLocParser glParser;
    SAMFileHeader header;
    GenomeLoc lastProcessedRegion = null;

    public SingleSampleConsensusReadCompressor(final String sampleName,
                                               final int readContextSize,
                                               final GenomeLocParser glParser,
                                               final int minBpForRunningConsensus,
                                               final int maxReadsAtVariableSites) {
        this.readContextSize = readContextSize;
        this.glParser = glParser;
        this.minBpForRunningConsensus = minBpForRunningConsensus;
        this.maxReadsAtVariableSites = maxReadsAtVariableSites;
        this.reducedReadGroup = createReducedReadGroup(sampleName);
    }

    /**
     * Helper function to create a read group for these reduced reads
     * @param sampleName
     * @return
     */
    private static final SAMReadGroupRecord createReducedReadGroup(final String sampleName) {
        SAMReadGroupRecord rg = new SAMReadGroupRecord(sampleName + RG_POSTFIX);
        rg.setSample(sampleName);
        return rg;
    }

    public SAMReadGroupRecord getReducedReadGroup() {
        return reducedReadGroup;
    }

    // ------------------------------------------------------------------------------------------
    //
    // public interface functions
    //
    // ------------------------------------------------------------------------------------------

    /**
     * @{inheritDoc}
     */
    @Override
    public Iterable<SAMRecord> addAlignment( SAMRecord read ) {
        if ( contig == null )
            contig = read.getReferenceName();
        if ( ! read.getReferenceName().equals(contig) )
            throw new ReviewedStingException("ConsensusRead system doesn't support multiple contig processing right now");

        if ( header == null )
            header = read.getHeader();

        if ( ! waitingReads.isEmpty() && read.getAlignmentStart() < waitingReads.peek().getAlignmentStart() )
            throw new ReviewedStingException(
                    String.format("Adding read %s starting at %d before current queue head start position %d",
                            read.getReadName(), read.getAlignmentStart(), waitingReads.peek().getAlignmentStart()));

        Collection<SAMRecord> result = Collections.emptyList();
        if ( retryTimer == 0 ) {
            if ( chunkReadyForConsensus(read) ) {
                result = consensusReads(false);
            }
        } else {
            //logger.info("Retry: " + retryTimer);
            retryTimer--;
        }

        if ( ! read.getDuplicateReadFlag() && ! read.getNotPrimaryAlignmentFlag() && ! read.getReadUnmappedFlag() )
            waitingReads.add(read);

        return result;
    }

    @Override
    public Iterable<SAMRecord> close() {
        return consensusReads(true);
    }

    // ------------------------------------------------------------------------------------------
    //
    // private implementation functions
    //
    // ------------------------------------------------------------------------------------------

    private boolean chunkReadyForConsensus(SAMRecord read) {
        if ( ! waitingReads.isEmpty() ) {
            SAMRecord firstRead = waitingReads.iterator().next();
            int refStart = firstRead.getAlignmentStart();
            int refEnd = read.getAlignmentStart();
            int size = refEnd - refStart;
            return size > minBpForRunningConsensus;
        } else
            return false;
    }

//
//    public void writeConsensusBed(PrintStream bedOut) {
//        for ( ConsensusSite site : calculateConsensusSites(waitingReads) ) {
//            GenomeLoc loc = site.getLoc();
//            bedOut.printf("%s\t%d\t%d\t%s%n", loc.getContig(), loc.getStart()-1, loc.getStop(), site.counts);
//        }
//    }

    private Collection<SAMRecord> consensusReads(boolean useAllRemainingReads) {
        if ( ! waitingReads.isEmpty() ) {
            logger.info("Calculating consensus reads");
            List<ConsensusSite> sites = calculateConsensusSites(waitingReads, useAllRemainingReads, lastProcessedRegion);
            List<ConsensusSpan> rawSpans = calculateSpans(sites);
            List<ConsensusSpan> spans = useAllRemainingReads ? rawSpans : excludeFinalSpan(rawSpans);

            if ( ! spans.isEmpty() ) {
                lastProcessedRegion = spannedRegion(spans);
                logger.info("Processing region: " + lastProcessedRegion);
                updateWaitingReads(sites, spans);
                return consensusReadsFromSitesAndSpans(sites, spans);
            } else {
                logger.info("Danger, spans is empty, may experience poor performance at: " + spannedRegion(rawSpans));
                retryTimer = CYCLES_BEFORE_RETRY;
                return Collections.emptyList();
            }
        } else {
            return Collections.emptyList();
        }
    }

    private static final List<ConsensusSpan> excludeFinalSpan(List<ConsensusSpan> rawSpans) {
        logger.info("Dropping final, potentially incomplete span: " + rawSpans.get(rawSpans.size()-1));
        return rawSpans.subList(0, rawSpans.size() - 1);
    }

    private static final GenomeLoc spannedRegion(List<ConsensusSpan> spans) {
        GenomeLoc region = spans.get(0).loc;
        for ( ConsensusSpan span : spans )
            region = region.merge(span.loc);
        return region;
    }

    private void updateWaitingReads(List<ConsensusSite> sites, List<ConsensusSpan> spans) {
        ConsensusSpan lastSpan = spans.get(spans.size() - 1);
        Set<SAMRecord> unprocessedReads = new HashSet<SAMRecord>();

        for ( ConsensusSite site : sites.subList(lastSpan.getOffsetFromStartOfSites() + 1, sites.size()) ) {
            for ( PileupElement p : site.getOverlappingReads() )
                unprocessedReads.add(p.getRead());
        }

        logger.info(String.format("Updating waiting reads: old=%d reads, new=%d reads", waitingReads.size(), unprocessedReads.size()));
        waitingReads = new LinkedList<SAMRecord>(ReadUtils.coordinateSortReads(new ArrayList<SAMRecord>(unprocessedReads)));
    }

    private List<ConsensusSite> expandVariableSites(List<ConsensusSite> sites) {
        for ( ConsensusSite site : sites )
            site.setMarkedType(ConsensusSpan.Type.CONSERVED);

        for ( int i = 0; i < sites.size(); i++ ) {
            ConsensusSite site = sites.get(i);
            if ( ! site.isStrongConsensus(MAX_FRACTION_DISAGREEING_BASES) ) {
                int start = Math.max(i - readContextSize, 0);
                int stop = Math.min(sites.size(), i + readContextSize + 1);
                for ( int j = start; j < stop; j++ ) {
                    // aggressive tagging -- you are only conserved if you are never variable
                    sites.get(j).setMarkedType(ConsensusSpan.Type.VARIABLE);
                }
            }
        }

        return sites;
    }

    private List<ConsensusSpan> calculateSpans(List<ConsensusSite> rawSites) {
        List<ConsensusSite> sites = expandVariableSites(rawSites);
        List<ConsensusSpan> spans = new ArrayList<ConsensusSpan>();
        int start = 0;

        // our first span type is the type of the first site
        ConsensusSpan.Type consensusType = sites.get(0).getMarkedType();
        while ( start < sites.size() ) {
            ConsensusSpan span = findSpan(sites, start, consensusType);

            if ( span == null ) // we are done
                return spans;
            else {
                spans.add(span);
                start += span.size();
            }

            consensusType = ConsensusSpan.Type.otherType(consensusType);
        }

        return spans;
    }

    private ConsensusSpan findSpan(List<ConsensusSite> sites, int start, ConsensusSpan.Type consensusType) {
        int refStart = sites.get(0).getLoc().getStart();

        for ( int end = start + 1; end < sites.size(); end++ ) {
            ConsensusSite site = sites.get(end);
            boolean conserved = site.getMarkedType() == ConsensusSpan.Type.CONSERVED;
            if ( (consensusType == ConsensusSpan.Type.CONSERVED && ! conserved) ||
                 (consensusType == ConsensusSpan.Type.VARIABLE && conserved) ||
                 end + 1 == sites.size() ) { // we are done with the complete interval
                GenomeLoc loc = glParser.createGenomeLoc(contig, start+refStart, end+refStart-1);
                return new ConsensusSpan(refStart, loc, consensusType);
            }
        }

        return null; // couldn't find anything
    }


    private List<ConsensusSite> calculateConsensusSites(Collection<SAMRecord> reads, boolean useAllRemainingReads, GenomeLoc lastProcessedRegion) {
        SAMRecord firstRead = reads.iterator().next();

        int minStart = lastProcessedRegion == null ? -1 : lastProcessedRegion.getStop() + 1;
        int refStart = Math.max(firstRead.getAlignmentStart(), minStart);
        int refEnd = furtherestEnd(reads);

        logger.info("Calculating sites for region " + refStart + " to " + refEnd);
        // set up the consensus site array
        List<ConsensusSite> consensusSites = new ArrayList<ConsensusSite>();
        int len = refEnd - refStart + 1;
        for ( int i = 0; i < len; i++ ) {
            int l = refStart + i;
            GenomeLoc loc = glParser.createGenomeLoc(contig, l, l);
            consensusSites.add(new ConsensusSite(loc, i));
        }

        for ( SAMRecord read : reads ) {
            for ( RefPileupElement p : RefPileupElement.walkRead(read, refStart) ) {
                // add to the consensus at this site
                consensusSites.get(p.getRefOffset()).addOverlappingRead(p);
            }
        }

        return consensusSites;
    }

    private List<SAMRecord> consensusReadsFromSitesAndSpans(List<ConsensusSite> sites, List<ConsensusSpan> spans) {
        List<SAMRecord> reads = new ArrayList<SAMRecord>();

        for ( ConsensusSpan span : spans ) {
            logger.info("Span is " + span);
            if ( span.isConserved() )
                reads.addAll(conservedSpanReads(sites, span));
            else
                reads.addAll(downsample(variableSpanReads(sites, span)));
        }

        return reads;
    }

    // todo -- should be smart -- should take reads in some priority order
    // todo -- by length, and by strand, ideally.  Perhaps alternating by strand
    // todo -- in order of length?
    private Collection<SAMRecord> downsample(Collection<SAMRecord> reads) {
        if ( reads.size() > maxReadsAtVariableSites ) {
            List<SAMRecord> readArray = new ArrayList<SAMRecord>(reads);
            Collections.shuffle(readArray, GenomeAnalysisEngine.getRandomGenerator());
            return readArray.subList(0, maxReadsAtVariableSites);
        } else {
            return reads;
        }
    }

    private List<SAMRecord> conservedSpanReads(List<ConsensusSite> sites, ConsensusSpan span) {
        byte[] bases = new byte[span.size()];
        byte[] quals = new byte[span.size()];

        for ( int i = 0; i < span.size(); i++ ) {
            int refI = i + span.getOffsetFromStartOfSites();
            ConsensusSite site = sites.get(refI);
            if ( site.getMarkedType() == ConsensusSpan.Type.VARIABLE )
                throw new ReviewedStingException("Variable site included in consensus: " + site);
            final int count = site.counts.countOfMostCommonBase();
            final byte base = count == 0 ? (byte)'N' : site.counts.baseWithMostCounts();
            bases[i] = base;
            quals[i] = QualityUtils.boundQual(count, (byte)64);
        }

        SAMRecord consensus = new SAMRecord(header);
        consensus.setAttribute("RG", reducedReadGroup.getId());
        consensus.setAttribute(ReadUtils.REDUCED_READ_QUALITY_TAG, Integer.valueOf(REDUCED_READ_BASE_QUALITY));
        consensus.setReferenceName(contig);
        consensus.setReadName(String.format("%s.read.%d", reducedReadGroup.getId(), consensusCounter++));
        consensus.setReadPairedFlag(false);
        consensus.setReadUnmappedFlag(false);
        consensus.setCigarString(String.format("%dM", span.size()));
        consensus.setAlignmentStart(span.getGenomeStart());
        consensus.setReadBases(bases);
        consensus.setBaseQualities(quals);
        consensus.setMappingQuality(60);

//        if ( INVERT && PRINT_CONSENSUS_READS )
//            for ( SAMRecord read : consensusReads )
//                finalDestination.addAlignment(read);

        return Collections.singletonList(consensus);
    }

    private Collection<SAMRecord> variableSpanReads(List<ConsensusSite> sites, ConsensusSpan span) {
        Set<SAMRecord> reads = new HashSet<SAMRecord>();

        // todo -- this code is grossly inefficient, as it checks each variable read at each site in the span
        for ( int i = 0; i < span.size(); i++ ) {
            int refI = i + span.getOffsetFromStartOfSites();
            ConsensusSite site = sites.get(refI);
            for ( PileupElement p : site.getOverlappingReads() ) {
                SAMRecord read = clipReadToSpan(p.getRead(), span);
                if ( keepClippedReadInVariableSpan(p.getRead(), read) )
                    reads.add(read);
            }
        }

        return reads;
    }

    private final static boolean keepClippedReadInVariableSpan(SAMRecord originalRead, SAMRecord variableRead) {
        int originalReadLength = originalRead.getReadLength();
        int variableReadLength = variableRead.getReadLength();

        return variableReadLength >= MIN_BASES_IN_VARIABLE_SPAN_TO_INCLUDE_READ;
//        &&
//                ((1.0 * variableReadLength) / originalReadLength) >= MIN_FRACT_BASES_FOR_VARIABLE_READ;
    }

    private SAMRecord clipReadToSpan(SAMRecord read, ConsensusSpan span) {
        ReadClipper clipper = new ReadClipper(read);
        int spanStart = span.getGenomeStart();
        int spanEnd = span.getGenomeStop();
        int readLen = read.getReadLength();

        for ( RefPileupElement p : RefPileupElement.walkRead(read) ) {
            if ( p.getRefOffset() == spanStart && p.getOffset() != 0 ) {
                clipper.addOp(new ClippingOp(0, p.getOffset() - 1));
            }

            if ( p.getRefOffset() == spanEnd && p.getOffset() != readLen - 1 ) {
                clipper.addOp(new ClippingOp(p.getOffset() + 1, readLen - 1));
            }
        }

        SAMRecord softClipped = clipper.clipRead(VARIABLE_READ_REPRESENTATION);
        return ReadUtils.hardClipSoftClippedBases(softClipped);
    }

    private static int furtherestEnd(Collection<SAMRecord> reads) {
        int end = -1;
        for ( SAMRecord read : reads ) {
            end = Math.max(end, read.getAlignmentEnd());
        }
        return end;
    }
}

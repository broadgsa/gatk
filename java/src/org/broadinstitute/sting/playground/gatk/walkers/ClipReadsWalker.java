package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.BasicPileup;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;

import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.regex.Matcher;
import java.io.File;

import net.sf.samtools.util.StringUtil;

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
 * This walker prints out the reads from the BAM files provided to the traversal engines.
 * It also supports the command line option '-outputBamFile filname', which outputs all the
 * reads to a specified BAM file
 * The walker now also optionally filters reads based on command line options.
 */
@Requires({DataSource.READS})
public class ClipReadsWalker extends ReadWalker<ClipReadsWalker.ReadClipper, ClipReadsWalker.ClippingData> {
    /** an optional argument to dump the reads out to a BAM file */
    @Argument(fullName = "outputBam", shortName = "ob", doc = "Write output to this BAM filename instead of STDOUT", required = false)
    SAMFileWriter outputBam = null;

    @Argument(fullName = "", shortName = "STD", doc="FOR DEBUGGING ONLY", required = false)
    boolean toStandardOut = false;    

    @Argument(fullName = "qTrimmingThreshold", shortName = "QT", doc="", required = false)
    int qTrimmingThreshold = -1;

    @Argument(fullName = "cyclesToTrim", shortName = "CT", doc="String of the form 1-10,20-30 indicating machine cycles to clip from the reads", required = false)
    String cyclesToClipArg = null;

    @Argument(fullName = "clipSequencesFile", shortName = "XF", doc="Remove sequences within reads matching these sequences", required = false)
    String clipSequenceFile = null;

    @Argument(fullName = "clipSequence", shortName = "X", doc="Remove sequences within reads matching this sequence", required = false)
    String[] clipSequencesArgs = null;

//    @Argument(fullName = "onlyClipFirstSeqMatch", shortName = "ESC", doc="Only clip the first occurrence of a clipping sequence, rather than all subsequences within a read that match", required = false)
//    boolean onlyClipFirstSeqMatch = false;

    @Argument(fullName = "clipRepresentation", shortName = "CR", doc="How should we actually clip the bases?", required = false)
    ClippingRepresentation clippingRepresentation = ClippingRepresentation.WRITE_NS;

    /**
     * List of sequence that should be clipped from the reads
     */
    List<SeqToClip> sequencesToClip = new ArrayList<SeqToClip>();

    /**
     * List of cycle start / stop pairs (0-based, stop is included in the cycle to remove) to clip from the reads
     */
    List<Pair<Integer, Integer>> cyclesToClip = null;

    /**
     * The initialize function.
     */
    public void initialize() {
        if ( qTrimmingThreshold >= 0 ) {
            logger.info(String.format("Creating Q-score clipper with threshold %d", qTrimmingThreshold));
        }

        //
        // Initialize the sequences to clip
        //
        if ( clipSequencesArgs != null ) {
            int i = 0;
            for ( String toClip : clipSequencesArgs ) {
                i++;
                ReferenceSequence rs = new ReferenceSequence("CMDLINE-"+i, -1, StringUtil.stringToBytes(toClip));
                addSeqToClip(rs.getName(), rs.getBases());
            }
        }

        if ( clipSequenceFile != null ) {
            ReferenceSequenceFile rsf = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(clipSequenceFile));
            
            while ( true ) {
                ReferenceSequence rs = rsf.nextSequence();
                if ( rs == null )
                    break;
                else {
                    addSeqToClip(rs.getName(), rs.getBases());
                }
            }
        }



        //
        // Initialize the cycle ranges to clip
        //
        if ( cyclesToClipArg != null ) {
            cyclesToClip = new ArrayList<Pair<Integer, Integer>>();
            for ( String range : cyclesToClipArg.split(",") ) {
                try {
                    String[] elts = range.split("-");
                    int start = Integer.parseInt(elts[0]) - 1;
                    int stop = Integer.parseInt(elts[1]) - 1;

                    if ( start < 0 ) throw new Exception();
                    if ( stop < start ) throw new Exception();

                    logger.info(String.format("Creating cycle clipper %d-%d", start, stop));
                    cyclesToClip.add(new Pair<Integer, Integer>(start, stop));
                } catch ( Exception e ) {
                    throw new RuntimeException("Badly formatted cyclesToClip argument: " + cyclesToClipArg);
                }
            }
        }
    }

    /**
     * Helper function that adds a seq with name and bases (as bytes) to the list of sequences to be clipped
     * 
     * @param name
     * @param bases
     */
    private void addSeqToClip(String name, byte[] bases) {
        SeqToClip clip = new SeqToClip(name, StringUtil.bytesToString(bases));
        sequencesToClip.add(clip);
        logger.info(String.format("Creating sequence clipper %s: %s/%s", clip.name, clip.seq, clip.revSeq));
    }

    /**
     * The reads map function.
     * @param ref the reference bases that correspond to our read, if a reference was provided
     * @param read the read itself, as a SAMRecord
     * @return the ReadClipper object describing what should be done to clip this read
     */
    public ReadClipper map( char[] ref, SAMRecord read ) {
        ReadClipper clipper = new ReadClipper(read);

        //
        // run all three clipping modules
        //
        clipBadQualityScores(clipper);
        clipCycles(clipper);
        clipSequences(clipper);

        return clipper;
    }

    /**
     * clip sequences from the reads that match all of the sequences in the global sequencesToClip variable.
     * Adds ClippingOps for each clip to clipper.
     *
     * @param clipper
     */
    private void clipSequences(ReadClipper clipper) {
        if ( sequencesToClip != null ) {                // don't bother if we don't have any sequences to clip
            SAMRecord read = clipper.getRead();

            for ( SeqToClip stc : sequencesToClip ) {
                // we have a pattern for both the forward and the reverse strands
                Pattern pattern = read.getReadNegativeStrandFlag() ? stc.revPat : stc.fwdPat;
                String bases = read.getReadString();
                Matcher match = pattern.matcher(bases);

                // keep clipping until match.find() says it can't find anything else
                boolean found = true;   // go through at least once
                while ( found ) {
                    found = match.find();
                    //System.out.printf("Matching %s against %s/%s => %b%n", bases, stc.seq, stc.revSeq, found);
                    if ( found ) {
                        int start = match.start();
                        int stop = match.end() - 1;
                        ClippingOp op = new ClippingOp(ClippingType.MATCHES_CLIP_SEQ, start, stop, stc.seq);
                        clipper.addOp(op);
                    }
                }
            }
        }
    }

    /**
     * Convenence function that takes a read and the start / stop clipping positions based on the forward
     * strand, and returns start/stop values appropriate for the strand of the read.
     *
     * @param read
     * @param start
     * @param stop
     * @return
     */
    private Pair<Integer, Integer> strandAwarePositions( SAMRecord read, int start, int stop ) {
        if ( read.getReadNegativeStrandFlag() )
            return new Pair<Integer, Integer>(  read.getReadLength() - stop - 1, read.getReadLength() - start - 1 );
        else
            return new Pair<Integer, Integer>(start, stop);
    }

    /**
     * clip bases at cycles between the ranges in cyclesToClip by adding appropriate ClippingOps to clipper.
     *
     * @param clipper
     */
    private void clipCycles(ReadClipper clipper) {
        if ( cyclesToClip != null ) {
            SAMRecord read = clipper.getRead();

            for ( Pair<Integer, Integer> p : cyclesToClip ) {   // iterate over each cycle range
                int cycleStart = p.first;
                int cycleStop = p.second;

                if ( cycleStart < read.getReadLength() ) {
                    // only try to clip if the cycleStart is less than the read's length
                    if ( cycleStop >= read.getReadLength() )
                        // we do tolerate [for convenience) clipping when the stop is beyond the end of the read
                        cycleStop = read.getReadLength() - 1;

                    Pair<Integer, Integer> startStop = strandAwarePositions( read, cycleStart, cycleStop );
                    int start = startStop.first;
                    int stop = startStop.second;

                    ClippingOp op = new ClippingOp(ClippingType.WITHIN_CLIP_RANGE, start, stop, null);
                    clipper.addOp(op);
                }
            }
        }
    }

    /**
     * Clip bases from the read in clipper from
     *
     * argmax_x{ \sum{i = x + 1}^l (qTrimmingThreshold - qual)
     *
     * to the end of the read.  This is blatantly stolen from BWA.
     *
     * Walk through the read from the end (in machine cycle order) to the beginning, calculating the
     * running sum of qTrimmingThreshold - qual.  While we do this, we track the maximum value of this
     * sum where the delta > 0.  After the loop, clipPoint is either -1 (don't do anything) or the
     * clipping index in the read (from the end).
     *
     * @param clipper
     */
    private void clipBadQualityScores(ReadClipper clipper) {
        SAMRecord read = clipper.getRead();
        int readLen = read.getReadBases().length;
        byte[] quals = read.getBaseQualities();


        int clipSum = 0, lastMax = -1, clipPoint = -1; // -1 means no clip
        for ( int i = readLen - 1; i >= 0; i-- ) {
            int baseIndex = read.getReadNegativeStrandFlag() ? readLen - i - 1 : i;
            byte qual = quals[baseIndex];
            clipSum += (qTrimmingThreshold - qual);
            if ( clipSum >= 0 && ( clipSum >= lastMax ) ) {
                lastMax = clipSum;
                clipPoint = baseIndex;
            }
        }

        if ( clipPoint != -1 ) {
            int start = read.getReadNegativeStrandFlag() ? 0 : clipPoint;
            int stop = read.getReadNegativeStrandFlag() ? clipPoint : readLen - 1;
            clipper.addOp(new ClippingOp(ClippingType.LOW_Q_SCORES, start, stop, null));
        }
    }

    /**
     * reduceInit is called once before any calls to the map function.  We use it here to setup the output
     * bam file, if it was specified on the command line
     * @return 
     */
    public ClippingData reduceInit() {
        return new ClippingData(outputBam, sequencesToClip);
    }

    public ClippingData reduce( ReadClipper clipper, ClippingData data ) {
        if (data.output != null) {
            data.output.addAlignment(clipper.clipRead(clippingRepresentation));
        } else if ( toStandardOut ) {
            out.println(clipper.clipRead(clippingRepresentation).format());
        }

        data.nTotalReads++;
        data.nTotalBases += clipper.getRead().getReadLength();
        if ( clipper.wasClipped() ) {
            data.nClippedReads++;
            for ( ClippingOp op : clipper.getOps() ) {
                switch ( op.type ) {
                    case LOW_Q_SCORES:
                        data.incNQClippedBases( op.getLength() );
                        break;
                    case WITHIN_CLIP_RANGE:
                        data.incNRangeClippedBases( op.getLength() );
                        break;
                    case MATCHES_CLIP_SEQ:
                        data.incSeqClippedBases( (String)op.extraInfo, op.getLength() );
                        break;
                }
            }
        }

        return data;
    }

    public void onTraversalDone( ClippingData data ) {
        out.printf(data.toString());
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // utility classes
    //
    // --------------------------------------------------------------------------------------------------------------

    private class SeqToClip {
        String name;
        String seq, revSeq;
        Pattern fwdPat, revPat;

        public SeqToClip(String name, String seq) {
            this.name = name;
            this.seq = seq;
            this.fwdPat = Pattern.compile(seq, Pattern.CASE_INSENSITIVE);
            this.revSeq = BaseUtils.simpleReverseComplement(seq);
            this.revPat = Pattern.compile(revSeq, Pattern.CASE_INSENSITIVE);
        }
    }

    /**
     * What is the type of a ClippingOp?
     */
    private enum ClippingType {
        LOW_Q_SCORES,
        WITHIN_CLIP_RANGE,
        MATCHES_CLIP_SEQ
    }

    /**
     * Represents a clip on a read.  It has a type (see the enum) along with a start and stop in the bases
     * of the read, plus an option extraInfo (useful for carrying info where needed).
     *
     * Also holds the critical apply function that actually execute the clipping operation on a provided read,
     * according to the wishes of the supplid ClippingAlgorithm enum.
     */
    private class ClippingOp {
        public ClippingType type;
        public int start, stop; // inclusive
        public Object extraInfo = null;

        public ClippingOp(ClippingType type, int start, int stop, Object extraInfo ) {
            this.type = type;
            this.start = start;
            this.stop = stop;
            this.extraInfo = extraInfo;
        }

        public int getLength() { return stop - start + 1; }

        /**
         * Clips the bases in clippedRead according to this operation's start and stop.  Uses the clipping
         * representation used is the one provided by algorithm argument.
         *
         * @param algorithm
         * @param clippedRead
         */
        public void apply(ClippingRepresentation algorithm, SAMRecord clippedRead) {
            switch ( algorithm ) {
                case WRITE_NS:
                    for ( int i = start; i <= stop; i++ )
                        clippedRead.getReadBases()[i] = 'N';
                    break;
                case WRITE_Q0S:
                    for ( int i = start; i <= stop; i++ )
                        clippedRead.getBaseQualities()[i] = 0;
                    break;
                case SOFTCLIP_BASES:
                    throw new RuntimeException("Softclipping of bases not yet implemented.");
            }
        }
    }

    /**
     * How should we represent a clipped bases in a read?
     */
    private enum ClippingRepresentation {
        WRITE_NS,           // change the bases to Ns
        WRITE_Q0S,          // change the quality scores to Q0
        SOFTCLIP_BASES      // change cigar string to S, but keep bases
    }

    /**
     * A simple collection of the clipping operations to apply to a read along with its read
     *
     */
    public class ReadClipper {
        SAMRecord read;
        List<ClippingOp> ops = null;

        /**
         * We didn't do any clipping work on this read, just leave everything as a default
         * @param read
         */
        public ReadClipper(final SAMRecord read) {
            this.read = read;
        }

        /**
         * Add another clipping operation to apply to this read
         * @param op
         */
        public void addOp( ClippingOp op ) {
            if ( ops == null ) ops = new ArrayList<ClippingOp>();
            ops.add(op);
        }

        public List<ClippingOp> getOps() { return ops; }
        public boolean wasClipped() { return ops != null; }
        public SAMRecord getRead() { return read; }


        /**
         * Return a new read corresponding to this.read that's been clipped according to ops, if any are present.
         *
         * @param algorithm
         * @return
         */
        public SAMRecord clipRead(ClippingRepresentation algorithm) {
            if ( ops == null )
                return getRead();
            else {
                try {
                    SAMRecord clippedRead = (SAMRecord)read.clone();
                    for ( ClippingOp op : getOps() ) {
                        op.apply(algorithm, clippedRead);
                    }
                    return clippedRead;
                } catch (CloneNotSupportedException e) {
                    throw new RuntimeException(e); // this should never happen
                }
            }
        }
    }

    public class ClippingData {
        public SAMFileWriter output = null;

        public long nTotalReads = 0;
        public long nTotalBases = 0;
        public long nClippedReads = 0;
        public long nClippedBases = 0;
        public long nQClippedBases = 0;
        public long nRangeClippedBases = 0;
        public long nSeqClippedBases = 0;

        HashMap<String, Long> seqClipCounts = new HashMap<String, Long>();

        public ClippingData(SAMFileWriter output, List<SeqToClip> clipSeqs) {
            this.output = output;
            for ( SeqToClip clipSeq : clipSeqs ) {
                seqClipCounts.put(clipSeq.seq, 0L);
            }
        }

        public void incNQClippedBases( int n ) {
            nQClippedBases += n;
            nClippedBases += n;
        }

        public void incNRangeClippedBases( int n ) {
            nRangeClippedBases += n;
            nClippedBases += n;
        }

        public void incSeqClippedBases( final String seq, int n ) {
            nSeqClippedBases += n;
            nClippedBases += n;
            seqClipCounts.put(seq, seqClipCounts.get(seq) + n);
        }

        public String toString() {
            StringBuilder s = new StringBuilder();

            s.append(Utils.dupString('-', 80) + "\n");
            s.append(String.format("Number of examined reads              %d%n", nTotalReads));
            s.append(String.format("Number of clipped reads               %d%n", nClippedReads));
            s.append(String.format("Percent of clipped reads              %.2f%n", (100.0*nClippedReads) / nTotalReads));
            s.append(String.format("Number of examined bases              %d%n", nTotalBases));
            s.append(String.format("Number of clipped bases               %d%n", nClippedBases));
            s.append(String.format("Percent of clipped bases              %.2f%n", (100.0*nClippedBases) / nTotalBases));
            s.append(String.format("Number of quality-score clipped bases %d%n", nQClippedBases));
            s.append(String.format("Number of range clipped bases         %d%n", nRangeClippedBases));
            s.append(String.format("Number of sequence clipped bases      %d%n", nSeqClippedBases));

            for ( Map.Entry<String, Long> elt : seqClipCounts.entrySet() ) {
                s.append(String.format("  %8d clip sites matching %s%n", elt.getValue(), elt.getKey() ));
            }

            s.append(Utils.dupString('-', 80) + "\n");
            return s.toString();
        }
    }
}
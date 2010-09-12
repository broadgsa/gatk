/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.sam.AlignmentUtils;
import net.sf.samtools.SAMRecord;

import java.util.*;
import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Jan 3, 2010
 * Time: 1:58:38 PM
 * To change this template use File | Settings | File Templates.
 */
public class DSBWalkerV3 extends ReadWalker<Integer,Integer> {
    @Output
    PrintStream out;

    @Argument(fullName="windowSize",shortName="W",doc="Size of the sliding window",required=true)
    int WINDOW_SIZE = 100;
    @Argument(fullName="enrichmentCutoff",shortName="E",doc="Report windows with enrichment (signal/control) above this cutoff",required=true)
    double ENRICHMENT_CUTOFF = 5.0;
    @Argument(fullName="minSignal",shortName="ms",doc="Do not report windows with signal lower than this value "+
                "(this cutoff is secondary to enrichmentCutoff and guards against windows where control signal is 0 or too low,"+
                "so that control*enrichmentCutoff is too low to be convincing)",required=true)
    int MIN_SIGNAL = 10;
    @Argument(fullName="coverageFactor",shortName="cf",doc="Total number of uniquely mapped signal reads/total number of uniquely mapped control reads",required=false)
    double COVERAGE_FACTOR=1.0;
    @Argument(fullName="coverageFactorNU",shortName="cfnu",doc="Total number of non-uniquely mapped signal reads/total number of non-uniquely mapped control reads",required=false)
    double COVERAGE_FACTOR_NU=1.0;

    private Set<String> signalReadGroups; // we are going to remember which read groups are stimulated tagged and which are unstimulated untagged in order to be able
    private Set<String> controlReadGroups ; // to properly assign the reads coming from a merged stream

    private GenomeLoc currentWindow = null;
    private String currentContig = "chrM";


    private LinkedList<SAMRecord> readsInSignalWindow = null;
    private LinkedList<SAMRecord> readsInControlWindow = null;

    private WindowStats signalCountsInCurrWindow = new WindowStats();
    private WindowStats controlCountsInCurrWindow = new WindowStats();


    // following variables are used by emitWindow to buffer adjacent windows
    private int MERGE_CUTOFF = -1;

    private long regionStart = -1;
    private long lastWindowStart = -1;
    private int addedSinceLastEmit = 0; // how many sliding window steps where buffered since the last emit (i.e. since the last window that really passed significance criteria)
    // buffered read count stats for the windows inside the currently held merged print region:
    private List<WindowStats> signalReadCountsBuffer = new ArrayList<WindowStats>(1000);
    private List<WindowStats> controlReadCountsBuffer = new ArrayList<WindowStats>(1000);


    /** Clears buffered reads and all counts. DOES NOT clear buffered print region */
    private void resetWindows() {
        readsInSignalWindow.clear();
        readsInControlWindow.clear();
        signalCountsInCurrWindow.clear();
        controlCountsInCurrWindow.clear();
    }

    private void addSignal(SAMRecord read) {
        readsInSignalWindow.add(read);
        signalCountsInCurrWindow.addRead(read);
    }

    private void addControl(SAMRecord read) {
        readsInControlWindow.add(read);
        controlCountsInCurrWindow.addRead(read);
    }

    /** Discard signal reads that start strictly before the specified position and
     * update associated counts
     * @param pos
     */
    private void purgeSignal(long pos) {
        Iterator<SAMRecord> it = readsInSignalWindow.iterator();
        while ( it.hasNext() ) {
            SAMRecord r = it.next();
            if ( r.getAlignmentStart() >= pos ) return; // we are done

            // read starts before pos: discard it and update the counts:
            signalCountsInCurrWindow.removeRead(r);
            it.remove();
        }
    }

    /** Discard signal reads that start strictly before the specified position and
     * update associated counts
     * @param pos
     */
    private void purgeControl(long pos) {
        Iterator<SAMRecord> it = readsInControlWindow.iterator();
        while ( it.hasNext() ) {
            SAMRecord r = it.next();
            if ( r.getAlignmentStart() >= pos ) return; // we are done

            // read starts before pos: discard it and update the counts:
            controlCountsInCurrWindow.removeRead(r);
            it.remove();
        }
    }

    private void resetWindowMergingBuffer(long start) {
        regionStart = start;
        lastWindowStart = start;
        signalReadCountsBuffer.clear();
        controlReadCountsBuffer.clear();
        signalReadCountsBuffer.add(signalCountsInCurrWindow.clone());
        controlReadCountsBuffer.add(controlCountsInCurrWindow.clone());
    }

    /** Delayed print: the window starting at 'start' will be added to the print buffer; if the window is close enough
     * to the current contents if the buffer, the addition will result in merging the window with the buffer;
     * otherwise, the old contents of the buffer will be printed and the buffer will be re-initialized with new window.
     * It is assumed that counters are in synch with the start position passed to this method.
     * @param start
     */
    private void emitWindow(long start) {
        //        System.out.println("Emitting at "+start);

        if ( regionStart == -1 ) { // we did not keep any region so far; initialize the buffer and return, will print later
            resetWindowMergingBuffer(start);
            addedSinceLastEmit = 0;
            return;
        }

        if ( start > lastWindowStart + MERGE_CUTOFF ) {
            // this loop is a dummy: we have already cleared those unneeded
            // counts in shiftWindows(); stays here to avoid generating bugs later
            // if we change something in shiftWindows()
            for ( ; addedSinceLastEmit > 0 ; addedSinceLastEmit-- ) {
                signalReadCountsBuffer.remove(signalReadCountsBuffer.size()-1);
                controlReadCountsBuffer.remove(controlReadCountsBuffer.size()-1);
            }
            printRegion();
            resetWindowMergingBuffer(start);
            return;
        }

        // the current window is too close to the previous one: we have to merge;
        // NOTE: if window is too close, bufferAccepts() returned true, so the counts are already
        // added.
        lastWindowStart = start;
        addedSinceLastEmit = 0;
//        signalReadCountsBuffer.add(uniqueSignalReads);
//        controlReadCountsBuffer.add(uniqueControlReads);

    }

    private boolean bufferAccepts(long pos) {
        return ( regionStart != -1 && pos <= lastWindowStart+MERGE_CUTOFF);
    }

   

    private void printRegion() {
        if ( regionStart == -1 ) return;

        long regionStop = lastWindowStart+WINDOW_SIZE-1;

        double[] tmpEnrU = new double[signalReadCountsBuffer.size()];
        int[] tmpSignalU = new int[signalReadCountsBuffer.size()];
        int[] tmpControlU = new int[signalReadCountsBuffer.size()];
        double[] tmpEnrNU = new double[signalReadCountsBuffer.size()];
        int[] tmpSignalNU = new int[signalReadCountsBuffer.size()];
        int[] tmpControlNU = new int[signalReadCountsBuffer.size()];

        double[] tmpFWDSignalFracU = new double[signalReadCountsBuffer.size()];
        double[] tmpFWDControlFracU = new double[signalReadCountsBuffer.size()];
        double[] tmpFWDSignalFracNU = new double[signalReadCountsBuffer.size()];
        double[] tmpFWDControlFracNU = new double[signalReadCountsBuffer.size()];

        int lastInd = signalReadCountsBuffer.size() - 1;

        //        out.println("Size="+signalReadCountsBuffer.size()+":");
 
        for ( int i = 0 ; i <= lastInd  ; i++ ) {
            tmpEnrU[i]= ( ((double) signalReadCountsBuffer.get(i).uniqueReads) / (controlReadCountsBuffer.get(i).uniqueReads+1.0 ) ) / COVERAGE_FACTOR ;

            tmpSignalU[i] = signalReadCountsBuffer.get(i).uniqueReads;
            tmpControlU[i] = controlReadCountsBuffer.get(i).uniqueReads;

            tmpEnrNU[i]= ( ((double) signalReadCountsBuffer.get(i).nonUniqueReads) / (controlReadCountsBuffer.get(i).nonUniqueReads+1.0 ) ) / COVERAGE_FACTOR_NU ;

            tmpSignalNU[i] = signalReadCountsBuffer.get(i).nonUniqueReads;
            tmpControlNU[i] = controlReadCountsBuffer.get(i).nonUniqueReads;

            tmpFWDSignalFracU[i] = signalReadCountsBuffer.get(i).uniqueReads > 0 ? ( ((double)signalReadCountsBuffer.get(i).uniqueFWDReads) / signalReadCountsBuffer.get(i).uniqueReads ) : 0.5;
            tmpFWDControlFracU[i] = controlReadCountsBuffer.get(i).uniqueReads > 0 ? ( ((double)controlReadCountsBuffer.get(i).uniqueFWDReads) / controlReadCountsBuffer.get(i).uniqueReads ) : 0.5;
            tmpFWDSignalFracNU[i] = signalReadCountsBuffer.get(i).nonUniqueReads > 0 ? ( ((double)signalReadCountsBuffer.get(i).nonUniqueFWDReads) / signalReadCountsBuffer.get(i).nonUniqueReads ) : 0.5;
            tmpFWDControlFracNU[i] = controlReadCountsBuffer.get(i).nonUniqueReads > 0 ? ( ((double)controlReadCountsBuffer.get(i).nonUniqueFWDReads) / controlReadCountsBuffer.get(i).nonUniqueReads ) : 0.5;
        }

        Arrays.sort(tmpEnrU);
        Arrays.sort(tmpSignalU);
        Arrays.sort(tmpControlU);

        Arrays.sort(tmpEnrNU);
        Arrays.sort(tmpSignalNU);
        Arrays.sort(tmpControlNU);

        Arrays.sort(tmpFWDSignalFracU);
        Arrays.sort(tmpFWDControlFracU);
        Arrays.sort(tmpFWDSignalFracNU);
        Arrays.sort(tmpFWDControlFracNU);


        out.print(currentContig+":"+regionStart+"-"+regionStop+"\t"+
                (regionStop-regionStart+1) +"\t"+
                "signal_unique:"+ tmpSignalU[0]+"-"+ tmpSignalU[lastInd/2]+"-"+ tmpSignalU[lastInd]+"\t"+
                  "control_unique:"+ tmpControlU[0]+"-"+ tmpControlU[lastInd/2]+"-"+ tmpControlU[lastInd]);

        out.printf("\tsignal_fwd_frac_unique:%.1f-%.1f-%.1f",tmpFWDSignalFracU[0],tmpFWDSignalFracU[lastInd/2],tmpFWDSignalFracU[lastInd]);
        out.printf("\tcontrol_fwd_frac_unique:%.1f-%.1f-%.1f",tmpFWDControlFracU[0],tmpFWDControlFracU[lastInd/2],tmpFWDControlFracU[lastInd]);

        out.print("\tsignal_nonnunique:"+ tmpSignalNU[0]+"-"+ tmpSignalNU[lastInd/2]+"-"+ tmpSignalNU[lastInd]+"\t"+
                  "control_nonunique:"+ tmpControlNU[0]+"-"+ tmpControlNU[lastInd/2]+"-"+ tmpControlNU[lastInd]);

        out.printf("\tsignal_fwd_frac_nonunique:%.1f-%.1f-%.1f",tmpFWDSignalFracNU[0],tmpFWDSignalFracNU[lastInd/2],tmpFWDSignalFracNU[lastInd]);
        out.printf("\tcontrol_fwd_frac_nonunique:%.1f-%.1f-%.1f",tmpFWDControlFracNU[0],tmpFWDControlFracNU[lastInd/2],tmpFWDControlFracNU[lastInd]);

        out.printf("\tnorm_enrichment_unique:%.2f-%.2f-%.2f",tmpEnrU[0],tmpEnrU[lastInd/2],tmpEnrU[lastInd]);
        out.printf("\tnorm_enrichment_nonunique:%.2f-%.2f-%.2f",tmpEnrNU[0],tmpEnrNU[lastInd/2],tmpEnrNU[lastInd]);

  //      if ( minUniqueSignalStrandBalance > 0.75 || minUniqueSignalStrandBalance < 0.25 ) out.print("\tS_U_STRAND_FILTER");
        out.println();

        regionStart = -1; // to indicate that there is nothing left to print, the buffer is empty
        //      System.exit(1);
    }

    private void updateWindowMergingBuffer(long i) {
        if ( bufferAccepts(i) ) {
            // we are not too far away from last window added to the buffer that actually passed significance criteria;
            // in this case we have to keep buffering since another significant window may be encountered soon
            //            System.out.println("Updating buffer at "+i+" with "+ uniqueSignalReads);
            signalReadCountsBuffer.add(signalCountsInCurrWindow.clone());
            controlReadCountsBuffer.add(controlCountsInCurrWindow.clone());
            addedSinceLastEmit++;
        } else {
            // we are too far from the last significant window; if another significant window comes later, it will not
            // be merged into this region but will start a new one. In this case we have to erase all the counts we have been
            // saving since the last significant window (the latter is where the current region is going to end!)
            for ( ; addedSinceLastEmit > 0 ; addedSinceLastEmit-- ) {
                signalReadCountsBuffer.remove(signalReadCountsBuffer.size()-1);
                controlReadCountsBuffer.remove(controlReadCountsBuffer.size()-1);
            }
            printRegion(); // print current region right away, why not? next significant window will start new region for sure.
        }

    }

    private void shiftWindows(long pos) {
        // we shift windows when there is a read that does not fit into the current window.
        // the position, to which the shift is performed, is the first position such that the new read
        // can be accomodated. Hence we can safely slide up to pos, only discarding reads that go out of scope -
        // we are guaranteed that there will be no new reads to add until we reach pos.


        for ( long i = currentWindow.getStart() ; i < pos ; i++ ) {
//            if ( readsInSignalWindow.size() == 0 ) {
//                i = pos-1;
//                continue;
//            };

//            if ( readsInSignalWindow.getFirst().getAlignmentStart() > i ) {
//                i = readsInSignalWindow.getFirst().getAlignmentStart() - 1; // jump directly to next read position
//                continue;
//            }

            purgeSignal(i);        // remove all the reads that start before current position i (and update all the counters)
            purgeControl(i);

            updateWindowMergingBuffer(i);

            if ( ( controlCountsInCurrWindow.uniqueReads + 1 ) * ENRICHMENT_CUTOFF < MIN_SIGNAL ) {
                 // too few control reads
                if ( signalCountsInCurrWindow.uniqueReads >= MIN_SIGNAL ) {
                    // emit signal only if it is higher that hard cut-off:
                    emitWindow(i); // print current window (print can be buffered and delayed!)
                }
            } else {
                // enough control reads;
                // check for actual enrichment:
                if ( ((double) signalCountsInCurrWindow.uniqueReads) / (controlCountsInCurrWindow.uniqueReads+1.0) > ENRICHMENT_CUTOFF ) {
                    emitWindow(i); // print current window (print can be buffered and delayed!)
                }
            }

        }

        // we emitted intermediate windows up to pos-1 as/if needed and purged everything that starts before pos-1
        // now we have to purge everything that starts before pos and return (no emitting yet, as we are about to add a read upon return):

        purgeSignal(pos);
        purgeControl(pos);

        currentWindow = GenomeLocParser.createGenomeLoc(currentWindow.getContigIndex(),pos,pos+WINDOW_SIZE-1);
    }

    @Override
    public void initialize() {
        int nSams = getToolkit().getArguments().samFiles.size();

        if ( nSams != 2 ) {
            out.println("ERROR: two input bam files (signal and backround control) must be specified");
            System.exit(1);
        }
        List<Set<String>> readGroupSets = getToolkit().getMergedReadGroupsByReaders();
        signalReadGroups = readGroupSets.get(0);
//        System.out.println(signalReadGroups.size()+" read groups in signal");
        controlReadGroups = readGroupSets.get(1);
//        System.out.println(controlReadGroups.size()+" read groups in control");

        currentWindow = GenomeLocParser.createGenomeLoc(0,1,WINDOW_SIZE);
        readsInSignalWindow = new LinkedList<SAMRecord>();
        readsInControlWindow = new LinkedList<SAMRecord>();

        MERGE_CUTOFF = WINDOW_SIZE;
        ENRICHMENT_CUTOFF *= COVERAGE_FACTOR;
        currentContig = getToolkit().getSAMFileHeader().getSequenceDictionary().getSequence(0).getSequenceName();
    }


    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {

        if ( AlignmentUtils.isReadUnmapped(read) ) return 0;

        if ( read.getReferenceIndex() > currentWindow.getContigIndex() ) {
            printRegion(); // print all we had on the previous contig

            currentWindow = GenomeLocParser.createGenomeLoc(read.getReferenceIndex(),
                    read.getAlignmentStart(),
                    read.getAlignmentStart()+WINDOW_SIZE-1);
            currentContig = read.getReferenceName();
            resetWindows();
        } else {
            // we are on the same contig
            if ( read.getAlignmentEnd() > currentWindow.getStop() ) {
                // can not accomodate the read inside the current window - shift!
                //                System.out.println("read ends at "+read.getAlignmentEnd()+" window ends at "+currentWindow.getStop()+ " shifting to "+ (currentWindow.getStart() + ( read.getAlignmentEnd() - currentWindow.getStop() )) +" ("+uniqueSignalReads+"/"+uniqueControlReads+")");

                // while shifting the window, the following method will issue (delayed) print commands for
                // all intermediate windows that pass significance criteria:
                shiftWindows(currentWindow.getStart() + ( read.getAlignmentEnd() - currentWindow.getStop() ));
            }
            // now the read will fit into the window
        }

        // at this point we are guaranteed that the read will fit into the window

        if ( signalReadGroups.contains( read.getReadGroup().getReadGroupId() ) ) {
            addSignal(read);
        } else if ( controlReadGroups.contains( read.getReadGroup().getReadGroupId() )) {
            addControl(read);
        } else {
            throw new UserException.MalformedBam(read, "Read "+read + " belongs to unrecognized read group");
        }
        return 1;
    }


    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return value+sum;  //To change body of implemented methods use File | Settings | File Templates.
    }


    /** Auxiliary class that encapsulates the task of monitoring counts of various read traits in some set of reads
     *  (for instance, reads in the current window). Counted traits include uniquely/non-uniquely mapped reads,
     * forward-strand aligned reads etc. 
     */
    class WindowStats implements Cloneable {
        public int uniqueReads = 0;
        public int nonUniqueReads = 0;
        public int uniqueFWDReads = 0;
        public int nonUniqueFWDReads = 0;

        /** Reset all counts to 0 */
        public void clear() {
            uniqueReads = nonUniqueReads = uniqueFWDReads = nonUniqueFWDReads = 0;
        }

        /** Examines the read and increments the counts for all the monitored traits observed in this read. */
        public void addRead(SAMRecord r) {
            if ( r.getMappingQuality() == 0 ) {
                // nonunique
                nonUniqueReads++;
                if ( ! r.getReadNegativeStrandFlag() ) nonUniqueFWDReads++;
            } else {
                // unique
                uniqueReads++;
                if ( ! r.getReadNegativeStrandFlag() ) uniqueFWDReads++;
            }
        }

        /** Examines the read and decrements the counts for all the monitored traits observed in this read. */
        public void removeRead(SAMRecord r) {
            if ( r.getMappingQuality() == 0 ) {
                // nonunique
                nonUniqueReads--;
                if ( ! r.getReadNegativeStrandFlag() ) nonUniqueFWDReads--;
            }
            else {
                // unique
                uniqueReads--;
                if ( ! r.getReadNegativeStrandFlag() ) uniqueFWDReads--;
            }

        }

        /** allocates new object, copies this object into it, and returns the copy */
        public WindowStats clone() {
            WindowStats ret = new WindowStats();
            ret.uniqueReads = this.uniqueReads;
            ret.nonUniqueReads = this.nonUniqueReads;
            ret.uniqueFWDReads = this.uniqueFWDReads;
            ret.nonUniqueFWDReads = this.nonUniqueFWDReads;
            return ret;
        }
    }
}

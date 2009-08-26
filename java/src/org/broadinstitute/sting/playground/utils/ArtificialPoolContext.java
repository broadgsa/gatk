package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.gatk.walkers.genotyper.SingleSampleGenotyper;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.GenotypeCall;

import java.io.PrintWriter;
import java.util.Set;
import java.util.List;
import java.util.LinkedList;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMFileWriter;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Aug 26, 2009
 * Time: 11:37:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class ArtificialPoolContext {
    private PrintWriter writerToAuxiliaryFile;
    private SAMFileWriter writerToSamFile;
    private SingleSampleGenotyper ssg;
    private List<Set<String>> readGroupSets;
    private long[] runningCoverage;
    private RefMetaDataTracker refTracker;
    private ReferenceContext refContext;
    private AlignmentContext aliContext;

    public ArtificialPoolContext() {
        readGroupSets = null;
        writerToAuxiliaryFile = null;
        writerToSamFile = null;
        ssg = null;
        refTracker = null;
        aliContext = null;
        refContext = null;
        runningCoverage = null;
    }

    public ArtificialPoolContext(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        refTracker = tracker;
        refContext = ref;
        aliContext = context;
        readGroupSets = null;
        writerToAuxiliaryFile = null;
        writerToSamFile=null;
        ssg = null;
        runningCoverage = null;
    }

    public ArtificialPoolContext(PrintWriter pw, SAMFileWriter sw, SingleSampleGenotyper g, List<Set<String>> rgs, long [] runcvg, RefMetaDataTracker rt, ReferenceContext rc, AlignmentContext ac) {
        writerToAuxiliaryFile = pw;
        writerToSamFile = sw;
        ssg = g;
        readGroupSets = rgs;
        runningCoverage = runcvg;
        refTracker = rt;
        refContext = rc;
        aliContext = ac;
    }

    public void setAuxWriter(PrintWriter writer) {
        writerToAuxiliaryFile = writer;
    }

    public void setSingleSampleGenotyper(SingleSampleGenotyper typer) {
        ssg = typer;
    }

    public void initializeSSG() {
        ssg.initialize();
    }

    public void setReadGroupSets(List<Set<String>> rgSets) {
        readGroupSets = rgSets;
    }

    public void setRefMetaDataTracker(RefMetaDataTracker tracker) {
        refTracker = tracker;
    }

    public void setReferenceContext(ReferenceContext ref) {
        refContext = ref;
    }

    public void setAlignmentContext(AlignmentContext context) {
        aliContext = context;
    }

    public void setRunningCoverage(long[] estimate) {
        runningCoverage = estimate;
    }

    public void setSAMFileWriter(SAMFileWriter writer) {
        writerToSamFile = writer;
    }

    public int getTotalNumberOfPeople() {
        return readGroupSets.size();
    }

    public RefMetaDataTracker getRefMetaDataTracker() {
        return refTracker;
    }

    public ReferenceContext getReferenceContext() {
        return refContext;
    }

    public AlignmentContext getAlignmentContext() {
        return aliContext;
    }

    public PrintWriter getWriterToAuxiliaryFile() {
        return writerToAuxiliaryFile;
    }

    public SingleSampleGenotyper getSingleSampleGenotyper() {
        return ssg;
    }

    public List<Set<String>> getReadGroupSets() {
        return readGroupSets;
    }

    public long[] getRunningCoverage() {
        return runningCoverage;
    }

    public SAMFileWriter getSAMFileWriter() {
        return writerToSamFile;
    }

    public List<SAMRecord> getReads() {
        List<SAMRecord> reads;
        if(aliContext == null) {
            reads=null;
        } else {
            reads = aliContext.getReads();
        }

        return reads;
    }

    public List<Integer> getOffsets() {
        List<Integer> offsets;
        if(aliContext == null) {
            offsets = null;
        } else {
            offsets = aliContext.getOffsets();
        }

        return offsets;
    }

    public List<SAMRecord> getNewReads() {
        List<SAMRecord> newReads;

        if(aliContext == null) {
            newReads = null;
        } else {
            newReads = new LinkedList<SAMRecord>();
            List<SAMRecord> allReads = aliContext.getReads();
            List<Integer> allOffsets = aliContext.getOffsets();

            for(int iter = 0; iter < allReads.size(); iter++) {
                if(allOffsets.get(iter) == 0) {
                    newReads.add(allReads.get(iter));
                }
            }
        }

        return newReads;
    }

    public Pair<List<SAMRecord>[],List<Integer>[]> splitByGroup(List<SAMRecord> unsplitReads, List<Integer> unsplitOffsets) {
        List<SAMRecord>[] readsSplitByGroup;
        List<Integer> [] offsetsSplitByGroup;

        if(unsplitReads != null && readGroupSets != null) {
            readsSplitByGroup = new ArrayList[this.getTotalNumberOfPeople()];
            if(unsplitOffsets != null) {
                offsetsSplitByGroup = new ArrayList[this.getTotalNumberOfPeople()];
            }
            else {
                offsetsSplitByGroup = null;
            }
            int listSize = unsplitReads.size();
            for(int element = 0; element < listSize; element++) {
                SAMRecord read = unsplitReads.get(element);
                for(int groupNumber = 0; groupNumber < this.getTotalNumberOfPeople(); groupNumber++) {
                    if(readGroupSets.get(groupNumber).contains((String) read.getAttribute("RG"))) {
                        readsSplitByGroup[groupNumber].add(read);
                        if(offsetsSplitByGroup != null) {
                            offsetsSplitByGroup[groupNumber].add(unsplitOffsets.get(element));
                        }
                        break;
                    }
                }
            }
        } else {
            readsSplitByGroup = null;
            offsetsSplitByGroup = null; // compiler complains without these lines
        }

        return new Pair(readsSplitByGroup,offsetsSplitByGroup);
    }

    public List<SAMRecord>[] splitReadsByGroup(List<SAMRecord> unsplitReads) {
        return (this.splitByGroup(unsplitReads,null)).first;
    }


    // Static methods follow

    public static ArtificialPoolContext mapReduceMerge(ArtificialPoolContext mapContext, ArtificialPoolContext reduceContext) {
        return new ArtificialPoolContext(reduceContext.getWriterToAuxiliaryFile(),reduceContext.getSAMFileWriter(),
                reduceContext.getSingleSampleGenotyper(), reduceContext.getReadGroupSets(), reduceContext.getRunningCoverage(),
                mapContext.getRefMetaDataTracker(),mapContext.getReferenceContext(),mapContext.getAlignmentContext());
    }

    public static Pair<List<SAMRecord>[],List<Integer>> sampleReadsAndOffsets(List<SAMRecord>[] reads, List<Integer>[] offsets, double[] propEstGlobal) {
        double[] samplingRate = calculateSamplingRateFromGlobalEstimate(propEstGlobal);
        List<SAMRecord>[] sampledReads = new ArrayList[reads.length];
        List<Integer>[] sampledOffsets;
        if(offsets != null){
            sampledOffsets = new ArrayList[offsets.length];
        } else {
            sampledOffsets = null;
        }

        for(int group = 0; group < reads.length; group++) {
            for(int readNumber = 0; readNumber < reads[group].size(); readNumber++) {
                if(Math.random() < samplingRate[group]) {
                    sampledReads[group].add(reads[group].get(readNumber));
                    if(sampledOffsets != null) {
                        sampledOffsets[group].add(offsets[group].get(readNumber));
                    }
                }
            }
        }

        return new Pair(sampledReads,sampledOffsets);
    }

    public String genotypeAndConfidenceToString(int group, String spacer) {
        GenotypeCall call = this.getGenotypeCall(group);
        return (call.getGenotypes() + spacer + call.getConfidenceScore().toString());
    }

    public GenotypeCall getGenotypeCall(int group) {
        AlignmentContext alicon = this.getAlignmentContext();
        Pair<List<SAMRecord>[],List<Integer>[]> byGroupSplitPair = this.splitByGroup(alicon.getReads(),alicon.getOffsets());
        return ssg.map(this.getRefMetaDataTracker(),this.getReferenceContext(),
                new AlignmentContext(this.getAlignmentContext().getLocation(), byGroupSplitPair.first[group],byGroupSplitPair.second[group]));
    }

    public static List<SAMRecord>[] sampleReads(List<SAMRecord>[] reads, double[] propEstGlobal) {
        return (sampleReadsAndOffsets(reads, null, propEstGlobal)).first;
    }

    public static double[] calculateSamplingRateFromGlobalEstimate(double[] ratios) {
        double min = ratios[0];
        for(double ratio : ratios) {
            if(ratio < min) {
                min = ratio;
            }
        }
        double[] samplingRate = new double[ratios.length];
        // now divide by minimum
        for(int j = 0; j < ratios.length; j++) {
            samplingRate[j] = ratios[j]/min;
        }

        return samplingRate;
    }

}

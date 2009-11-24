package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;
                                       
import java.util.List;
import java.util.Collection;
import java.io.IOException;
import java.io.File;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import net.sf.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.util.Interval;
import edu.mit.broad.picard.util.OverlapDetector;
import edu.mit.broad.picard.directed.IntervalList;

@By(DataSource.REFERENCE)
public class HybSelPerformanceWalker extends LocusWalker<Integer, HybSelPerformanceWalker.TargetInfo> {
    @Argument(fullName="min_mapq", shortName="mmq", required=false, doc="Minimum mapping quality of reads to consider")
    public Integer MIN_MAPQ = 1;

    @Argument(fullName="include_duplicates", shortName="idup", required=false, doc="consider duplicate reads")
    public boolean INCLUDE_DUPLICATE_READS = false;

    @Argument(fullName="free_standing_distance", shortName="fsd", required=false, doc="minimum distance to next interval to consider freestanding")
    public Integer FREE_STANDING_DISTANCE = 500;

    @Argument(fullName="booster", required=false, doc="interval list of booster baits")
    public File BOOSTER_FILE;

    @Argument(fullName="booster_distance", required=false, doc="distance up to which a booster can affect a target")
    public Integer BOOSTER_DISTANCE = 100; // how far away can a booster be to "hit" its target?

    @Argument(fullName="bait_quantity", required=false, doc="interval list of baits with quantity/concentration (obtainied via sequencing, Nanostring)")
    public File BAIT_QUANT_FILE;

    @Argument(fullName="refseq", shortName="refseq",
            doc="Name of RefSeq transcript annotation file. If specified, intervals will be specified with gene names", required=false)
    String REFSEQ_FILE = null;

    private SeekableRODIterator<rodRefSeq> refseqIterator=null;

    public static class TargetInfo {
        public int counts = 0;

        // did at least two reads hit this target
        public boolean hitTwice = false;

        public int positionsOver2x = 0;
        public int positionsOver10x = 0;
        public int positionsOver20x = 0;
        public int positionsOver30x = 0;
    }

//    @Argument(fullName="suppressLocusPrinting",required=false,defaultValue="false")
//    public boolean suppressPrinting;

    @Override
    public void initialize() {
        if ( REFSEQ_FILE != null ) {
            ReferenceOrderedData<rodRefSeq> refseq = new ReferenceOrderedData<rodRefSeq>("refseq",
                    new java.io.File(REFSEQ_FILE), rodRefSeq.class);

            refseqIterator = refseq.iterator();
            logger.info("Using RefSeq annotations from "+REFSEQ_FILE);
        }

        if ( refseqIterator == null ) logger.info("No annotations available");

    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<SAMRecord> reads = context.getReads();

        int depth = 0;
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag()) { continue; }
            if (read.getReadUnmappedFlag()) { continue; }
            if (!INCLUDE_DUPLICATE_READS && read.getDuplicateReadFlag()) { continue; }
            if (read.getMappingQuality() < MIN_MAPQ) { continue; }
            depth++;
        }

        return depth;
    }

    /**
     * Return true if your walker wants to reduce each interval separately.  Default is false.
     *
     * If you set this flag, several things will happen.
     *
     * The system will invoke reduceInit() once for each interval being processed, starting a fresh reduce
     * Reduce will accumulate normally at each map unit in the interval
     * However, onTraversalDone(reduce) will be called after each interval is processed.
     * The system will call onTraversalDone( GenomeLoc -> reduce ), after all reductions are done,
     *   which is overloaded here to call onTraversalDone(reduce) for each location
     */
    public boolean isReduceByInterval() {
        return true;
    }

    public TargetInfo reduceInit() { return new TargetInfo(); }

    public TargetInfo reduce(Integer depth, TargetInfo sum) {
        sum.counts += depth;
        if (depth >= 2) { sum.hitTwice = true; sum.positionsOver2x++;}
        if (depth >= 10) { sum.positionsOver10x++; }
        if (depth >= 20) { sum.positionsOver20x++; }
        if (depth >= 30) { sum.positionsOver30x++; }
        return sum;
    }

    public void onTraversalDone(TargetInfo result) {
    }

    @Override
    public void onTraversalDone(List<Pair<GenomeLoc, TargetInfo>> results) {
        out.println("location\tlength\tgc\tavg_coverage\tnormalized_coverage\thit_twice\tfreestanding\tboosted\tbases_over_2x\tbases_over_10x\tbases_over_20x\tbases_over_30x\tgene_name\tbait_quantity\tadjacent_baits");

        // first zip through and build an overlap detector of all the intervals, so later
        // we can calculate if this interval is free-standing
        OverlapDetector<Interval> od = new OverlapDetector<Interval>(-FREE_STANDING_DISTANCE,0);
        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            Interval interval = makeInterval(target);
            od.addLhs(interval, interval);
        }

        OverlapDetector<Interval> booster = new OverlapDetector<Interval>(-BOOSTER_DISTANCE,0);
        if (BOOSTER_FILE != null) {
            IntervalList il = IntervalList.fromFile(BOOSTER_FILE);
            List<Interval> l = il.getUniqueIntervals();
            booster.addAll(l, l);
        }

        // Create a interval detector that will give us the bait quantity specified in
        // an optional bait quantity file
        OverlapDetector<Interval> bait_quant = new OverlapDetector<Interval>(0,0);
        if (BAIT_QUANT_FILE != null) {
            IntervalList il = IntervalList.fromFile(BAIT_QUANT_FILE);
            List<Interval> l = il.getIntervals();
            bait_quant.addAll(l, l);
        }

        // Create a detector of adjacent baits 
        OverlapDetector<Interval> adjacent_bait_detector = new OverlapDetector<Interval>(-1,-1);
        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            Interval interval = makeInterval(target);
            adjacent_bait_detector.addLhs(interval, interval);
        }

        // now zip through and calculate the total average coverage
        long totalCoverage = 0;
        long basesConsidered = 0;
        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            TargetInfo ti = pair.getSecond();

            // as long as it was hit twice, count it
            if(ti.hitTwice) {
                long length = target.getStop() - target.getStart() + 1;
                totalCoverage += ti.counts;
                basesConsidered += length;
            }
        }
        double meanTargetCoverage = totalCoverage / basesConsidered;


        for(Pair<GenomeLoc, TargetInfo> pair : results) {
            GenomeLoc target = pair.getFirst();
            TargetInfo ti = pair.getSecond();
            long length = target.getStop() - target.getStart() + 1;

            double avgCoverage = ((double)ti.counts / (double)length);
            double normCoverage = avgCoverage / meanTargetCoverage;

            // calculate gc for the target
            double gc = calculateGC(target);

            // if there is more than one hit on the overlap detector, it's not freestanding
            Interval targetInterval = makeInterval(target);

            Collection<Interval> hits = od.getOverlaps(targetInterval);
            boolean freestanding = (hits.size() == 1);

            boolean boosted = (booster.getOverlaps(targetInterval).size() > 0);

            // look up the gene name info
            String geneName = getGeneName(target);

            Collection<Interval> bait_quant_hits = bait_quant.getOverlaps(targetInterval);
            String bait_quant_string = (bait_quant_hits.size() == 1) ? bait_quant_hits.iterator().next().getName() : "0";
            if (bait_quant_hits.size() <= 1) { out.printf("Warning: multiple bait quantity intervals detected; perhaps bait quantity interval lengths don't match primary interval list specified with -L"); }
            int adjacent_baits = adjacent_bait_detector.getOverlaps(targetInterval).size() - 1;

            out.printf("%s:%d-%d\t%d\t%6.4f\t%6.4f\t%6.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%d\n",
                       target.getContig(), target.getStart(), target.getStop(), length, gc,
                        avgCoverage, normCoverage, ((ti.hitTwice)?1:0), ((freestanding)?1:0), ((boosted)?1:0),
                        ti.positionsOver2x,
                        ti.positionsOver10x,
                        ti.positionsOver20x,
                        ti.positionsOver30x,
                        geneName,
                        bait_quant_string,
                        adjacent_baits
                    );


        }
    }

    private Interval makeInterval(GenomeLoc target) {
        return new Interval(target.getContig(), (int) target.getStart(), (int) target.getStop());
    }

    private String getGeneName(GenomeLoc target) {
        if (refseqIterator == null) { return "UNKNOWN"; }

        RODRecordList<rodRefSeq> annotationList = refseqIterator.seekForward(target);
        if (annotationList == null) { return "UNKNOWN"; }
        
        for(rodRefSeq rec : annotationList) {
            if ( rec.overlapsExonP(target) ) {
                return rec.getGeneName();
            }
        }

        return "UNKNOWN";

    }

    IndexedFastaSequenceFile seqFile = null;

    private double calculateGC(GenomeLoc target) {
        try {
            if (seqFile == null) {
                seqFile = new IndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
            }
            ReferenceSequence refSeq = seqFile.getSubsequenceAt(target.getContig(),target.getStart(), target.getStop());


            int gcCount = 0;
            for(char base : StringUtil.bytesToString(refSeq.getBases()).toCharArray()) {
                if (base == 'C' || base == 'c' || base == 'G' || base == 'g') { gcCount++; }
            }
            return ( (double) gcCount ) / ((double) refSeq.getBases().length);

        } catch (IOException ioe) {
            throw new RuntimeException(ioe);
        }

    }
}
package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.IntervalRod;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Oct 20, 2009
 * Time: 5:00:53 PM
 * To change this template use File | Settings | File Templates.
 */

@By(DataSource.REFERENCE)
public class CoverageAcrossIntervalsWalker extends LocusWalker<Pair<Integer, Integer>, CoverageAcrossIntervalsWalker.Coverage> {
    /* Accumulates coverage across hybrid selection bait intervals to assess effect of bait adjacency.
     Requires input bait intervals that have an overhang beyond the actual bait interval to capture coverage data at these points.
     Outputs R parseable file that has all data in lists and then does some basic plotting.
      */
    @Argument(fullName="include_duplicates", shortName="idup", required=false, doc="consider duplicate reads")
    public boolean INCLUDE_DUPLICATE_READS = false;

    @Argument(fullName="min_mapq", shortName="mmq", required=false, doc="Minimum mapping quality of reads to consider")
    public Integer MIN_MAPQ = 1;

    @Argument(fullName="min_len", shortName="min", required=false, doc="Minimum length of interval to consider")
    public Integer MIN_LEN = Integer.MIN_VALUE;

    @Argument(fullName="max_len", shortName="max", required=false, doc="Maximum length of interval to consider")
    public Integer MAX_LEN = 5000;

    @Argument(fullName="max_bait_count", shortName="mbc", required=false, doc="Maximum number of baits to consider")
    public Integer MAX_BAIT_COUNT = 24;

    @Argument(fullName="free_standing_distance", shortName="fsd", required=false, doc="minimum distance to next interval to consider freestanding")
    public Integer FREE_STANDING_DISTANCE = 500;

    public static class Coverage {
        // Container for carrying positive and negative strand coverage
        public ArrayList<Integer> pos = new ArrayList<Integer>();
        public ArrayList<Integer> neg = new ArrayList<Integer>();
        public void addCoveragePoints(Pair<Integer,Integer> pn) {
            pos.add(pn.first);
            neg.add(pn.second);
        }
    }

    public Pair<Integer, Integer> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<SAMRecord> reads = context.getReads();

        IntervalRod intervalROD = (IntervalRod)tracker.lookup("interval", null);
        GenomeLoc interval = intervalROD == null ? null : intervalROD.getLocation();
        if (interval == null) { throw new StingException("No intervals at locus; should not happen"); }
        int offset = (int)(context.getPosition() - interval.getStart());

        int depth[] = new int[2];
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag()) { continue; }
            if (read.getReadUnmappedFlag()) { continue; }
            if (!INCLUDE_DUPLICATE_READS && read.getDuplicateReadFlag()) { continue; }
            if (read.getMappingQuality() < MIN_MAPQ) { continue; }

            depth[read.getReadNegativeStrandFlag() ? 0 : 1]++;
        }

        return new Pair<Integer,Integer>(depth[0],depth[1]);
    }

    public Coverage reduceInit() { return new Coverage(); }
    public Coverage reduce(Pair<Integer, Integer> depths, Coverage cov) {
        cov.addCoveragePoints(depths);
        return cov;
    }

    /**
     * Return true if your walker wants to reduce each interval separately.  Default is false.
     * If you set this flag, several things will happen.
     * The system will invoke reduceInit() once for each interval being processed, starting a fresh reduce
     * Reduce will accumulate normally at each map unit in the interval
     * However, onTraversalDone(reduce) will be called after each interval is processed.
     * The system will call onTraversalDone( GenomeLoc -> reduce ), after all reductions are done,
     *   which is overloaded here to call onTraversalDone(reduce) for each location
     */
    public boolean isReduceByInterval() {
        return true;
    }

    public void onTraversalDone(List<Pair<GenomeLoc, Coverage>> results) {
        int TOTAL_BAIT_WINDOW = MAX_BAIT_COUNT*120+FREE_STANDING_DISTANCE*2;
        int depths[][][] = new int[2][MAX_BAIT_COUNT+1][TOTAL_BAIT_WINDOW]; // Indices: pos/neg strand count, bait count, position in bait
        int bait_count[] = new int[MAX_BAIT_COUNT+1];

        for (Pair<GenomeLoc, Coverage> pair : results) {
            GenomeLoc loc = pair.first;
            //if (loc.size() < MIN_LEN || loc.size() > MAX_LEN) { continue; }

            long interval_width = loc.size() - FREE_STANDING_DISTANCE * 2;
            int baits = (int)(interval_width / 120);
            //out.format("#Interval_width: %d Bait_count: %d\n", interval_width, baits);
            if (interval_width % 120 != 0) { continue; }
            if (interval_width > MAX_BAIT_COUNT*120) { continue; }

            // This target is good; count it
            bait_count[baits]++;
            Coverage pn_cov = pair.second;
            for (int pn=0; pn<2; pn++) {
                ArrayList<Integer> cov = pn == 0 ? pn_cov.neg : pn_cov.pos;
                for (int i=0; i<cov.size(); i++) {
                    depths[pn][baits][i] += cov.get(i);
                }
            }
        }

        // Write out the R script that plots all of the coverage across intervals
        // All of the data is contained in variables b? where ? corresponds
        // to 
        Boolean plot_begun = false;
        out.format("pcounts <- NULL\n");
        out.format("ncounts <- NULL\n");
        for (int baits=MAX_BAIT_COUNT; baits>0; baits--) {
            int norm_depth_width = baits*120+FREE_STANDING_DISTANCE*2;
            float norm_depths[] = new float[norm_depth_width];
            for (int pn=0; pn<2; pn++) {
                for (int i=0; i<norm_depth_width; i++) {
                    norm_depths[i] = (float)depths[pn][baits][i]/bait_count[baits];
                }
                //out.format("Baits/target: %s\n", baits);
                //out.format("Depths: %s\n", Arrays.toString(norm_depths));
                //out.format("Number of baits: %d\n", bait_count[baits]);
                char pn_char = pn == 0 ? 'n' : 'p';
                String pn_color = pn == 0 ? "black" : "red";
                if (bait_count[baits] > 0) {
                    String depth_str = Arrays.toString(norm_depths);
                    depth_str = depth_str.substring(1,depth_str.length()-1);
                    out.format("b%c%d <- c(%s)\n", pn_char, baits, depth_str);
                    if (plot_begun) {
                        out.format("points(b%c%d, type='l', col=\"%s\")\n", pn_char, baits, pn_color);
                    }else{
                        out.format("plot(b%c%d, type='l', col=\"%s\", xlab='Position from start of 1st bait', ylab='Coverage')\n", pn_char, baits, pn_color);
                        plot_begun = true;
                    }
                }
                out.format("%ccounts <- c(%ccounts,%d)\n", pn_char, pn_char, bait_count[baits]);
            }
        }
    }

}


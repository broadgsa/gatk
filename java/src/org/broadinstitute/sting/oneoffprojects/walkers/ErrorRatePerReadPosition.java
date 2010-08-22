package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;

import java.util.List;
import java.util.HashMap;
import java.util.HashSet;
import java.io.PrintStream;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

/**
 * Computes the read error rate per position in read (in the original 5'->3' orientation that the read had coming off the machine)
 */
public class ErrorRatePerReadPosition extends LocusWalker<Integer, Integer> {
    @Output PrintStream out;
    @Argument(fullName="min_base_quality_score", shortName="mbq", doc="Minimum base quality required to consider a base for calling (default: 0)", required=false) public Integer MIN_BASE_QUAL = 0;
    @Argument(fullName="min_mapping_quality_score", shortName="mmq", doc="Minimum read mapping quality required to consider a read for calling (default: 0)", required=false) public Integer MIN_MAPPING_QUAL = 0;

    private HashMap<String, int[]> mismatches;
    private HashMap<String, int[]> counts;
    private HashMap<String, int[]> quals;
    private HashMap<String, HashSet<Integer>> readLengthsPerReadGroup;
    private HashSet<Integer> readLengths;
    private int readLength = 10000;

    public void initialize() {
        mismatches = new HashMap<String, int[]>();
        counts = new HashMap<String, int[]>();
        quals = new HashMap<String, int[]>();
        readLengthsPerReadGroup = new HashMap<String, HashSet<Integer>>();
        readLengths = new HashSet<Integer>();

        for (SAMReadGroupRecord rg : this.getToolkit().getSAMFileHeader().getReadGroups()) {
            int[] mm = new int[readLength];
            int[] cm = new int[readLength];
            int[] qm = new int[readLength];

            mismatches.put(rg.getReadGroupId(), mm);
            counts.put(rg.getReadGroupId(), cm);
            quals.put(rg.getReadGroupId(), qm);
        }
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        List<Integer> offsets = context.getOffsets();
        List<SAMRecord> reads = context.getReads();

        for (int i = 0; i < offsets.size(); i++) {
            int offset = offsets.get(i);

            if (reads.get(i).getMappingQuality() >= MIN_MAPPING_QUAL &&
                reads.get(i).getBaseQualities()[offset] >= MIN_BASE_QUAL) {

                char readBase = reads.get(i).getReadString().charAt(offset);

                int refIndex = ref.getBaseIndex();
                int readIndex = BaseUtils.simpleBaseToBaseIndex(readBase);

                if (!reads.get(i).getReadNegativeStrandFlag() && (!reads.get(i).getReadPairedFlag() || reads.get(i).getFirstOfPairFlag())) {
                    String keyName = reads.get(i).getReadGroup().getReadGroupId();

                    mismatches.get(keyName)[offset] += (refIndex != readIndex) ? 1 : 0;
                    counts.get(keyName)[offset]++;
                    quals.get(keyName)[offset] += reads.get(i).getBaseQualities()[offset];

                    int readLength = reads.get(i).getReadLength();
                    if (!readLengthsPerReadGroup.containsKey(keyName)) {
                        readLengthsPerReadGroup.put(keyName, new HashSet<Integer>());
                    }

                    readLengthsPerReadGroup.get(keyName).add(readLength);
                    readLengths.add(readLength);
                }
            }
        }

        return null;
    }

    public Integer reduceInit() {
        return null;
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        String[] rgs = mismatches.keySet().toArray(new String[1]);

        out.printf("position");
        for (String rg : rgs) { out.printf("\t%s", rg); }
        //out.printf("\tmean\tsd\tmin\tmax");

        //for (int readLength : readLengths) {
        //    out.printf("\t%dreadLength", readLength);
        //}

        out.println();

        for (int i = 0; i < readLength; i++) {
            boolean print = false;
            String row = "";

            row += String.format("%d", i);

            double rsum = 0.0;
            double min = 1.00;
            double max = 0.00;
            
            for (String rg : rgs) {
                double value = ((double) mismatches.get(rg)[i])/((double) counts.get(rg)[i]);
                //double value = ((double) quals.get(rg)[i])/((double) counts.get(rg)[i]);

                if (Double.isInfinite(value) || Double.isNaN(value)) {
                    value = 0.0;
                } else {
                    print = true;
                }

                row += String.format("\t%f", value);

                rsum += value;
                if (value > max) { max = value; }
                if (value < min) { min = value; }
            }

            double mean = rsum/rgs.length;

            double squareDeviationSum = 0.0;
            for (String rg : rgs) {
                double value = ((double) mismatches.get(rg)[i])/((double) counts.get(rg)[i]);

                if (Double.isInfinite(value) || Double.isNaN(value)) {
                    value = 0.0;
                } else {
                    print = true;
                }

                squareDeviationSum += Math.pow(value - mean, 2.0);
            }

            double sd = Math.sqrt(squareDeviationSum/rgs.length);

            //row += String.format("\t%f\t%f\t%f\t%f", mean, sd, min, max);

            row += String.format("%n");

            if (print) {
                out.print(row);
            }
        }
    }
}

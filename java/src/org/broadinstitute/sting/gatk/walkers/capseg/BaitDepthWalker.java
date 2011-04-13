package org.broadinstitute.sting.gatk.walkers.capseg;

import net.sf.samtools.SAMRecord;
import org.broad.tribble.bed.FullBEDFeature;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;
import java.util.*;

/**
 * This tool captures baits depth (as reads per bait) across all baits and
 * lanes.  The output is one line per bait, with columns for each of the
 * lanes.  Reads which overlap two or more baits are added to all baits that
 * the read overlaps.
 *
 */
public class BaitDepthWalker extends ReadWalker<ReadStats, ReadStats> {

    // what sample we should assign to our unknown reads
    private static final String UNKNOWN_SAMPLE_NAME = "unknown";

    // what separator to use
    private static final String SEPERATOR = ",";

    // the current output location
    @Output(doc="File to which baits and coverage information for lanes should be written",required=true)
    PrintStream out;

    String currentBait = null;

    // our current bait name
    String baitName = null;

    // create a coverage mapping by readGroup
    HashMap<String,Integer> currentCoverage = new HashMap<String,Integer>();

    // set of all of then read groups
    Set<String> readGroups = new LinkedHashSet<String>();

    // set containing the baits and their positions, used for metrics during and after the run
    HashMap<String,GenomeLoc> baitMappings = new HashMap<String,GenomeLoc>();

    // a structure to hold reads that align to other baits (other than the one we're currently processing)
    // baits should be contiguous, so reads should not be left in this queue once we've passed the appropriate bait
    HashMap<String,ArrayList<String>> recordsForOtherBaits = new HashMap<String,ArrayList<String>>();

    // a mapping of the read group to the bam
    HashMap<String,String> readGroupToBam = new HashMap<String,String>();

    /**
     * get all of the read groups from the BAM headers, initializing our maps and arrays with this information
     */
    @Override
    public void initialize() {
        int bamNumber = 0;
        for (Set<String> set : getToolkit().getMergedReadGroupsByReaders())
            for (String smp : set) {
                readGroups.add(smp);
                readGroupToBam.put(smp,String.valueOf(bamNumber));
                bamNumber++;
            }
        out.print("#bait,location");
        for (String st : readGroups) {
            out.print(SEPERATOR);
            out.print(st);
            currentCoverage.put(st,0);
        }
        out.println();
    }

    /**
     * get the read, determine which baits it falls into, and update our matrix
     * @param ref the reference bases
     * @param read the sequencer read
     * @param metaDataTracker our metadata, used to pull in the BED file
     * @return a read stat object, which tracks how many reads fell into baits, outside of baits, etc
     */
    @Override
    public ReadStats map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker metaDataTracker) {
        // setup the return value
        ReadStats stats = new ReadStats();

        // get a list of the features
        List<GATKFeature> features = metaDataTracker.getAllCoveringRods();
        if (features.size() == 0) {
            stats.addRead(0);
        } else {
            LinkedHashSet<String> coveringBaits = new LinkedHashSet<String>();
            for (GATKFeature feature : features) {
                if (feature.getUnderlyingObject() instanceof FullBEDFeature &&
                        !((FullBEDFeature)feature.getUnderlyingObject()).getName().contains("ooster")) {
                    String name = ((FullBEDFeature)feature.getUnderlyingObject()).getName();
                    coveringBaits.add(name);
                    if (!baitMappings.containsKey(name)) {
                        baitMappings.put(name,feature.getLocation());
                    }

                }
            }
            if (coveringBaits.size() == 1)
                addToBaits(coveringBaits.iterator().next(),read.getReadGroup().getReadGroupId());
            else {
                List<String> baits = new ArrayList<String>(coveringBaits);
                for (String bait : baits)
                    if (bait.equals(currentBait))
                        addToBaits(bait,read.getReadGroup().getReadGroupId());
                    else {
                        if (!recordsForOtherBaits.containsKey(bait))
                            recordsForOtherBaits.put(bait,new ArrayList<String>());
                        recordsForOtherBaits.get(bait).add(read.getReadGroup().getReadGroupId());
                    }

            }
            stats.addRead(coveringBaits.size());
        }
        return stats;
    }

    /**
     * add a read to a particular bait name
     * @param baitName the bait name
     * @param readGroup the read group (lane) the read came from
     */
    public void addToBaits(String baitName, String readGroup) {
        // if we haven't seen the bait name yet, create a hashmap for it and dump the currentContents to disk
        if (currentBait == null || !currentBait.equals(baitName)) {
            if (currentBait != null) {
                out.print(currentBait + SEPERATOR);
                currentBait = baitName;
                out.print(baitMappings.get(baitName));
                for (String rg : readGroups) {
                    out.print(SEPERATOR + currentCoverage.get(rg));
                    currentCoverage.put(rg,0);
                }

                out.println();
                if (recordsForOtherBaits.containsKey(baitName)) {
                    ArrayList<String> records = recordsForOtherBaits.get(baitName);
                    recordsForOtherBaits.remove(baitName);
                    for (String rec : records)
                        addToBaits(baitName,rec);
                }

            }
            currentBait = baitName;

        }
        if (!readGroups.contains(readGroup))
            throw new IllegalArgumentException("Novel (unseen) sample name " + readGroup + " is not in the header of any of the BAM files");
        else
            currentCoverage.put(readGroup, currentCoverage.get(readGroup) + 1);

    }

    @Override
    public ReadStats reduceInit() {
        return new ReadStats();
    }

    @Override
    public ReadStats reduce(ReadStats value, ReadStats sum) {
        return sum.add(value);
    }

    /**
     * emit the last bait, any remaining baits, and the read stats
     * @param result the reads stats object
     */
    public void onTraversalDone(ReadStats result) {
        out.print(baitName + SEPERATOR);
        out.print(baitMappings.get(baitName));
        for (String rg : this.readGroups)
            out.print(SEPERATOR + currentCoverage.get(rg));
        out.println();
        /*for (Map.Entry<String,ArrayList<String>> entry : recordsForOtherBaits.entrySet()) {
            out.print(entry.getKey() + SEPERATOR);
            out.print(baitMappings.get(baitName));
            for (String rg : this.readGroups) {
                int index = entry.getValue().indexOf(rg);
                ///out.print(SEPERATOR + (index < 0 ? 0 : entry.getValue().));
            }
            out.println();
        }*/

        currentBait = baitName;
        currentCoverage.clear();
        out.flush();
        out.close();
        logger.info("[REDUCE RESULT] Traversal result is: " + result);
    }

}

// a small helper class for storing read statistics
class ReadStats {
    int numberInBait = 0;
    int numberOutsideOfBait = 0;
    int baitToReadCount = 0; // the number of baits we've seen, used to determine the bait/read ratio
    int readCount = 0;

    public void addRead(int baitCount) {
        baitToReadCount += baitCount;
        numberInBait += (baitCount > 0) ? 1 : 0;
        numberOutsideOfBait += (baitCount > 0) ? 0 : 1;
        readCount++;
    }

    public String toString() {
        return "Number in baits = " + numberInBait + ", number outside of baits = " + numberOutsideOfBait + ", baits per read" + (baitToReadCount/readCount);
    }

    public org.broadinstitute.sting.gatk.walkers.capseg.ReadStats add(org.broadinstitute.sting.gatk.walkers.capseg.ReadStats stat) {
        this.numberInBait += stat.numberInBait;
        this.numberOutsideOfBait += stat.numberOutsideOfBait;
        this.baitToReadCount += stat.baitToReadCount;
        this.readCount += stat.readCount;
        return this;
    }
}
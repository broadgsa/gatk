package org.broadinstitute.sting.oneoffprojects.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.ReadMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.utils.GATKFeature;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.commandline.Output;

import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.io.PrintStream;

/**
 * validate the rods for reads
 */
public class ValidateRODForReads extends ReadWalker<Integer, Integer> {
    // a mapping of the position to the count of rods
    HashMap<GenomeLoc, Integer> map = new LinkedHashMap<GenomeLoc, Integer>();

    @Output
    private PrintStream out;

    @Override
    public Integer map(ReferenceContext ref, SAMRecord read, ReadMetaDataTracker tracker) {
        if (tracker != null) {
            Map<Long, Collection<GATKFeature>> mapping = tracker.getContigOffsetMapping();
            for (Map.Entry<Long, Collection<GATKFeature>> entry : mapping.entrySet()) {
                GenomeLoc location = GenomeLocParser.createGenomeLoc(read.getReferenceIndex(),entry.getKey());
                if (!map.containsKey(location)) {
                    map.put(location,0);
                }
                map.put(location,map.get(location)+1);
            }

            return mapping.size();
        }
        return 0;
    }

    @Override
    public Integer reduceInit() {
        return 0;
    }

    @Override
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer result) {
        out.println("[REDUCE RESULT] Traversal result is: " + result + " ROD entries seen");
        for (GenomeLoc location : map.keySet()) {
            out.println(location + " -> " + map.get(location));
        }
    }
}

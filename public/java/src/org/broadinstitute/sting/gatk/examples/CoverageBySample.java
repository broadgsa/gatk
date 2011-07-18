
package org.broadinstitute.sting.gatk.examples;

import net.sf.samtools.SAMReadGroupRecord;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

/**
 * Computes the coverage per sample.
 */
public class CoverageBySample extends LocusWalker<Integer, Integer> {
    @Output
    protected PrintStream out;    

    private HashSet<String> sampleNames = new HashSet<String>();

    public boolean requiresReads() { return true; }

    public void initialize() {

        List<SAMReadGroupRecord> read_groups = this.getToolkit().getSAMFileHeader().getReadGroups();

        for ( SAMReadGroupRecord record : read_groups ) {
            String sample = record.getSample();
            if ( sample != null )
                sampleNames.add(sample);
        } 
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        HashMap<String, Integer> depthBySample = new HashMap<String, Integer>();
        for ( String sample : sampleNames )
            depthBySample.put(sample, 0);

        ReadBackedPileup pileup = context.getPileup();
        for ( PileupElement p : pileup ) {

            SAMReadGroupRecord readGroup = p.getRead().getReadGroup();
            if ( readGroup == null )
                continue;

            String sample = readGroup.getSample();
            if ( sample != null ) {
                int oldDepth = depthBySample.get(sample);
                depthBySample.put(sample, oldDepth + 1);
            }
        }

        for ( Map.Entry<String, Integer> sample : depthBySample.entrySet() ) {
            out.printf("  %s %8d%n", sample.getKey(), sample.getValue());
        }

        return 1;
    }


    public void onTraversalDone(Integer result) {
        out.println("Processed " + result + " loci.");
    }

    public Integer reduceInit() {
		return 0;
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }
}

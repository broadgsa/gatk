package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariationRod;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.SingleSampleGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.SSGenotypeCall;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import org.broadinstitute.sting.utils.genotype.Variation;

import java.util.Set;
import java.util.List;
import java.util.ArrayList;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Sep 1, 2009
 * Time: 11:04:55 AM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE, DataSource.REFERENCE_BASES, DataSource.READS},referenceMetaData={@RMD(name="child",type= VariationRod.class)})
@Allows({DataSource.READS, DataSource.REFERENCE})
@ReadFilters(ZeroMappingQualityReadFilter.class)
//, @RMD(name="parent1",type= VariationRod.class), @RMD(name="parent2",type= VariationRod.class)})

public class DeNovoSNPWalker extends RefWalker<String, Integer>{

    SingleSampleGenotyper SSG;
    private List<Set<String>> readGroupSets;

    public void initialize() {
         SSG = new SingleSampleGenotyper();
         SSG.initialize();

         readGroupSets = getToolkit().getMergedReadGroupsByReaders();
     }

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        Variation child = (Variation)tracker.lookup("child", null);
        Variation dbsnp = (Variation)tracker.lookup("dbSNP", null);
        if (child != null) {
            if (child.isSNP() && child.getNegLog10PError() > 5) { // BTR > 5

                List<SAMRecord> reads = context.getReads();
                List<Integer> offsets = context.getOffsets();

                List<SAMRecord> parent1_reads = new ArrayList<SAMRecord>();
                List<SAMRecord> parent2_reads = new ArrayList<SAMRecord>();
                List<Integer> parent1_offsets = new ArrayList<Integer>();
                List<Integer> parent2_offsets = new ArrayList<Integer>();

                assert( reads.size() == offsets.size() ); // should be same number or we're in trouble
                int num_reads = reads.size();
                for (int i=0; i<num_reads; i++) {
                    SAMRecord read = reads.get(i);
                    Integer offset = offsets.get(i);
                    String readGroup = (String)read.getAttribute("RG");
                    if (readGroupSets.get(0).contains(readGroup)) {
                        parent1_reads.add(read);
                        parent1_offsets.add(offset);
                    } else if (readGroupSets.get(1).contains(readGroup)) {
                        parent2_reads.add(read);
                        parent2_offsets.add(offset);
                    }
                }

                AlignmentContext parent1_subContext = new AlignmentContext(context.getLocation(), parent1_reads, parent1_offsets);
                SSGenotypeCall parent1 = SSG.map(tracker, ref, parent1_subContext);

                AlignmentContext parent2_subContext = new AlignmentContext(context.getLocation(), parent2_reads, parent2_offsets);
                SSGenotypeCall parent2 = SSG.map(tracker, ref, parent2_subContext);

                if (!parent1.isVariant(parent1.getReference()) &&
                    parent1.getNegLog10PError() > 5 &&
                    !parent2.isVariant(parent2.getReference()) &&
                    parent2.getNegLog10PError() > 5
                ) {

                    double sumConfidences = 0.5 * (0.5 * child.getNegLog10PError() +
                            Math.min(parent1.getNegLog10PError(), parent2.getNegLog10PError()));

                    out.format("%s\t", child.getLocation().getContig());
                    out.format("%s\t", child.getLocation().getStart());
                    out.format("%.4f\t", sumConfidences);
                    out.format("%.4f\t", child.getNegLog10PError());
                    out.format("%.4f\t", parent1.getNegLog10PError());
                    out.format("%.4f\t", parent2.getNegLog10PError());
                    out.format("%s\t", dbsnp != null);

                    out.format ("%s\t", child.toString());
                    out.format ("%s\t", parent1.toString());
                    out.format ("%s", parent2.toString());
                    if (dbsnp != null)
                        out.format ("\tDBSNP\t:%s", dbsnp.toString());
                    out.println();
                }
            }
        }

        return "";
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce(String line, Integer a) {
        return 1;
    }

    public void onTraversalDone(Integer result) {} // Don't print the reduce result
}

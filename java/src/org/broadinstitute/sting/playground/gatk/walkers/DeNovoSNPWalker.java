package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariationRod;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.VariantCallContext;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Implements a rudimentary algorithm for calling SNPs that are de novo in that they appear in an individual but not in
 * its parents.
 */

@By(DataSource.REFERENCE)
@Requires(value={DataSource.REFERENCE, DataSource.REFERENCE_BASES, DataSource.READS},referenceMetaData={@RMD(name="child",type= VariationRod.class)})
@Allows({DataSource.READS, DataSource.REFERENCE})
@ReadFilters(ZeroMappingQualityReadFilter.class)
//, @RMD(name="parent1",type= VariationRod.class), @RMD(name="parent2",type= VariationRod.class)})

public class DeNovoSNPWalker extends RefWalker<String, Integer>{
/**
 * Implements a rudimentary algorithm for calling SNPs that are de novo in that they appear in an individual but not in
 * its parents. Using BAM files corresponding to parents and child, it calls UnifiedGenotyper directly and outputs a
 * confidence for positions being de novo SNPs.
 * */

    UnifiedGenotyper UG;
    private List<Set<String>> readGroupSets;

    public void initialize() {
         UG = new UnifiedGenotyper();
         UG.initialize();

         readGroupSets = getToolkit().getMergedReadGroupsByReaders();
     }

    public String map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        VariationRod child = tracker.lookup("child",VariationRod.class);
        VariationRod dbsnp = tracker.lookup("dbSNP",VariationRod.class);

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
                VariantCallContext parent1 = UG.map(tracker, ref, parent1_subContext);

                AlignmentContext parent2_subContext = new AlignmentContext(context.getLocation(), parent2_reads, parent2_offsets);
                VariantCallContext parent2 = UG.map(tracker, ref, parent2_subContext);

                if ( parent1 != null && parent1.vc != null && parent2 != null && parent2.vc != null ) {
                    Genotype parent1call = parent1.vc.getGenotype(0);
                    Genotype parent2call = parent2.vc.getGenotype(0);

                    if (parent1call.isHomRef() &&
                        parent1call.getNegLog10PError() > 5 &&
                        !parent2call.isHomRef() &&
                        parent2call.getNegLog10PError() > 5) {

                        double sumConfidences = 0.5 * (0.5 * child.getNegLog10PError() +
                                Math.min(parent1call.getNegLog10PError(), parent2call.getNegLog10PError()));

                        out.format("%s\t", child.getLocation().getContig());
                        out.format("%s\t", child.getLocation().getStart());
                        out.format("%.4f\t", sumConfidences);
                        out.format("%.4f\t", child.getNegLog10PError());
                        out.format("%.4f\t", parent1call.getNegLog10PError());
                        out.format("%.4f\t", parent2call.getNegLog10PError());
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
        }

        return "";
    }

    public Integer reduceInit() { return 0; }
    public Integer reduce(String line, Integer a) {
        return 1;
    }

    public void onTraversalDone(Integer result) {} // Don't print the reduce result
}

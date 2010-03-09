package org.broadinstitute.sting.playground.gatk.walkers.contamination;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.utils.NamedTable;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

import cern.jet.stat.Probability;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * FindContaminatingReadGroupsWalker lists read groups in a single-sample BAM file that appear
 * to be contaminants (meaning a read group that's not actually associated with the sample) by searching
 * for evidence of systematic underperformance at likely homozygous-variant sites.
 *
 * @author Kiran Garimella
 */
public class FindContaminatingReadGroupsWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="balance", shortName="bal", doc="The expected alternate allele balance for homozygous-variant sites", required=false)
    private Double BALANCE = 0.95;

    @Argument(fullName="limit", shortName="lim", doc="The pValue limit for which a read group will be deemed to be a contaminant", required=false)
    private Double LIMIT = 1e-9;

    @Argument(fullName="scaleForSample", shortName="scale", doc="the scale by which the pvalue limit should reduce for testing samples directly. "+
              "E.g. if a sample has three 1e-3 read groups, pvalue is 1e-9 -- significant; so the scale should reduce by some multiplicative factor"+
              "For each read group associated with the sample. Defaults to 1e-4 [1e-9 for 1 RG, 1e-13 for 2 RG, 1e-17 for 3, etc]", required=false)
    private Double SCALE = 1e-4;

    private UnifiedGenotyperEngine ug;
    private NamedTable altTable;
    private final double EPSILON = 1e-20;

    public void initialize() {
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        uac.CONFIDENCE_THRESHOLD = 50;
        ug = new UnifiedGenotyperEngine(getToolkit(), uac);

        altTable = new NamedTable();
    }

    /**
     * Identify likely homozygous-variant sites that are called as
     * heterozygous, so that we can isolate our inspection to these sites.
     *
     * @param tracker  the meta-data tracker
     * @param ref      information regarding the reference
     * @param context  information regarding the reads
     * @return true if this site is a suspicious het, false if otherwise
     */
    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        int altCount = 0;
        int totalCount = 0;

        ReadBackedPileup pileup = context.getBasePileup();
        int refIndex = BaseUtils.simpleBaseToBaseIndex(ref.getBase());

        for (byte base : pileup.getBases() ) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex((char) base);

            if (baseIndex != refIndex) {
                altCount++;
            }
            totalCount++;
        }

        double altBalance = ((double) altCount)/((double) totalCount);

        if (altBalance > 0.70) {
            VariantCallContext ugResult = ug.runGenotyper(tracker, ref, context);

            if (ugResult != null && ugResult.vc != null && ugResult.vc.getNSamples() > 0) {
                return ugResult.vc.getGenotype(0).isHet();
            }
        }

        return false;
    }

    /**
     * For each read group represented in the pileup, determine the fraction of bases supporting the alternate allele
     *
     * @param tracker  the meta-data tracker
     * @param ref      information regarding the reference
     * @param context  information regarding the reads
     * @return 1
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        NamedTable alleleCounts = new NamedTable();

        int refIndex = BaseUtils.simpleBaseToBaseIndex(ref.getBase());
        String colName = String.format("%s.%d", context.getContig(), context.getPosition());

        for (int i = 0; i < context.size(); i++) {
            SAMRecord read = context.getReads().get(i);
            int offset = context.getOffsets().get(i);

            SAMReadGroupRecord rg = read.getReadGroup();
            int alleleIndex = BaseUtils.simpleBaseToBaseIndex((char) read.getReadBases()[offset]);

            alleleCounts.increment(rg.getReadGroupId(), (alleleIndex == refIndex) ? "ref" : "alt");
        }

        for (String rg : alleleCounts.getRowNames()) {
            double altCount = alleleCounts.get(rg, "alt");
            double refCount = alleleCounts.get(rg, "ref");

            altTable.set(rg, colName, altCount / (altCount + refCount));
        }

        return 1;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    public Integer reduceInit() {
        return null;
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return accumulator with result of the map taken into account.
     */
    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    /**
     * Perform the t-test and list the read groups that are significant underperformers.
     *
     * @param result  the number of suspicious sites we're inspecting (this argument is ignored)
     */
    public void onTraversalDone(Integer result) {
        //out.println("readgroup\tpvalue\tstatus\tbalances");
        out.printf("%-10s\t%-13s\t%-10s\t%-10s%n", "readgroup", "pvalue", "status", "balances");

        HashMap<String,Double> pvalByReadGroup = new HashMap<String,Double>();
        for (String rg : altTable.getRowNames()) {
            String balances = "";

            // Compute mean
            double sum = 0.0, total = 0.0;

            for (String locus : altTable.getColumnNames()) {
                double value = altTable.get(rg, locus);

                sum += value;
                total += 1.0;

                balances += String.format("%2.2f,", value);
            }

            double mean = sum/total;

            // Compute stdev
            double squareSumOfMeanDifferences = 0.0;

            for (String locus : altTable.getColumnNames()) {
                double value = altTable.get(rg, locus);

                squareSumOfMeanDifferences += Math.pow(value - mean, 2.0);
            }

            double stdev = Math.sqrt(squareSumOfMeanDifferences/total);

            // Compute standard error of the mean (SEM)
            double sem = stdev/Math.sqrt(total);

            // Compute test statistic t
            double t = (mean - BALANCE) / sem;

            // Degrees of freedom
            double dof = total - 1.0;

            // Compute pValue
            double pValue = Probability.studentT(dof, t);
            pValue = pValue < EPSILON ? EPSILON : pValue;
            pvalByReadGroup.put(rg,pValue);

            //out.printf("%s\t%e\t%s\t[%s]\n", rg, pValue, (pValue < LIMIT ? "aberrant" : "nominal"), balances);
            out.printf("%-10s\t%-13s\t%-10s\t[%-10s]\n",
                       rg,
                       String.format("%e", pValue),
                       (pValue < LIMIT ? "aberrant" : "nominal"),
                       balances);

            logger.debug(rg);
        }

        out.printf("%n%n%s%n","SECTION ON BADLY CONTAMINATED SAMPLES");
        out.printf("%s\t%s\t%s\t%s%n","sample","p-value","status","info");

        HashMap<String,List<String>> samplesToReadGroups = new HashMap<String,List<String>>();
        for ( SAMReadGroupRecord rec : getToolkit().getSAMFileHeader().getReadGroups() ) {
            if ( samplesToReadGroups.containsKey(rec.getSample()) ) {
                samplesToReadGroups.get(rec.getSample()).add(rec.getReadGroupId());
            } else {
                ArrayList<String> newList = new ArrayList<String>();
                newList.add(rec.getReadGroupId());
                samplesToReadGroups.put(rec.getSample(),newList);
            }
        }

        for ( String sample : samplesToReadGroups.keySet() ) {
            double p_value = 1;
            double limit = LIMIT;
            boolean containsAberrantReads = false;
            for ( String rg : samplesToReadGroups.get(sample) ) {
                double rg_pval = ( pvalByReadGroup.get(rg) == null ? 1 : pvalByReadGroup.get(rg) );
                p_value = p_value*rg_pval;
                containsAberrantReads = containsAberrantReads || rg_pval < LIMIT;
                limit = limit*SCALE;
                logger.debug(rg);
            }

            out.printf("%s\t%-13s\t%s\t%s%n", sample, String.format("%e",p_value), ( p_value < limit ? "aberrant" : "nominal"), ( containsAberrantReads ? "contains_aberrant_RG" : "no_aberrant_RG"));
        }
    }
}

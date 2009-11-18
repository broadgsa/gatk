package org.broadinstitute.sting.playground.gatk.walkers.contamination;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedArgumentCollection;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCalculationModel;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeLocusData;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.utils.NamedTable;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMReadGroupRecord;

import java.util.List;

import cern.jet.stat.Probability;

/**
 * FindContaminatingReadGroupsWalker lists read groups in a single-sample BAM file that appear
 * to be contaminants by searching for evidence of systematic underperformance at likely
 * homozygous-variant sites.  First, sites that are likely homozygous-variant but are called
 * as heterozygous are identified.  Next, per each site and read group, we compute the proportion
 * of bases in the pileup supporting an alternate allele.  Finally, a one-sample, left-tailed
 * t-test is performed with the null hypothesis being that the alternate allele distribution has
 * a mean of 0.95 and the alternate hypothesis being that the true mean is statistically
 * significantly less than expected.
 *
 * @author Kiran Garimella
 */
public class FindContaminatingReadGroupsWalker extends LocusWalker<Integer, Integer> {
    @Argument(fullName="balance", shortName="bal", doc="The expected alternate allele balance for homozygous-variant sites", required=false)
    private Double BALANCE = 0.95;

    @Argument(fullName="limit", shortName="lim", doc="The pValue limit for which a read group will be deemed to be a contaminant", required=false)
    private Double LIMIT = 1e-9;

    private UnifiedArgumentCollection uac;
    private UnifiedGenotyper ug;
    private NamedTable altTable;

    public void initialize() {
        uac = new UnifiedArgumentCollection();
        uac.genotypeModel = GenotypeCalculationModel.Model.EM_POINT_ESTIMATE;
        uac.CONFIDENCE_THRESHOLD = 50;

        ug = new UnifiedGenotyper();
        ug.initialize();
        ug.setUnifiedArgumentCollection(uac);

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

        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(), context);
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
            Pair<List<Genotype>, GenotypeLocusData> ugResult = ug.map(tracker, ref, context);

            if (ugResult != null && ugResult.first != null) {
                return ugResult.first.get(0).isHet();
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

        for (int i = 0; i < context.numReads(); i++) {
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

            //out.printf("%s\t%e\t%s\t[%s]\n", rg, pValue, (pValue < LIMIT ? "aberrant" : "nominal"), balances);
            out.printf("%-10s\t%-13s\t%-10s\t[%-10s]\n",
                       rg,
                       String.format("%e", pValue),
                       (pValue < LIMIT ? "aberrant" : "nominal"),
                       balances);
        }
    }
}

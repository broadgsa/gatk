package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.VariantContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.utils.BaseUtils;
import net.sf.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.List;


public class VECFisherStrand implements VariantExclusionCriterion {
    private double pvalueLimit = 0.00001;
    private double pValue;
    private boolean exclude;
    private ArrayList<Double> factorials = new ArrayList<Double>();

    public void initialize(String arguments) {
        if (arguments != null && !arguments.isEmpty()) {
            pvalueLimit = Double.valueOf(arguments);
        }
        factorials.add(0.0);  // add fact(0)
    }

    public void compute(VariantContextWindow contextWindow) {
        VariantContext context = contextWindow.getContext();
        rodVariants variant = context.getVariant();
        int allele1 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(0));
        int allele2 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(1));

        if (allele1 != allele2) {
            exclude = strandTest(context.getReferenceContext().getBase(), context.getAlignmentContext(useZeroQualityReads()), allele1, allele2, pvalueLimit, null);
        } else {
            exclude = false;
        }
    }

    public double inclusionProbability() {
        // A hack for now until this filter is actually converted to an empirical filter
        return exclude ? 0.0 : 1.0;
    }

    public String getStudyHeader() {
        return "FisherStrand("+pvalueLimit+")\tpvalue";
    }

    public String getStudyInfo() {
        return (exclude ? "fail" : "pass") + "\t" + pValue;
    }

    public boolean useZeroQualityReads() { return false; }

    public boolean strandTest(char ref, AlignmentContext context, int allele1, int allele2, double threshold, StringBuffer out) {
        int[][] table = getContingencyTable(context, allele1, allele2);
        if ( !variantIsHet(table) )
            return false;

        updateFactorials(table);

        double pCutoff = computePValue(table);
        //printTable(table, pCutoff);

        pValue = pCutoff;
        while (rotateTable(table)) {
            double pValuePiece = computePValue(table);

            //printTable(table, pValuePiece);

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        table = getContingencyTable(context, allele1, allele2);

        while (unrotateTable(table)) {
            double pValuePiece = computePValue(table);

            //printTable(table, pValuePiece);

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        //System.out.printf("P-cutoff: %f\n", pCutoff);
        //System.out.printf("P-value: %f\n\n", pValue);

        // optionally print out the pvalue and the alternate allele counts
        if ( out != null ) {
            int refBase = BaseUtils.simpleBaseToBaseIndex(ref);
            table = getContingencyTable(context, allele1, allele2);
            if ( allele1 != refBase )
                out.append(pValue + "\t" + table[0][0] + "\t" + table[0][1] + "\t");
            else
                out.append(pValue + "\t" + table[1][0] + "\t" + table[1][1] + "\t");
        }

        return pValue < threshold;
    }

    private static boolean variantIsHet(int[][] table) {
        return ((table[0][1] != 0 || table[0][1] != 0) && (table[1][0] != 0 || table[1][1] != 0));
    }

    private static void printTable(int[][] table, double pValue) {
        System.out.printf("%d %d; %d %d : %f\n", table[0][0], table[0][1], table[1][0], table[1][1], pValue);
    }

    private static boolean rotateTable(int[][] table) {
        table[0][0] -= 1;
        table[1][0] += 1;

        table[0][1] += 1;
        table[1][1] -= 1;

        return (table[0][0] >= 0 && table[1][1] >= 0) ? true : false;
    }

    private static boolean unrotateTable(int[][] table) {
        table[0][0] += 1;
        table[1][0] -= 1;

        table[0][1] -= 1;
        table[1][1] += 1;

        return (table[0][1] >= 0 && table[1][0] >= 0) ? true : false;
    }

    private double computePValue(int[][] table) {

        int[] rowSums = { sumRow(table, 0), sumRow(table, 1) };
        int[] colSums = { sumColumn(table, 0), sumColumn(table, 1) };
        int N = rowSums[0] + rowSums[1];

        // calculate in log space so we don't die with high numbers
        double pCutoff = factorials.get(rowSums[0])
                         + factorials.get(rowSums[1])
                         + factorials.get(colSums[0])
                         + factorials.get(colSums[1])
                         - factorials.get(table[0][0])
                         - factorials.get(table[0][1])
                         - factorials.get(table[1][0])
                         - factorials.get(table[1][1])
                         - factorials.get(N);
        return Math.exp(pCutoff);

        //pCutoff *= (double) Arithmetic.factorial(rowSums[0]);
        //pCutoff *= (double) Arithmetic.factorial(rowSums[1]);
        //pCutoff *= (double) Arithmetic.factorial(colSums[0]);
        //pCutoff *= (double) Arithmetic.factorial(colSums[1]);
        //pCutoff /= (double) Arithmetic.factorial(N);
        //
        //for (int i = 0; i < 2; i++) {
        //    for (int j = 0; j < 2; j++) {
        //        pCutoff /= (double) Arithmetic.factorial(table[i][j]);
        //    }
        //}
        //
        //return pCutoff;
    }

    private static int sumRow(int[][] table, int column) {
        int sum = 0;
        for (int r = 0; r < table.length; r++) {
            sum += table[r][column];
        }

        return sum;
    }

    private static int sumColumn(int[][] table, int row) {
        int sum = 0;
        for (int c = 0; c < table[row].length; c++) {
            sum += table[row][c];
        }

        return sum;
    }

    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @param context  the context for the locus
     * @param allele1  information for the called variant
     * @param allele2  information for the called variant
     * @return a 2x2 contingency table
     */
    private static int[][] getContingencyTable(AlignmentContext context, int allele1, int allele2) {

        int[][] table = new int[2][2];

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        for (int readIndex = 0; readIndex < reads.size(); readIndex++) {
            SAMRecord read = reads.get(readIndex);
            int offset = offsets.get(readIndex);

            int readAllele = BaseUtils.simpleBaseToBaseIndex(read.getReadString().charAt(offset));
            boolean isFW = !read.getReadNegativeStrandFlag();

            if (readAllele == allele1 || readAllele == allele2) {
                int row = (readAllele == allele1) ? 0 : 1;
                int column = isFW ? 0 : 1;

                table[row][column]++;
            }
        }

        return table;
    }

    private void updateFactorials(int[][] table) {
        // calculate in log space so we don't die with high numbers
        int max = table[0][0] + table[0][1] + table[1][0] + table[1][1];
        for (int i = factorials.size(); i <= max; i++)
            factorials.add(factorials.get(i - 1) + Math.log(i));
    }
}

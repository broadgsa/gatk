package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.rodVariants;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import net.sf.samtools.SAMRecord;

import java.util.List;

import cern.jet.math.Arithmetic;

public class VECFisherStrand implements VariantExclusionCriterion {
    public void initialize(String arguments) {
        //if (arguments != null && !arguments.isEmpty()) {}
    }

    public boolean exclude(char ref, LocusContext context, rodVariants variant) {
        int allele1 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(0));
        int allele2 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(1));

        if (allele1 != allele2) {
            int[][] table = getContingencyTable(context, variant);

            double pCutoff = computePValue(table);
            //printTable(table, pCutoff);

            double pValue = 0.0;
            while (rotateTable(table)) {
                double pValuePiece = computePValue(table);

                //printTable(table, pValuePiece);

                if (pValuePiece <= pCutoff) {
                    pValue += pValuePiece;
                }
            }

            table = getContingencyTable(context, variant);

            while (unrotateTable(table)) {
                double pValuePiece = computePValue(table);

                //printTable(table, pValuePiece);

                if (pValuePiece <= pCutoff) {
                    pValue += pValuePiece;
                }
            }

            //System.out.printf("P-cutoff: %f\n", pCutoff);
            //System.out.printf("P-value: %f\n\n", pValue);

            return pValue < 0.05;
        }

        return false;
    }

    private void printTable(int[][] table, double pValue) {
        System.out.printf("%d %d; %d %d : %f\n", table[0][0], table[0][1], table[1][0], table[1][1], pValue);
    }

    private boolean rotateTable(int[][] table) {
        table[0][0] -= 1;
        table[1][0] += 1;

        table[0][1] += 1;
        table[1][1] -= 1;

        return (table[0][0] >= 0 && table[1][1] >= 0) ? true : false;
    }

    private boolean unrotateTable(int[][] table) {
        table[0][0] += 1;
        table[1][0] -= 1;

        table[0][1] -= 1;
        table[1][1] += 1;

        return (table[0][1] >= 0 && table[1][0] >= 0) ? true : false;
    }

    private double computePValue(int[][] table) {
        double pCutoff = 1.0;

        int[] rowSums = { sumRow(table, 0), sumRow(table, 1) };
        int[] colSums = { sumColumn(table, 0), sumColumn(table, 1) };
        int N = rowSums[0] + rowSums[1];

        pCutoff *= (double) Arithmetic.factorial(rowSums[0]);
        pCutoff *= (double) Arithmetic.factorial(rowSums[1]);
        pCutoff *= (double) Arithmetic.factorial(colSums[0]);
        pCutoff *= (double) Arithmetic.factorial(colSums[1]);
        pCutoff /= (double) Arithmetic.factorial(N);

        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                pCutoff /= (double) Arithmetic.factorial(table[i][j]);
            }
        }
        return pCutoff;
    }

    private int sumRow(int[][] table, int column) {
        int sum = 0;
        for (int r = 0; r < table.length; r++) {
            sum += table[r][column];
        }

        return sum;
    }

    private int sumColumn(int[][] table, int row) {
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
     * @param variant  information for the called variant
     * @return a 2x2 contingency table
     */
    private int[][] getContingencyTable(LocusContext context, rodVariants variant) {

        int[][] table = new int[2][2];

        int allele1 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(0));
        int allele2 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(1));

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
}

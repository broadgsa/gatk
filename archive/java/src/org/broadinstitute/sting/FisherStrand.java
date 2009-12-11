package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;
import net.sf.samtools.SAMRecord;
import cern.jet.math.Arithmetic;

import java.util.List;


public class FisherStrand implements VariantAnnotation {

    public String annotate(ReferenceContext ref, ReadBackedPileup pileup, Variation variation) {

        if ( !(variation instanceof VariantBackedByGenotype) )
            return null;
        final List<Genotype> genotypes = ((VariantBackedByGenotype)variation).getGenotypes();
        if ( genotypes == null || genotypes.size() == 0 )
            return null;

        // this test doesn't make sense for homs
        Genotype genotype = VariantAnnotator.getFirstHetVariant(genotypes);
        if ( genotype == null )
            return null;

        final String genotypeStr = genotype.getBases().toUpperCase();
        if ( genotypeStr.length() != 2 )
            return null;

        int allele1 = BaseUtils.simpleBaseToBaseIndex(genotypeStr.charAt(0));
        int allele2 = BaseUtils.simpleBaseToBaseIndex(genotypeStr.charAt(1));

        Double pvalue = strandTest(pileup, allele1, allele2);
        if ( pvalue == null )
            return null;

        // use Math.abs to prevent -0's
        return String.format("%.1f", Math.abs(QualityUtils.phredScaleErrorRate(pvalue)));
    }

    public String getKeyName() { return "FisherStrand"; }

    public String getDescription() { return "FisherStrand,1,Float,\"Phred-scaled p-value Using Fisher's Exact Test to Detect Strand Bias\""; }

    public boolean useZeroQualityReads() { return false; }

    private Double strandTest(ReadBackedPileup pileup, int allele1, int allele2) {
        int[][] table = getContingencyTable(pileup, allele1, allele2);
        if ( !variantIsHet(table) )
            return null;

        double pCutoff = computePValue(table);
        //printTable(table, pCutoff);

        double pValue = pCutoff;
        while (rotateTable(table)) {
            double pValuePiece = computePValue(table);

            //printTable(table, pValuePiece);

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        table = getContingencyTable(pileup, allele1, allele2);

        while (unrotateTable(table)) {
            double pValuePiece = computePValue(table);

            //printTable(table, pValuePiece);

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        //System.out.printf("P-cutoff: %f\n", pCutoff);
        //System.out.printf("P-value: %f\n\n", pValue);

       return pValue;
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
        double pCutoff = Arithmetic.logFactorial(rowSums[0])
                         + Arithmetic.logFactorial(rowSums[1])
                         + Arithmetic.logFactorial(colSums[0])
                         + Arithmetic.logFactorial(colSums[1])
                         - Arithmetic.logFactorial(table[0][0])
                         - Arithmetic.logFactorial(table[0][1])
                         - Arithmetic.logFactorial(table[1][0])
                         - Arithmetic.logFactorial(table[1][1])
                         - Arithmetic.logFactorial(N);
        return Math.exp(pCutoff);
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
     * @param pileup   the pileup for the locus
     * @param allele1  information for the called variant
     * @param allele2  information for the called variant
     * @return a 2x2 contingency table
     */
    private static int[][] getContingencyTable(ReadBackedPileup pileup, int allele1, int allele2) {

        int[][] table = new int[2][2];

        List<SAMRecord> reads = pileup.getReads();
        List<Integer> offsets = pileup.getOffsets();

        for (int readIndex = 0; readIndex < reads.size(); readIndex++) {
            SAMRecord read = reads.get(readIndex);
            int offset = offsets.get(readIndex);

	    // skip over deletion sites
	    if ( offset == -1 )
		continue;

            int readAllele = BaseUtils.simpleBaseToBaseIndex((char)read.getReadBases()[offset]);
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

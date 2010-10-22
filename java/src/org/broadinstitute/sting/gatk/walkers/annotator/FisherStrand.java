/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.annotator;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.WorkInProgressAnnotation;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import cern.jet.math.Arithmetic;

import java.util.List;
import java.util.Map;
import java.util.Arrays;
import java.util.HashMap;


public class FisherStrand implements InfoFieldAnnotation, WorkInProgressAnnotation {
    private static final String REFFWD = "REFFWD";
    private static final String REFREV = "REFREV";
    private static final String ALTFWD = "ALTFWD";
    private static final String ALTREV = "ALTREV";
    private static final String FS = "FS";

    public Map<String, Object> annotate(RefMetaDataTracker tracker, ReferenceContext ref, Map<String, StratifiedAlignmentContext> stratifiedContexts, VariantContext vc) {
        if ( ! vc.isVariant() || vc.isFiltered() || ! vc.isBiallelic() || ! vc.isSNP() )
            return null;

        int[][] table = getContingencyTable(stratifiedContexts, vc.getReference(), vc.getAlternateAllele(0));
        Double pvalue = pValueForContingencyTable(table);
        if ( pvalue == null )
            return null;

        // use Math.abs to prevent -0's
        Map<String, Object> map = new HashMap<String, Object>();
        int phredPValue = (int)Math.round(QualityUtils.phredScaleErrorRate(pvalue));
        addInt(map, REFFWD, table[0][0]);
        addInt(map, REFREV, table[0][1]);
        addInt(map, ALTFWD, table[1][0]);
        addInt(map, ALTREV, table[1][1]);
        addInt(map, FS, phredPValue);
        return map;
    }

    private static void addInt(Map<String, Object> map, String key, int value)  {
        map.put(key, String.format("%d", value));
    }

    public List<String> getKeyNames() {
        return Arrays.asList(REFFWD,REFREV,ALTFWD,ALTREV,FS);
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(
                new VCFInfoHeaderLine(REFFWD, 1, VCFHeaderLineType.Integer, "Count of bases with REF allele, forward strand"),
                new VCFInfoHeaderLine(REFREV, 1, VCFHeaderLineType.Integer, "Count of bases with REF allele, reverse strand"),
                new VCFInfoHeaderLine(ALTFWD, 1, VCFHeaderLineType.Integer, "Count of bases with ALT allele, forward strand"),
                new VCFInfoHeaderLine(ALTREV, 1, VCFHeaderLineType.Integer, "Count of bases with ALT allele, reverse strand"),
                new VCFInfoHeaderLine(FS, 1, VCFHeaderLineType.Integer, "Integer-rounded Phred-scaled p-value using Fisher's exact test to detect strand bias"));
    }

    private Double pValueForContingencyTable(int[][] originalTable) {
        int [][] table = copyContingencyTable(originalTable);

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

        table = copyContingencyTable(originalTable);
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

    private static int [][] copyContingencyTable(int [][] t) {
        int[][] c = new int[2][2];

        for ( int i = 0; i < 2; i++ )
            for ( int j = 0; j < 2; j++ )
                c[i][j] = t[i][j];

        return c;
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

    private static double computePValue(int[][] table) {

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
     * @return a 2x2 contingency table
     */
    private static int[][] getContingencyTable(Map<String, StratifiedAlignmentContext> stratifiedContexts, Allele ref, Allele alt) {
        int[][] table = new int[2][2];

        for ( Map.Entry<String, StratifiedAlignmentContext> sample : stratifiedContexts.entrySet() ) {
            for (PileupElement p : sample.getValue().getContext(StratifiedAlignmentContext.StratifiedContextType.COMPLETE).getBasePileup()) {
                if ( p.isDeletion() ) // ignore deletions
                    continue;

                if ( p.getRead().getMappingQuality() < 20 || p.getQual() < 20 )
                    continue; // todo -- fixme, should take filtered context!

                Allele base = Allele.create(p.getBase(), false);
                boolean isFW = !p.getRead().getReadNegativeStrandFlag();

                boolean matchesRef = ref.equals(base, true);
                boolean matchesAlt = alt.equals(base, true);
                if ( matchesRef || matchesAlt ) {
                    int row = matchesRef ? 0 : 1;
                    int column = isFW ? 0 : 1;

                    table[row][column]++;
                }
            }
        }

        return table;
    }
}

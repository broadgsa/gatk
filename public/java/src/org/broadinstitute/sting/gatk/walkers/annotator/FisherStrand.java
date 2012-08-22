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

import cern.jet.math.Arithmetic;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.ActiveRegionBasedAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.AnnotatorCompatible;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.InfoFieldAnnotation;
import org.broadinstitute.sting.gatk.walkers.annotator.interfaces.StandardAnnotation;
import org.broadinstitute.sting.gatk.walkers.genotyper.IndelGenotypeLikelihoodsCalculationModel;
import org.broadinstitute.sting.gatk.walkers.genotyper.PerReadAlleleLikelihoodMap;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLineType;
import org.broadinstitute.sting.utils.codecs.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.*;


/**
 * Phred-scaled p-value using Fisher's Exact Test to detect strand bias (the variation
 * being seen on only the forward or only the reverse strand) in the reads? More bias is
 * indicative of false positive calls.  Note that the fisher strand test may not be
 * calculated for certain complex indel cases or for multi-allelic sites.
 */
public class FisherStrand extends InfoFieldAnnotation implements StandardAnnotation, ActiveRegionBasedAnnotation {
    private static final String FS = "FS";
    private static final double MIN_PVALUE = 1E-320;
    public Map<String, Object> annotate(final RefMetaDataTracker tracker,
                                        final AnnotatorCompatible walker,
                                        final ReferenceContext ref,
                                        final Map<String, AlignmentContext> stratifiedContexts,
                                        final VariantContext vc,
                                        final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap) {
        if ( !vc.isVariant() )
            return null;

        int[][] table;

        if (vc.isSNP() && stratifiedContexts != null) {
            table = getSNPContingencyTable(stratifiedContexts, vc.getReference(), vc.getAltAlleleWithHighestAlleleCount());
        }
        else if (stratifiedPerReadAlleleLikelihoodMap != null) {
            // either SNP with no alignment context, or indels: per-read likelihood map needed
            table = getContingencyTable(stratifiedPerReadAlleleLikelihoodMap, vc.getReference(), vc.getAltAlleleWithHighestAlleleCount());
        }
        else
        // for non-snp variants, we  need per-read likelihoods.
        // for snps, we can get same result from simple pileup
            return null;

        if (table == null)
            return null;

        Double pvalue = Math.max(pValueForContingencyTable(table), MIN_PVALUE);
        if ( pvalue == null )
            return null;

        Map<String, Object> map = new HashMap<String, Object>();
        map.put(FS, String.format("%.3f", QualityUtils.phredScaleErrorRate(pvalue)));
        return map;
    }

    public List<String> getKeyNames() {
        return Arrays.asList(FS);
    }

    public List<VCFInfoHeaderLine> getDescriptions() {
        return Arrays.asList(
            new VCFInfoHeaderLine(FS, 1, VCFHeaderLineType.Float, "Phred-scaled p-value using Fisher's exact test to detect strand bias"));
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

        return (table[0][0] >= 0 && table[1][1] >= 0);
    }

    private static boolean unrotateTable(int[][] table) {
        table[0][0] += 1;
        table[1][0] -= 1;

        table[0][1] -= 1;
        table[1][1] += 1;

        return (table[0][1] >= 0 && table[1][0] >= 0);
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
    private static int[][] getContingencyTable( final Map<String, PerReadAlleleLikelihoodMap> stratifiedPerReadAlleleLikelihoodMap,
                                                final Allele ref, final Allele alt) {
        int[][] table = new int[2][2];

        for (PerReadAlleleLikelihoodMap maps : stratifiedPerReadAlleleLikelihoodMap.values() ) {
            for (Map.Entry<GATKSAMRecord,Map<Allele,Double>> el : maps.getLikelihoodReadMap().entrySet()) {
                if ( el.getKey().isReducedRead() ) // ignore reduced reads
                    continue;
                final boolean matchesRef = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue()).equals(ref,true);
                final boolean matchesAlt = PerReadAlleleLikelihoodMap.getMostLikelyAllele(el.getValue()).equals(alt,true);

                if ( !matchesRef && !matchesAlt )
                    continue;

                boolean isFW = el.getKey().getReadNegativeStrandFlag();

                int row = matchesRef ? 0 : 1;
                int column = isFW ? 0 : 1;

                table[row][column]++;
            }
        }

        return table;
    }
    /**
     Allocate and fill a 2x2 strand contingency table.  In the end, it'll look something like this:
     *             fw      rc
     *   allele1   #       #
     *   allele2   #       #
     * @return a 2x2 contingency table
     */
    private static int[][] getSNPContingencyTable(Map<String, AlignmentContext> stratifiedContexts, Allele ref, Allele alt) {
        int[][] table = new int[2][2];

        for ( Map.Entry<String, AlignmentContext> sample : stratifiedContexts.entrySet() ) {
            for (PileupElement p : sample.getValue().getBasePileup()) {
                if ( ! RankSumTest.isUsableBase(p, false) || p.getRead().isReducedRead() ) // ignore deletions and reduced reads
                    continue;

                Allele base = Allele.create(p.getBase(), false);
                boolean isFW = !p.getRead().getReadNegativeStrandFlag();

                final boolean matchesRef = ref.equals(base, true);
                final boolean matchesAlt = alt.equals(base, true);
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

/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.VariantEvalWalker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.VCFConstants;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.*;

import java.util.*;

@Analysis(description = "Evaluation summary for multi-allelic variants")
public class MultiallelicSummary extends VariantEvaluator { // implements StandardEval {
    final protected static Logger logger = Logger.getLogger(MultiallelicSummary.class);

    public enum Type {
        SNP, INDEL
    }

    // basic counts on various rates found
    @DataPoint(description = "Number of processed loci")
    public long nProcessedLoci = 0;

    @DataPoint(description = "Number of SNPs")
    public int nSNPs = 0;
    @DataPoint(description = "Number of multi-allelic SNPs")
    public int nMultiSNPs = 0;
    @DataPoint(description = "% processed sites that are multi-allelic SNPs", format = "%.5f")
    public double processedMultiSnpRatio = 0;
    @DataPoint(description = "% SNP sites that are multi-allelic", format = "%.3f")
    public double variantMultiSnpRatio = 0;

    @DataPoint(description = "Number of Indels")
    public int nIndels = 0;
    @DataPoint(description = "Number of multi-allelic Indels")
    public int nMultiIndels = 0;
    @DataPoint(description = "% processed sites that are multi-allelic Indels", format = "%.5f")
    public double processedMultiIndelRatio = 0;
    @DataPoint(description = "% Indel sites that are multi-allelic", format = "%.3f")
    public double variantMultiIndelRatio = 0;

    @DataPoint(description = "Number of Transitions")
    public int nTi = 0;
    @DataPoint(description = "Number of Transversions")
    public int nTv = 0;
    @DataPoint(description = "Overall TiTv ratio", format = "%.2f")
    public double TiTvRatio = 0;

    @DataPoint(description = "Multi-allelic SNPs partially known")
    public int knownSNPsPartial = 0;
    @DataPoint(description = "Multi-allelic SNPs completely known")
    public int knownSNPsComplete = 0;
    @DataPoint(description = "Multi-allelic SNP Novelty Rate")
    public String SNPNoveltyRate = "NA";

    //TODO -- implement me
    //@DataPoint(description = "Multi-allelic Indels partially known")
    public int knownIndelsPartial = 0;
    //@DataPoint(description = "Multi-allelic Indels completely known")
    public int knownIndelsComplete = 0;
    //@DataPoint(description = "Multi-allelic Indel Novelty Rate")
    public String indelNoveltyRate = "NA";

    @DataPoint(description="Histogram of allele frequencies for most common alternate allele")
    AFHistogram AFhistogramMax = new AFHistogram();

    @DataPoint(description="Histogram of allele frequencies for less common alternate alleles")
    AFHistogram AFhistogramMin = new AFHistogram();

    /*
     * AF histogram table object
     */
    static class AFHistogram implements TableType {
        private Object[] colKeys, rowKeys = {"pairwise_AF"};
        private int[] AFhistogram;

        private static final double AFincrement = 0.01;
        private static final int numBins = (int)(1.00 / AFincrement);

        public AFHistogram() {
            colKeys = initColKeys();
            AFhistogram = new int[colKeys.length];
        }

        public Object[] getColumnKeys() {
            return colKeys;
        }

        public Object[] getRowKeys() {
            return rowKeys;
        }

        public Object getCell(int row, int col) {
            return AFhistogram[col];
        }

        private static Object[] initColKeys() {
            ArrayList<String> keyList = new ArrayList<String>(numBins + 1);
            for ( double a = 0.00; a <= 1.01; a += AFincrement ) {
                keyList.add(String.format("%.2f", a));
            }
            return keyList.toArray();
        }

        public String getName() { return "AFHistTable"; }

        public void update(final double AF) {
            final int bin = (int)(numBins * MathUtils.round(AF, 2));
            AFhistogram[bin]++;
       }
    }

    public void initialize(VariantEvalWalker walker) {}

    @Override public boolean enabled() { return true; }

    public int getComparisonOrder() {
        return 2;
    }

    public void update0(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        nProcessedLoci += context.getSkippedBases() + (ref == null ? 0 : 1);
    }

    public String update2(VariantContext eval, VariantContext comp, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( eval == null || eval.isMonomorphicInSamples() )
            return null;

        // update counts
        switch ( eval.getType() ) {
            case SNP:
                nSNPs++;
                if ( !eval.isBiallelic() ) {
                    nMultiSNPs++;
                    calculatePairwiseTiTv(eval);
                    calculateSNPPairwiseNovelty(eval, comp);
                }
                break;
            case INDEL:
                nIndels++;
                if ( !eval.isBiallelic() ) {
                    nMultiIndels++;
                    calculateIndelPairwiseNovelty(eval, comp);
                }
                break;
            default:
                throw new UserException.BadInput("Unexpected variant context type: " + eval);
        }
        updateAFhistogram(eval);
        
        return null; // we don't capture any interesting sites
    }

    private void calculatePairwiseTiTv(VariantContext vc) {
        for ( Allele alt : vc.getAlternateAlleles() ) {
            if ( VariantContextUtils.isTransition(vc.getReference(), alt) )
                nTi++;
            else
                nTv++;
        }
    }

    private void calculateSNPPairwiseNovelty(VariantContext eval, VariantContext comp) {
        if ( comp == null )
            return;

        int knownAlleles = 0;
        for ( Allele alt : eval.getAlternateAlleles() ) {
            if ( comp.getAlternateAlleles().contains(alt) )
                knownAlleles++;
        }

        if ( knownAlleles == eval.getAlternateAlleles().size() )
            knownSNPsComplete++;
        else if ( knownAlleles > 0 )
            knownSNPsPartial++;
    }

    private void calculateIndelPairwiseNovelty(VariantContext eval, VariantContext comp) {
    }

    private void updateAFhistogram(VariantContext vc) {

        final Object obj = vc.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY, null);
        if ( obj == null || !(obj instanceof List) )
            return;

        List<String> list = (List<String>)obj;
        ArrayList<Double> AFs = new ArrayList<Double>(list.size());
        for ( String str : list ) {
            AFs.add(Double.valueOf(str));
        }

        Collections.sort(AFs);
        AFhistogramMax.update(AFs.get(AFs.size()-1));
        for ( int i = 0; i < AFs.size() - 1; i++ )
            AFhistogramMin.update(AFs.get(i));
    }
    
    private final String noveltyRate(final int all, final int known) {
        final int novel = all - known;
        final double rate = (novel / (1.0 * all));
        return all == 0 ? "NA" : String.format("%.2f", rate);
    }

    public void finalizeEvaluation() {
        processedMultiSnpRatio = (double)nMultiSNPs / (double)nProcessedLoci;
        variantMultiSnpRatio = (double)nMultiSNPs / (double)nSNPs;
        processedMultiIndelRatio = (double)nMultiIndels / (double)nProcessedLoci;
        variantMultiIndelRatio = (double)nMultiIndels / (double)nIndels;

        TiTvRatio = (double)nTi / (double)nTv;

        SNPNoveltyRate = noveltyRate(nMultiSNPs, knownSNPsPartial + knownSNPsComplete);
        indelNoveltyRate = noveltyRate(nMultiSNPs, knownIndelsPartial + knownIndelsComplete);
    }
}

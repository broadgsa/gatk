/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.varianteval.evaluators;

import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;

/**
 * Created by knoblett on 9/15/15.
 */

@Analysis(description = "Metrics Collection")
public class MetricsCollection extends VariantEvaluator {

    @DataPoint(description = "The concordance rate from CompOverlap", format = "%.2f")
    public double concordantRate;
    @DataPoint(description = "Number of SNPs from IndelSummary", format = "%d")
    public int nSNPs;
    @DataPoint(description = "Number of SNP loci from CountVariants", format = "%d")
    public long nSNPloci;
    @DataPoint(description = "Number of indels from IndelSummary", format = "%d")
    public int nIndels;
    @DataPoint(description = "Number of indel loci from MultiallelicSummary", format = "%d")
    public int nIndelLoci;
    @DataPoint(description = "Insertion  to deletion ratio from IndelSummary")
    public String indelRatio;
    @DataPoint(description = "Insertion to deletion ratio from CountVariants", format = "%.2f")
    public double indelRatioLociBased;
    @DataPoint(description = "The transition to transversion ratio from TiTvVariantEvaluator", format = "%.2f")
    public double tiTvRatio;

    public int getComparisonOrder() {return 2;}

    public void setData(double concordantRate, int nSNPs, long nSNPloci, int nIndels, int nIndelLoci, String indelRatio, double indelRatioLociBased, double tiTvRatio){
        this.concordantRate = concordantRate;
        this.nSNPs = nSNPs;
        this.nSNPloci = nSNPloci;
        this.nIndels = nIndels;
        this.nIndelLoci = nIndelLoci;
        this.indelRatio = indelRatio;
        this.indelRatioLociBased = indelRatioLociBased;
        this.tiTvRatio = tiTvRatio;
    }
}

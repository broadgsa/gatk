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

import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.Analysis;
import org.broadinstitute.gatk.tools.walkers.varianteval.util.DataPoint;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

@Analysis(description = "Computes different estimates of theta based on variant sites and genotypes")
public class ThetaVariantEvaluator extends VariantEvaluator {
    @DataPoint(description = "Average heterozygosity at variant sites; note that missing genotypes are ignored when computing this value", format = "%.8f")
    public double avgHet = 0.0;
    @DataPoint(description = "Average pairwise differences at aligned sequences; averaged over both number of sequeneces and number of variant sites; note that missing genotypes are ignored when computing this value", format = "%.8f")
    public double avgAvgDiffs = 0.0;
    @DataPoint(description = "Sum of heterozygosity over all variant sites; divide this by total target to get estimate of per base theta", format = "%.8f")
    public double totalHet = 0.0;
    @DataPoint(description = "Sum of pairwise diffs over all variant sites; divide this by total target to get estimate of per base theta", format = "%.8f")
    public double totalAvgDiffs = 0.0;
    @DataPoint(description = "Theta for entire region estimated based on number of segregating sites; divide ths by total target to get estimate of per base theta", format = "%.8f")
    public double thetaRegionNumSites = 0.0;

    //helper variables
    double numSites = 0;

    public int getComparisonOrder() {
        return 1;
    }

    public void update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (vc == null || !vc.isSNP() || (getWalker().ignoreAC0Sites() && vc.isMonomorphicInSamples())) {
            return;
        }

        //this maps allele to a count
        ConcurrentMap<String, Integer> alleleCounts = new ConcurrentHashMap<String, Integer>();

        int numHetsHere = 0;
        int numGenosHere = 0;
        int numIndsHere = 0;

        for (final Genotype genotype : vc.getGenotypes()) {
            numIndsHere++;
            if (!genotype.isNoCall()) {
                //increment stats for heterozygosity
                if (genotype.isHet()) {
                    numHetsHere++;
                }

                numGenosHere++;
                //increment stats for pairwise mismatches

                for (Allele allele : genotype.getAlleles()) {
                    if (allele.isCalled()) {
                        String alleleString = allele.toString();
                        alleleCounts.putIfAbsent(alleleString, 0);
                        alleleCounts.put(alleleString, alleleCounts.get(alleleString) + 1);
                    }
                }
            }
        }
        if (numGenosHere > 0) {
            //only if have one called genotype at least
            this.numSites++;

            this.totalHet += numHetsHere / (double)numGenosHere;

            //compute based on num sites
            float harmonicFactor = 0;
            for (int i = 1; i <= numIndsHere; i++) {
                harmonicFactor += 1.0 / i;
            }
            this.thetaRegionNumSites += 1.0 / harmonicFactor;

            //now compute pairwise mismatches
            float numPairwise = 0;
            int numDiffs = 0;
            for (String allele1 : alleleCounts.keySet()) {
                int allele1Count = alleleCounts.get(allele1);

                for (String allele2 : alleleCounts.keySet()) {
                    if (allele1.compareTo(allele2) < 0) {
                        continue;
                    }
                    if (allele1 .compareTo(allele2) == 0) {
                        numPairwise += allele1Count * (allele1Count - 1) * .5;

                    }
                    else {
                        int allele2Count = alleleCounts.get(allele2);
                        numPairwise += allele1Count * allele2Count;
                        numDiffs += allele1Count * allele2Count;
                    }
                }
            }

            if (numPairwise > 0) {
                this.totalAvgDiffs += numDiffs / numPairwise;
            }
        }
    }

    @Override
    public void finalizeEvaluation() {

        if (this.numSites > 0) {

            this.avgHet = this.totalHet / this.numSites;
            this.avgAvgDiffs = this.totalAvgDiffs / this.numSites;

        }
    }
}
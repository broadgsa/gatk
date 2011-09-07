package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

@Analysis(description = "Computes different estimates of theta based on variant sites and genotypes")
public class ThetaVariantEvaluator extends VariantEvaluator {
    @DataPoint(description = "Average heterozygosity at variant sites; note that missing genotypes are ignored when computing this value")
    double avgHet = 0.0;
    @DataPoint(description = "Average pairwise differences at aligned sequences; averaged over both number of sequeneces and number of variant sites; note that missing genotypes are ignored when computing this value")
    double avgAvgDiffs = 0.0;
    @DataPoint(description = "Sum of heterozygosity over all variant sites; divide this by total target to get estimate of per base theta")
    double totalHet = 0.0;
    @DataPoint(description = "Sum of pairwise diffs over all variant sites; divide this by total target to get estimate of per base theta")
    double totalAvgDiffs = 0.0;
    @DataPoint(description = "Theta for entire region estimated based on number of segregating sites; divide ths by total target to get estimate of per base theta")
    double thetaRegionNumSites = 0.0;

    //helper variables
    double numSites = 0;

    public boolean enabled() {
        return true;
    }

    public int getComparisonOrder() {
        return 1;
    }

    public String update1(VariantContext vc, RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if (vc == null || !vc.isSNP() || !vc.hasGenotypes() || vc.isMonomorphic()) {
            return null; //no interesting sites
        }

        //this maps allele to a count
        ConcurrentMap<String, Integer> alleleCounts = new ConcurrentHashMap<String, Integer>();

        int numHetsHere = 0;
        float numGenosHere = 0;
        int numIndsHere = 0;

        for (Genotype genotype : vc.getGenotypes().values()) {
            numIndsHere++;
            if (!genotype.isNoCall()) {
                //increment stats for heterozygosity
                if (genotype.isHet()) {
                    numHetsHere++;
                }

                numGenosHere++;
                //increment stats for pairwise mismatches

                for (Allele allele : genotype.getAlleles()) {
                    if (allele.isNonNull() && allele.isCalled()) {
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

            this.totalHet += numHetsHere / numGenosHere;

            //compute based on num sites
            float harmonicFactor = 0;
            for (int i = 1; i <= numIndsHere; i++) {
                harmonicFactor += 1.0 / i;
            }
            this.thetaRegionNumSites += 1.0 / harmonicFactor;

            //now compute pairwise mismatches
            float numPairwise = 0;
            float numDiffs = 0;
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

        return null;
    }

    @Override
    public void finalizeEvaluation() {

        if (this.numSites > 0) {

            this.avgHet = this.totalHet / this.numSites;
            this.avgAvgDiffs = this.totalAvgDiffs / this.numSites;

        }
    }
}
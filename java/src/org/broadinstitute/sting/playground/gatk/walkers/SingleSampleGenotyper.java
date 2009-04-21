package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;

import org.broadinstitute.sting.playground.utils.*;
import org.broadinstitute.sting.playground.gatk.*;
import org.broadinstitute.sting.playground.gatk.walkers.*;

import net.sf.samtools.SAMRecord;

import java.util.List;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<AlleleFrequencyEstimate, Integer> 
{
    AlleleMetrics metrics;
    
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        return true;    // We are keeping all the reads
    }

    public boolean requiresReads()     { return true; }    

    public void initialize()
    {
        metrics = new AlleleMetrics("metrics.out");
    }

    protected class GenotypeLikelihoods
    {
        public double[] likelihoods;
        public String[] genotypes;

        GenotypeLikelihoods()
        {
            likelihoods = new double[10];
            genotypes   = new String[10];

            genotypes[0] = "AA";
            genotypes[1] = "AC";
            genotypes[2] = "AG";
            genotypes[3] = "AT";
            genotypes[4] = "CC";
            genotypes[5] = "CG";
            genotypes[6] = "CT";
            genotypes[7] = "GG";
            genotypes[8] = "GT";
            genotypes[9] = "TT";
        }

        void add(char ref, char read, byte qual)
        {
            double p_error = Math.pow(10.0, (double)qual / -10);  
            for (int i = 0; i < genotypes.length; i++)
            {
                likelihoods[i] += AlleleLikelihood(ref, read, genotypes[i], p_error);
            }
        }

        double AlleleLikelihood(char ref, char read, String genotype, double p_error)
        {
            char h1 = genotype.charAt(0);
            char h2 = genotype.charAt(1);

            double p_base;

            if      ((h1 == h2) && (h1 == read))                 { p_base = Math.log10(1-p_error); }
            else if ((h1 != h2) && (h1 == read) || (h2 == read)) { p_base = Math.log10(0.5 - (p_error/2.0)); }
            else                                                 { p_base = Math.log10(p_error); }

            return p_base;
        }

        public String[] sorted_genotypes;
        public double[] sorted_likelihoods;

        public void sort()
        {
            Integer[] permutation = Utils.SortPermutation(likelihoods);

            Integer[] reverse_permutation = new Integer[permutation.length];
            for (int i = 0; i < reverse_permutation.length; i++) { reverse_permutation[i] = permutation[(permutation.length-1) - i]; }

            sorted_genotypes      = Utils.PermuteArray(genotypes, reverse_permutation);
            sorted_likelihoods    = Utils.PermuteArray(likelihoods, reverse_permutation);
        }

        public String toString(char ref)
        {
            this.sort();
            String s = String.format("%s %f %f ", this.BestGenotype(), this.LodVsNextBest(), this.LodVsRef(ref));
            for (int i = 0; i < sorted_genotypes.length; i++)
            {
                if (i != 0) { s = s + " "; }
                s = s + sorted_genotypes[i] + ":" + String.format("%.2f", sorted_likelihoods[i]);
            }
            return s;
        }

        public void ApplyPrior(char ref, double p_alt)
        {
            for (int i = 0; i < genotypes.length; i++)
            {
                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref))
                {
                    // hom-ref
                    likelihoods[i] += Math.log10(1.0 - 1e-3);
                }
                else if ((genotypes[i].charAt(0) != ref) && (genotypes[i].charAt(1) != ref))
                {
                    // hom-nonref
                    likelihoods[i] += Math.log10(1e-5);
                }
                else
                {
                    // het
                    likelihoods[i] += Math.log10(1e-3);
                }
            }
            this.sort();
        }

        public double LodVsNextBest()
        {
            this.sort();
            return sorted_likelihoods[0] - sorted_likelihoods[1];
        }

        double ref_likelihood = Double.NaN;
        public double LodVsRef(char ref)
        {
            if ((this.BestGenotype().charAt(0) == ref) && (this.BestGenotype().charAt(1) == ref)) { ref_likelihood = sorted_likelihoods[0]; return LodVsNextBest(); }
            else
            {
	            for (int i = 0; i < genotypes.length; i++)
	            {
	                if ((genotypes[i].charAt(0) == ref) && (genotypes[i].charAt(1) == ref))
	                {
	                    ref_likelihood = likelihoods[i];
	                }
	            }
            }
            return sorted_likelihoods[0] - ref_likelihood;
        }

        public String BestGenotype()
        {
            this.sort();
            return this.sorted_genotypes[0];
        }

        public double BestPosterior()
        {
            this.sort();
            return this.sorted_likelihoods[0];
        }

        public AlleleFrequencyEstimate toAlleleFrequencyEstimate(GenomeLoc location, char ref, int depth, String bases, double[] posteriors)
        {
            double qhat  = Double.NaN;
            double qstar = Double.NaN;
            char   alt   = 'N';

                if ((sorted_genotypes[0].charAt(0) == ref) && (sorted_genotypes[0].charAt(1) == ref))
                {
                    // hom-ref
                    qhat  = 0.0;
                    qstar = 0.0; 
                    alt   = 'N';
                }
                else if ((sorted_genotypes[0].charAt(0) != ref) && (sorted_genotypes[0].charAt(1) != ref))
                {
                    // hom-nonref
                    likelihoods[0] += Math.log10(1e-5);
                    qhat = 1.0;
                    qstar = 1.0;
                    alt = sorted_genotypes[0].charAt(0);
                }
                else
                {
                    // het
                    likelihoods[0] += Math.log10(1e-3);
                    qhat  = 0.5;
                    qstar  = 0.5;

                    if (sorted_genotypes[0].charAt(0) != ref) { alt = sorted_genotypes[0].charAt(0); }
                    if (sorted_genotypes[0].charAt(1) != ref) { alt = sorted_genotypes[0].charAt(1); }
                }
                
            return new AlleleFrequencyEstimate(location, ref, alt, 2, qhat, qstar, this.LodVsRef(ref), this.LodVsNextBest(), sorted_likelihoods[0], ref_likelihood, depth, bases, (double[][])null, this.likelihoods); 
        }

    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) 
    {
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        String bases = "";
        String quals = "";

        ref = Character.toUpperCase(ref);

        String rodString = "";
        // Look up dbsnp priors
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() ) 
        {
            if ( datum != null ) 
            {
                if ( datum instanceof rodDbSNP)
                {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    rodString += dbsnp.toString();
                }
                else 
                {
                    rodString += datum.toSimpleString();
                }
            }
        }
        if ( rodString != "" )
            rodString = "[ROD: " + rodString + "]";

        // Accumulate genotype likelihoods
        GenotypeLikelihoods G = new GenotypeLikelihoods(); 
        for ( int i = 0; i < reads.size(); i++ ) 
        {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            bases += read.getReadString().charAt(offset);
            quals += read.getBaseQualityString().charAt(offset);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }
        G.ApplyPrior(ref, Double.NaN);

        System.out.printf("%s %s %s %s\n", context.getLocation(), ref, bases, G.toString(ref), rodString);

        AlleleFrequencyEstimate freq = G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods);

        metrics.nextPosition(freq, tracker);
        metrics.printMetricsAtLocusIntervals(1000);

        return freq;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(AlleleFrequencyEstimate value, Integer sum) 
    {
        return 0;
    }
}

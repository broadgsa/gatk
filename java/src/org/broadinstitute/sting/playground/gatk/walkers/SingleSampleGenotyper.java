package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.Utils;
import net.sf.samtools.SAMRecord;

import java.util.List;

// Draft single sample genotyper
// j.maguire 3-7-2009

public class SingleSampleGenotyper extends LocusWalker<Integer, Integer> {
    public boolean filter(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public boolean requiresReads()     { return true; }    

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

        public String toString()
        {
            Integer[] permutation       = Utils.SortPermutation(likelihoods);
            String[] sorted_genotypes   = Utils.PermuteArray(genotypes, permutation);
            double[] sorted_likelihoods = Utils.PermuteArray(likelihoods, permutation);

            String s = "";
            for (int i = sorted_genotypes.length-1; i >= 0; i--)
            {
                if (i != sorted_genotypes.length-1) { s = s + " "; }
                s = s + sorted_genotypes[i] + ":" + sorted_likelihoods[i];
            }
            return s;
        }

    }

    // Map over the org.broadinstitute.sting.gatk.LocusContext
    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {
        //System.out.printf("Reads %s:%d %d%n", context.getContig(), context.getPosition(), context.getReads().size());
        //for ( SAMRecord read : context.getReads() ) {
        //    System.out.println("  -> " + read.getReadName());
        //}

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        String bases = "";
        String quals = "";
        //String offsetString = "";

        // Look up hapmap and dbsnp priors
        String rodString = "";
        for ( ReferenceOrderedDatum datum : rodData ) 
        {
            if ( datum != null ) 
            {
                if ( datum instanceof rodDbSNP)
                {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    rodString += dbsnp.toMediumString();
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

        if ( context.getLocation().getStart() % 1 == 0 ) {
            //System.out.printf("%s: %s %s %s %s%n", context.getLocation(), ref, bases, quals, rodString);
            out.printf("%s %s %s %s\n", ref, bases, G.toString(), rodString);
        }

        return 1;
    }

    // Given result of map function
    public Integer reduceInit() { return 0; }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }
}

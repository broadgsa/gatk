package org.broadinstitute.sting.gatk.walkers;
//import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.utils.AlleleFrequencyEstimate;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.Arrays;

public class AlleleFrequencyWalker extends BasicLociWalker<AlleleFrequencyEstimate, Integer> {

    public AlleleFrequencyEstimate map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {

        // Convert context data into bases and 4-base quals

        // Set number of chromosomes, N, to 2 for now
        int N = 2;

        // Convert bases to CharArray
        int numReads = context.getReads().size(); //numReads();
        byte[] bases = new byte[numReads];
        String base_string = "";
        double[][] quals = new double[numReads][4];
        int refnum = nuc2num[ref];

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        for (int i =0; i < numReads; i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            bases[i] = read.getReadBases()[offset];
            base_string += read.getReadString().charAt(offset);
            int callednum = nuc2num[bases[i]];
            quals[i][callednum] = 1 - Math.pow(10.0, (double)read.getBaseQualities()[offset] / -10.0);

            // Set all nonref qual scores to their share of the remaining probality not "used" by the reference base's qual
            double nonref_quals = (1.0 - quals[i][callednum]) / 3;
            for (int b=0; b<4; b++)
                if (b != callednum)
                    quals[i][b] = nonref_quals;
        }

        // Count bases
        int[] base_counts = new int[4];
        for (byte b : bases)
            base_counts[nuc2num[b]]++;
        
        // Find alternate allele - 2nd most frequent non-ref allele
        // (maybe we should check for ties and eval both or check most common including quality scores)
        int altnum = -1; // start with first base numerical identity (0,1,2,3)
        double altcount = -1; // start with first base count

        for (int b=0; b<4; b++) {
            if ((b != refnum) && (base_counts[b] > altcount)) {
                altnum = b; altcount = base_counts[b]; 
            }
        }
        assert(altnum != -1);

        AlleleFrequencyEstimate alleleFreq = AlleleFrequencyEstimator(context.getLocation().toString(), N, bases, quals, refnum, altnum, base_string.length());

        // Print dbSNP data if its there
        if (false) {
            for ( ReferenceOrderedDatum datum : rodData ) {
                if ( datum != null && datum instanceof rodDbSNP) {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    //System.out.printf("  DBSNP %s on %s => %s%n", dbsnp.toSimpleString(), dbsnp.strand, Utils.join("/", dbsnp.getAllelesFWD()));
                    System.out.printf("ROD: %s ",dbsnp.toMediumString());
                }
            }
        }
        return alleleFreq;
    }

    public AlleleFrequencyEstimate AlleleFrequencyEstimator(String location, int N, byte[] bases, double[][] quals, int refnum, int altnum, int depth)
    {

        // q = hypothetical %nonref
        // qstar = true %nonref
        // G = N, qstar, alt
        // N = number of chromosomes
        // alt = alternate hypothesis base (most frequent nonref base)
        // ref = reference base

        // b = number of bases at locus

        double epsilon = 0; //  1e-2;
        double qstar;
        int qstar_N;

        double qstep = 0.01;
        double qend = 1.0 + qstep / 10; // makes sure we get to 1.0 even with rounding error of qsteps accumulating

        double best_qstar     = Math.log10(0);
        double best_qhat      = Math.log10(0);
        double best_posterior = Math.log10(0);

        double best_pDq = Math.log10(0);
        double best_pqG = Math.log10(0);
        double best_pG  = Math.log10(0);

        for (double q=0.0; q <= qend; q += qstep) // hypothetic allele balance that we sample over
        {
            long q_R = Math.round(q*bases.length);
            for (qstar = epsilon + ((1.0 - 2*epsilon)/N), qstar_N = 1; qstar <= 1.0; qstar += (1.0 - 2*epsilon)/N, qstar_N++) // qstar - true allele balances
            { 
                // for N=2: these are 0.0 + epsilon, 0.5, 1.0 - epsilon corresponding to reference, het-SNP, homo-SNP
                double pDq = P_D_q(bases, quals, q, refnum, altnum);
                double pqG = P_q_G(bases, N, q, qstar, q_R);
                double pG  = P_G(N, qstar_N); 
                double posterior = pDq + pqG + pG;

                if (posterior > best_posterior) 
                {
                    best_qstar = qstar;
                    best_qhat  = q;
                    best_posterior = posterior; 

                    best_pDq = pDq;
                    best_pqG = pqG;
                    best_pG  = pG;
                }
            }
        }

        double posterior_null_hyp = P_D_q(bases, quals, 0.0, refnum, altnum) + P_q_G(bases, N, 0.0, epsilon, 0) + P_G(N, 0);
        double LOD = best_posterior - posterior_null_hyp;

        AlleleFrequencyEstimate alleleFreq = new AlleleFrequencyEstimate(location,
                                                                         num2nuc[refnum], 
                                                                         num2nuc[altnum],
                                                                         N, 
                                                                         best_qhat, 
                                                                         best_qstar, 
                                                                         LOD,
                                                                         depth);
        return alleleFreq;
    }

    public class MixtureLikelihood implements Comparable<MixtureLikelihood> {
        public double posterior;
        public double qstar;
        public double qhat;

        MixtureLikelihood(double posterior, double qstar, double qhat) {
            this.posterior = posterior;
            this.qstar = qstar;
            this.qhat = qhat;
        }

        public int compareTo (MixtureLikelihood other) {
            // default sorting is REVERSE sorting on posterior prob; gives array with top post. prob. as first element
            if (posterior < other.posterior)      return 1;
            else if (posterior > other.posterior) return -1;
            else                                  return 0;
        }
    }

    static double P_D_q(byte[] bases, double[][]quals, double q, int refnum, int altnum) 
    {
        double p = 0.0;

        for (int i=0; i<bases.length; i++) 
        {
            p += Math.log10((1-q) * quals[i][refnum] + q * quals[i][altnum]);
        }

        return p;
    }

    static double P_q_G(byte [] bases, int N, double q, double qstar, long q_R) 
    {
        return Math.log10(binomialProb(q_R, bases.length, qstar));
    }

    static double P_G(int N, int qstar_N) 
    {
        // badly hard coded right now.
        if      (qstar_N == 0) { return Math.log10(0.999); }
        else if (qstar_N == N) { return Math.log10(1e-5);  }
        else                   { return Math.log10(1e-3);  }
    }

    static String genotypeTypeString(double q, int N){
        switch (Math.round((float)q*N)) {
            case (0):
                return "ref";
            case (1):
                return "het";
            case (2):
                return "hom";
        }
        return "ERROR";
    }

    static double binomialProb(long k, long n, double p) {
        // k - numebr of successes
        // n - number of Bernoulli trials
        // p - probability of success

        return (double)nchoosek(n, k) * Math.pow(p, k) * Math.pow(1-p, n-k);
    }

    static long nchoosek(long n, long k) {
        long m = n - k;
        if (k < m)
            k = m;

        long t = 1;
        for (long i = n, j = 1; i > k; i--, j++)
            t = t * i / j;
        return t;
    }

    public static String repeat(char letter, long count) {
        // Repeat a character count times
        String result = "";
        for (int i=0; i<count; i++) {
            result = result + letter;
        }
        return result;
    }


    public Integer reduceInit() { return 0; }

    public Integer reduce(AlleleFrequencyEstimate alleleFreq, Integer sum) 
    {
        if (alleleFreq.LOD >= 5)
        {
	        System.out.print(String.format("RESULT %s %c %c %f %f %f %d\n", 
	                                        alleleFreq.location,
	                                        alleleFreq.ref, 
	                                        alleleFreq.alt, 
	                                        alleleFreq.qhat, 
	                                        alleleFreq.qstar, 
	                                        alleleFreq.LOD, 
	                                        alleleFreq.depth));
        }
        return 0;
    }


    static int nuc2num[];
    static char num2nuc[];
    public AlleleFrequencyWalker() {
        nuc2num = new int[128];
        nuc2num['A'] = 0;
        nuc2num['C'] = 1;
        nuc2num['T'] = 2;
        nuc2num['G'] = 3;
        nuc2num['a'] = 0;
        nuc2num['c'] = 1;
        nuc2num['t'] = 2;
        nuc2num['g'] = 3;

        num2nuc = new char[4];
        num2nuc[0] = 'A';
        num2nuc[1] = 'C';
        num2nuc[2] = 'T';
        num2nuc[3] = 'G';
    }



    void print_base_qual_matrix(double [][]quals, int numReads) {
        // Print quals for debugging
        System.out.println("Base quality matrix");
        for (int b=0; b<4; b++) {
            System.out.print(num2nuc[b]);
            for (int i=0; i < numReads; i++){
                System.out.format(" %.4f", quals[i][b]);
            }
            System.out.println();
        }
    }

    public static void main(String[] args)
    {
        int N = 2;

        byte[] het_bases = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  
                        1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
        double[][] het_quals = {{0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.999, 0.001/3.0, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0},
                                {0.001/3.0, 0.999, 0.001/3.0, 0.001/3.0}};

        AlleleFrequencyWalker w = new AlleleFrequencyWalker();
        AlleleFrequencyEstimate estimate = w.AlleleFrequencyEstimator("null", N, het_bases, het_quals, 0, 1, 20);

        System.out.print(String.format("50/50 Het : %s %c %c %f %f %f %d %s\n", 
                                        "null", estimate.ref, estimate.alt, estimate.qhat, estimate.qstar, estimate.LOD, 20, "null"));

    }

}

















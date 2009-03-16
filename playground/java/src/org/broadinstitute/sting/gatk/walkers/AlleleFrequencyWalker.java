package org.broadinstitute.sting.gatk.walkers;

//import org.broadinstitute.sting.gatk.iterators.LocusIterator;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import net.sf.samtools.SAMRecord;

import java.util.List;

public class AlleleFrequencyWalker extends BasicLociWalker<Integer, Integer> {


    

    public Integer map(List<ReferenceOrderedDatum> rodData, char ref, LocusContext context) {

        // Convert context data into bases and 4-base quals

        // Set number of chromosomes, N, to 2 for now
        int N = 2;
        boolean debug = false;

        // Convert bases to CharArray
        int numReads = context.getReads().size(); //numReads();
        byte[] bases = new byte[numReads];
        double[][] quals = new double[numReads][4];
        int refnum = nuc2num[ref];

        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        for (int i =0; i < numReads; i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            bases[i] = read.getReadBases()[offset];
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
        assert(altnum > 0);

        if (bases.length > 0) { //altcount > 2) {
            System.out.format("Pos: %s | ref: %c | alt: %c | %2dA | %2dC | %2dT | %2dG | %2d total | ", context.getLocation(), num2nuc[refnum],
                    num2nuc[altnum], base_counts[0], base_counts[1], base_counts[2], base_counts[3], bases.length);

            if (debug) print_base_qual_matrix(quals, numReads);

            // Check if we really want to do this one
            AlleleFrequencyEstimator(N, bases, quals, refnum, altnum, debug);
        }
        if (debug) System.out.println();

        return 1;
    }

    static void AlleleFrequencyEstimator(int N, byte[] bases, double[][] quals, int refnum, int altnum, boolean debug) {

        // q = hypothetical %nonref
        // qstar = true %nonref
        // G = N, qstar, alt
        // N = number of chromosomes
        // alt = alternate hypothesis base (most frequent nonref base)
        // ref = reference base

        double epsilon = 0; //  1e-2;
        double qstar;
        if (debug) {
            System.out.format("%4s | ", "q ");
            for (qstar = epsilon; qstar <= 1.0; qstar += (1.0 - 2*epsilon)/N)
                System.out.format("%5s%12s %12s %12s | ", "qstar", "p(D|q)   ", "p(q|G)   ", "posterior  ");
            if (debug) System.out.println();
        }

        double highest_posterior = -1.0;
        double highest_qstar = -1.0;
        double highest_q = -1.0;

        double qstep = 0.01;
        double qend = 1.0 + qstep / 10; // makes sure we get to 1.0 even with rounding error of qsteps accumulating
        for (double q=0.0; q <= qend; q += qstep) {
            if (debug) System.out.format("%.2f | ", q);
            for (qstar = epsilon; qstar <= 1.0; qstar += (1.0 - 2*epsilon)/N) {
                double pDq = P_D_q(bases, quals, q, refnum, altnum);
                double pqG = P_q_G(bases, N, q, qstar);
                double pG = P_G(N, qstar, altnum);
                double posterior = pDq * pqG * pG;
                if (debug) System.out.format("%.2f %.10f %.10f %.10f | ", qstar, pDq, pqG, posterior);
                if (posterior > highest_posterior) {
                    highest_posterior = posterior;
                    highest_qstar = qstar;
                    highest_q = q;
                }
            }
            if (debug) System.out.println();
        }
        /*System.out.printf("Maximal likelihood posterior: %.6f\n", highest_posterior);
        System.out.printf("q that maximimizes posterior: %.6f\n", highest_q);
        System.out.printf("qstar used when calculating max q: %.6f\n", highest_qstar);*/
        System.out.printf("qhat %.2f | qstar %.2f | ", highest_q, highest_qstar);

        // Print out the called bases
        long numNonrefBases = Math.round(highest_q * N);
        long numRefBases = N - numNonrefBases;
        String genotype = repeat(num2nuc[refnum], numRefBases) + repeat(num2nuc[altnum], numNonrefBases);
        System.out.printf("gen: %s\n", genotype);
    }

    static double P_D_q(byte[] bases, double[][]quals, double q, int refnum, int altnum) {
        double p = 1.0;

        for (int i=0; i<bases.length; i++) {
            p *= (1-q) * quals[i][refnum] + q * quals[i][altnum];
        }

        return p;
    }

    static double P_q_G(byte [] bases, int N, double q, double qstar) {
        return binomialProb(Math.round(q*bases.length), bases.length, qstar);
    }

    static double P_G(int N, double qstar, int altnum) {
        if (N==2) {
            return p_G_N_2[ Math.round((float)(qstar * N)) ];
        }else{
            return 1.0;
        }
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
        String result = "";
        for (int i=0; i<count; i++) {
            result = result + letter;
        }
        return result;
    }


    public Integer reduceInit() { return 0; }

    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }


    static int nuc2num[];
    static char num2nuc[];
    static double p_G_N_2[]; // pop. gen. priors for N=2
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

        p_G_N_2 = new double[3];
        p_G_N_2[0] = 0.999;
        p_G_N_2[1] = 1e-3;
        p_G_N_2[2] = 1e-5;
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

}

















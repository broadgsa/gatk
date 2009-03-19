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
        int debug = 0;

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
        assert(altnum > 0);

        if (debug >= 1) System.out.format("Pos: %s ref: %c alt: %c ", context.getLocation(), num2nuc[refnum], num2nuc[altnum]);
        //System.out.format("%2dA %2dC %2dT %2dG %2d total", base_counts[0], base_counts[1], base_counts[2], base_counts[3], bases.length);

        if (debug >= 2) print_base_qual_matrix(quals, numReads);

        AlleleFrequencyEstimate alleleFreq = AlleleFrequencyEstimator(N, bases, quals, refnum, altnum, debug);

        // Print dbSNP data if its there
        if (debug >= 1) {
            for ( ReferenceOrderedDatum datum : rodData ) {
                if ( datum != null && datum instanceof rodDbSNP) {
                    rodDbSNP dbsnp = (rodDbSNP)datum;
                    //System.out.printf("  DBSNP %s on %s => %s%n", dbsnp.toSimpleString(), dbsnp.strand, Utils.join("/", dbsnp.getAllelesFWD()));
                    System.out.printf("ROD: %s ",dbsnp.toMediumString());
                }
            }
        }
        if (debug >= 1) System.out.format(" %s\n", base_string);
        if (debug >= 2) System.out.println();

        return alleleFreq;
    }

    AlleleFrequencyEstimate AlleleFrequencyEstimator(int N, byte[] bases, double[][] quals, int refnum, int altnum, int debug) {

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
        if (debug >= 2) {
            System.out.format("%4s ", "q ");
            for (qstar = epsilon; qstar <= 1.0; qstar += (1.0 - 2*epsilon)/N)
                System.out.format("%5s%10s %10s %10s ", "qstar", "p(D|q)   ", "p(q|G)   ", "posterior  ");
            System.out.println();
        }

        MixtureLikelihood[] bestMixtures = { new MixtureLikelihood(-1,-1,-1), new MixtureLikelihood(-1,-1,-1), new MixtureLikelihood(-1,-1,-1) };
        /*double[] highest_posterior = new double[] {-1.0, -1.0, -1.0};
        double[] highest_qstar = new double[] {-1.0, -1.0, -1.0};
        double[] highest_q = new double[] {-1.0, -1.0, -1.0};*/

        double qstep = 0.01;
        double qend = 1.0 + qstep / 10; // makes sure we get to 1.0 even with rounding error of qsteps accumulating
        for (double q=0.0; q <= qend; q += qstep) { // hypothetic allele balance that we sample over
            if (debug >= 2) System.out.format("%.2f ", q);

            long q_R = Math.round(q*bases.length);
            for (qstar = epsilon, qstar_N = 0; qstar <= 1.0; qstar += (1.0 - 2*epsilon)/N, qstar_N++) { // qstar - true allele balances
                // for N=2: these are 0.0 + epsilon, 0.5, 1.0 - epsilon corresponding to reference, het-SNP, homo-SNP
                //int qstar_N = Math.round((float)qstar*N);
                double pDq = P_D_q(bases, quals, q, refnum, altnum);
                double pqG = P_q_G(bases, N, q, qstar, q_R);
                double pG = P_G(N, qstar_N); //= P_G(N, qstar);
                double posterior = pDq * pqG * pG;
                if (debug >= 2) System.out.format("%.2f %10.7f %10.7f %10.7f ", qstar, Math.log10(pDq), Math.log10(pqG), Math.log10(posterior));
                if (posterior > bestMixtures[qstar_N].posterior)
                    bestMixtures[qstar_N] = new MixtureLikelihood(posterior, qstar, q);
            }
            if (debug >= 2) System.out.println();
        }

        // Calculate log-odds difference - must happen before mixtures are sorted by posterior prob.
        double logOddsVarRef = logOddsNonrefRef(bestMixtures);
        if (debug >= 1) System.out.printf("LOD(Var,Ref): %6.2f ", logOddsVarRef);

        // Print all mixture likelihoods in order of qstar
        for (MixtureLikelihood mix : bestMixtures) // outputs: "het 3.26 "
            if (debug >= 1) System.out.printf("%3s: %6.2f ", genotypeTypeString(mix.qstar, N), Math.log10(mix.posterior / (1-mix.posterior)));

        // Sort by posterior probability
        Arrays.sort(bestMixtures); // We only used best entry, so this could be a generic max function (optimization)
        if (debug >= 1) System.out.printf("qhat: %.2f qstar: %.2f ", bestMixtures[0].qhat, bestMixtures[0].qstar); //highest_q, highest_qstar);
        if (debug >= 1) System.out.printf("prob: %.2e ", bestMixtures[0].posterior);
        if (debug >= 1) System.out.printf("type: %3s ", genotypeTypeString(bestMixtures[0].qstar, N));

        double posterior_null_hyp = P_D_q(bases, quals, 0.0, refnum, altnum) * P_q_G(bases, N, 0.0, 0.0, 0) * P_G(N, 0);
        if (debug >= 1) System.out.printf("P(null): %.2e ", posterior_null_hyp);
        //double LOD = Math.log10(bestMixtures[0].posterior) - Math.log10(posterior_null_hyp);
        double LOD = logOdds(bestMixtures[0].posterior, posterior_null_hyp);
        if (debug >= 1) System.out.printf("LOD: %6.2f ", LOD);

        AlleleFrequencyEstimate alleleFreq = new AlleleFrequencyEstimate(num2nuc[refnum], num2nuc[refnum],
                N, bestMixtures[0].qhat, bestMixtures[0].qstar, logOddsVarRef);
        if (debug >= 1) System.out.printf("gen: %s", genotypeTypeString(bestMixtures[0].qstar, N));

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

    static double logOdds(double prob1, double prob2)  {
        return Math.log10(prob1 / (1-prob1)) - Math.log10(prob2 / (1-prob2));
    }

    static double logOddsNonrefRef(MixtureLikelihood[] mixtures) {
        // takes list of mixtureLikelihoods SORTED BY QSTAR (default ordering)
        // Works for N == 2 right now!
        assert (mixtures.length == 2);
        MixtureLikelihood bestNonref = ((mixtures[1].posterior > mixtures[2].posterior) ? mixtures[1] : mixtures[2]);
        return logOdds(bestNonref.posterior, mixtures[0].posterior);

        // Code to calculate LOD using qstar=0 qhat=0
        //double posterior_null_hyp = P_D_q(bases, quals, 0.0, refnum, altnum) * P_q_G(bases, N, 0.0, 0.0) * P_G(N, 0.0);
        //double LOD = Math.log10(highest_posterior) - Math.log10(posterior_null_hyp);
    }

    static double P_D_q(byte[] bases, double[][]quals, double q, int refnum, int altnum) {
        double p = 1.0;

        for (int i=0; i<bases.length; i++) {
            p *= (1-q) * quals[i][refnum] + q * quals[i][altnum];
        }

        return p;
    }

    static double P_q_G(byte [] bases, int N, double q, double qstar, long q_R) {
        return binomialProb(q_R, bases.length, qstar);
    }

    static double P_G(int N, int qstar_N) {
        if (N==2) {
            return p_G_N_2[ qstar_N ];
        }else{
            return 1.0;
        }
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

    public Integer reduce(AlleleFrequencyEstimate alleleFreq, Integer sum) {

        //System.out.printf("%s %.2f\n", alleleFreq.asString(), alleleFreq.logOddsVarRef);
        return 0;//value + sum;
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

















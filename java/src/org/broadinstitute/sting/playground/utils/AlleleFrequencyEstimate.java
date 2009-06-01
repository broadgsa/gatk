package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import java.util.Arrays;
import java.lang.Math;
import org.broadinstitute.sting.utils.GenomeLoc;

public class AlleleFrequencyEstimate {

        static 
        {
	        boolean assertsEnabled = false;
	        assert assertsEnabled = true; // Intentional side effect!!!
	        if (!assertsEnabled)
            {
                System.err.printf("\n\n\nERROR: You must run with asserts enabled. \"java -ea\".\n\n\n");
	            throw new RuntimeException("Asserts must be enabled!");
            }
        }

    //AlleleFrequencyEstimate();
    public GenomeLoc location;
    public char ref;
    public char alt;
    public int N;
    public double qhat;
    public double qstar;
    public double lodVsRef;
    public double lodVsNextBest;
    public double pBest;
    public double pRef;
    public int depth;
    public String notes;
    public String bases;
    public double[][] quals;
    public double[] posteriors;
	public String sample_name;

    GenomeLoc l;


    public AlleleFrequencyEstimate(GenomeLoc location, char ref, char alt, int N, double qhat, double qstar, double lodVsRef, double lodVsNextBest, double pBest, double pRef, int depth, String bases, double[][] quals, double[] posteriors, String sample_name)
    {
        if( Double.isNaN(lodVsRef)) { System.out.printf("%s: lodVsRef is NaN\n", location.toString()); }
        if( Double.isNaN(lodVsNextBest)) { System.out.printf("lodVsNextBest is NaN\n"); }
        if( Double.isNaN(qhat)) { System.out.printf("qhat is NaN\n"); }
        if( Double.isNaN(qstar)) { System.out.printf("qstar is NaN\n"); }
        if( Double.isNaN(pBest)) { System.out.printf("pBest is NaN\n"); }
        if( Double.isNaN(pRef)) { System.out.printf("pRef is NaN (%c %s)\n", ref, bases); }

        if( Double.isInfinite(lodVsRef)) 
        { 
            System.out.printf("lodVsRef is Infinite: %s %c %s\n", location.toString(), ref, bases); 
            for (int i = 0; i < posteriors.length; i++)
            {
                System.out.printf("POSTERIOR %d %f\n", i, posteriors[i]);
            }
        }
        if( Double.isInfinite(lodVsNextBest)) { System.out.printf("lodVsNextBest is Infinite\n"); }
        if( Double.isInfinite(qhat)) { System.out.printf("qhat is Infinite\n"); }
        if( Double.isInfinite(qstar)) { System.out.printf("qstar is Infinite\n"); }
        if( Double.isInfinite(pBest)) { System.out.printf("pBest is Infinite\n"); }
        if( Double.isInfinite(pRef)) { System.out.printf("pRef is Infinite\n"); }

        assert(! Double.isNaN(lodVsRef));
        assert(! Double.isNaN(lodVsNextBest));
        assert(! Double.isNaN(qhat));
        assert(! Double.isNaN(qstar));
        assert(! Double.isNaN(pBest));
        assert(! Double.isNaN(pRef));

        assert(! Double.isInfinite(lodVsRef));
        assert(! Double.isInfinite(lodVsNextBest));
        assert(! Double.isInfinite(qhat));
        assert(! Double.isInfinite(qstar));
        assert(! Double.isInfinite(pBest));
        assert(! Double.isInfinite(pRef));

        this.location = location;
        this.ref = ref;
        this.alt = alt;
        this.N = N;
        this.qhat = qhat;
        this.qstar = qstar;
        this.lodVsRef = lodVsRef;
        this.lodVsNextBest = lodVsNextBest;
		this.pBest = pBest;
		this.pRef = pRef;
        this.depth = depth;
        this.notes = "";
        this.bases = bases;
        this.quals = quals;
        this.posteriors = posteriors;
		this.sample_name = sample_name;
    }

    public boolean isREF() { return (this.lodVsRef <= -5.0); }
    public boolean isSNP() { return (this.lodVsRef >=  5.0); }

    /** Return the most likely genotype. */
    public String genotype()
    {
        int alt_count = (int)(qstar * N);
        int ref_count = N-alt_count;
        char[] alleles = new char[N];
        int i;
        for (i = 0; i < ref_count; i++) { alleles[i] = ref; }
        for (; i < N; i++) { alleles[i] = alt; }
        Arrays.sort(alleles);
        return new String(alleles);
    }

    public double emperical_allele_frequency()
    {
        return (double)Math.round((double)qstar * (double)N) / (double)N;
    }

    public double emperical_allele_frequency(int N)
    {
        return (double)Math.round((double)qstar * (double)N) / (double)N;
    }

    public String asGFFString()
    {
        String s = "";
        s += String.format("%s\tCALLER\tVARIANT\t%s\t%s\t%f\t.\t.\t",
                               location.getContig(),
                               location.getStart(),
                               location.getStart(),
                               lodVsRef);
        s += String.format("\t;\tSAMPLE %s", sample_name);
        s += String.format("\t;\tREF %c", ref);
        s += String.format("\t;\tALT %c", alt);
        s += String.format("\t;\tFREQ %f", qstar);
        s += String.format("\t;\tGENOTYPE %s", this.genotype());
        s += String.format("\t;\tDEPTH %d", depth);
        s += String.format("\t;\tLODvsREF %f", lodVsRef);
        s += String.format("\t;\tLODvsNEXTBEST %f", lodVsNextBest);
        s += String.format("\t;\tQHAT %f", qhat);
        s += String.format("\t;\tQSTAR %f", qstar);
        s += String.format("\t;\tBASES %s", bases);

        s += ";\n";

        // add quals.

        return s;
    }

    public static String asTabularStringHeader() 
	{
		return "location sample_name ref alt genotype qhat qstar lodVsRef lodVsNextBest depth bases";
	}

    public static String geliHeaderString() {
        return "#Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod dbSNP   AA      AC      AG      AT      CC      CG      CT      GG      GT      TT";
    }

    public String asGeliString() 
	{
        // #Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod dbSNP   AA      AC      AG      AT      CC      CG      CT      GG      GT      TT
        // chr1    7764136 A       48      99      CC      83.650421       9.18159         -92.83638       -18.367548      -96.91729       -96.614204      -9.185958       -23.33643       -23.033337      -101.282059     -101.583092     -101.279999

        // chr pos ref Nreads maxMapQ genotype BtrLod BtnbLod dbSNP AA      AC      AG      AT      CC      CG      CT      GG      GT      TT
        return String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0",
	                                        location.getContig(),
                                            location.getStart(),
											ref,
                                            depth,
                                            -1,
	                                        genotype(),
	                                        lodVsRef,
	                                        lodVsNextBest);
    }

    public String asTabularString() 
	{
        return String.format("%s %s %c %c %s %f %f %f %f %d %s",
	                                        location,
											sample_name,
	                                        ref,
	                                        alt,
											genotype(),
	                                        qhat,
	                                        qstar,
                                            lodVsRef,
                                            lodVsNextBest,
	                                        depth, 
											bases);
    }

    public String toString() { return asTabularString(); }

    public String asString() {
        // Print out the called bases
        // Notes: switched from qhat to qstar because qhat doesn't work at n=1 (1 observed base) where having a single non-ref
        //        base has you calculate qstar = 0.0 and qhat = 0.49 and that leads to a genotype predicition of AG according
        //        to qhat, but AA according to qstar.  This needs to be further investigated to see whether we really want
        //        to use qstar, but make N (number of chormosomes) switch to n (number of reads at locus) for n=1
        long numNonrefBases = Math.round(qstar * N);
        long numRefBases = N - numNonrefBases;
        if (ref < alt) { // order bases alphabetically
            return AlleleFrequencyWalker.repeat(ref, numRefBases) + AlleleFrequencyWalker.repeat(alt, numNonrefBases);
        }else{
            return AlleleFrequencyWalker.repeat(alt, numNonrefBases) + AlleleFrequencyWalker.repeat(ref, numRefBases);
        }
    }


    public double posterior()
    {
        return this.posteriors[(int)this.qstar * this.N];
    }
}

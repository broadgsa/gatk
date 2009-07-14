package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.playground.gatk.walkers.AlleleFrequencyWalker;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Arrays;

public class AlleleFrequencyEstimate {

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
    //public double[][] quals;
    public double[] posteriors;
	public String sample_name;
	public int n_ref;
	public int n_het;
	public int n_hom;

	public GenotypeLikelihoods genotypeLikelihoods = null;

    GenomeLoc l;


    public AlleleFrequencyEstimate(GenomeLoc location, char ref, char alt, int N, double qhat, double qstar, double lodVsRef, double lodVsNextBest, double pBest, double pRef, int depth, String bases, double[][] quals, double[] posteriors, String sample_name)
    {
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
        //this.quals = quals;
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
        //public double[] posteriors;
        return String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
	                                        location.getContig(),
                                            location.getStart(),
											ref,
                                            depth,
                                            -1,
	                                        genotype(),
	                                        lodVsRef,
	                                        lodVsNextBest,
                                            posteriors[0],
                                            posteriors[1],
                                            posteriors[2],
                                            posteriors[3],
                                            posteriors[4],
                                            posteriors[5],
                                            posteriors[6],
                                            posteriors[7],
                                            posteriors[8],
                                            posteriors[9]);
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

	public String asPoolTabularString()
	{
		return String.format("%s %c %c %f %f %f %s %f %d %d %d %d", 
								location,
								ref, 
								alt, 
								qstar, 
								pBest, 
								pRef, 
								"NA", 
								lodVsRef, 
								N, 
								n_ref, 
								n_het, 
								n_hom);
	}


    public double posterior()
    {
        return this.posteriors[(int)this.qstar * this.N];
    }

    public String callType() {
        // Returns a string indicating whether the call is homozygous reference, heterozygous SNP, or homozygous SNP

        String[] callTypeString = {"HomozygousSNP", "HeterozygousSNP", "HomozygousReference"};
        String genotype = genotype();
        int ref_matches = (genotype.charAt(0) == ref ? 1 : 0) + (genotype.charAt(1) == ref ? 1 : 0);
        return callTypeString[ref_matches];
    }

}

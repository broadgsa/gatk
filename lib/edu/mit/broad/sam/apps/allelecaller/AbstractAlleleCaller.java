package edu.mit.broad.sam.apps.allelecaller;

import edu.mit.broad.sam.SAMLocusIterator;
import edu.mit.broad.arachne.FastbReader;

import java.io.IOException;
import java.io.BufferedWriter;
import java.io.File;
import java.util.SortedSet;
import java.util.List;

/**
 * Base class for AlleleCallers.  Handles efficient access to the reference, output of data to a
 * standard file format, and application of priors
 */
public abstract class AbstractAlleleCaller {
    // writer for output
    private final BufferedWriter writer;

    // for providing access to reference data
    // TODO: replace with standard mechanism when defined/implemented
    private final FastbReader fastbReader;
    private String cachedChromName;
    private String cachedChrom;

    public AbstractAlleleCaller(final File fastbReference, final BufferedWriter writer) throws IOException {
        this.writer = writer;
        this.fastbReader = new FastbReader(fastbReference);
    }


    /**
     * emit allele calls to the writer specified in the constructor
     * 
     * @param li Locus to call
     */
    public void callAlleles(final SAMLocusIterator.LocusInfo li) throws IOException {

        // TODO: replace with standard mechanism when defined/implemented (making use of SAM Header)
        // make sure we have access to reference chrom information
        if (!li.getChrom().equals(cachedChromName)) {
            final int contig = translateChromToContig(li.getChrom());
            cachedChrom = null; // CRITICAL -- to allow for GC
            cachedChrom = fastbReader.readSequence(contig);
            cachedChromName = li.getChrom();
        }

        final char ref = cachedChrom.charAt(li.getPosition() - 1);


        // delegate to the specific implementation
        final SortedSet<GenotypeTheory> likelihoods = call(ref, li.getBasesAsString(), li.getQualities());


        final GenotypeTheory bestTheory = likelihoods.first();
        GenotypeTheory nextBestTheory = null;
        GenotypeTheory refTheory = null;
        final String refString = new String(new char[]{ref,ref});
        final DiploidGenotype refGenotype = DiploidGenotype.valueOf(refString);


        final StringBuilder theoryString = new StringBuilder();
        int k=0;
        for(final GenotypeTheory t : likelihoods) {
            if (k == 1) { nextBestTheory = t; }
            if (t.getGenotype() == refGenotype) { refTheory = t; }

            theoryString.append(t.getGenotype())
                        .append(":")
                        .append(String.format("%.2f",t.getLikelihood()))
                        .append(" ");
            k++;
        }

        final double btnb = bestTheory.getLikelihood() - nextBestTheory.getLikelihood();
        final double btr = bestTheory.getLikelihood() - refTheory.getLikelihood();

        final DiploidGenotype gt = likelihoods.first().getGenotype();

        final String type;
        if (!gt.isHet() && gt.getAllele1() == ref) {
            type = "homozygous";
        } else if (!gt.isHet() && gt.getAllele1() != ref) {
            type = "homozygous-SNP";
        } else {
            type = "heterozygous-SNP";
        }

        final String bases = li.getBasesAsString();
        int a = 0,c = 0,g = 0,t = 0;
        for(int i=0; i<bases.length(); i++) {
            if (bases.charAt(i) == 'A') { a++; }
            else if (bases.charAt(i) == 'C') { c++; }
            else if (bases.charAt(i) == 'G') { g++; }
            else if (bases.charAt(i) == 'T') { t++; }
            else { throw new RuntimeException("Unknown Base " + bases.charAt(i)); }
        }

        writer.write(
                translateChromToContig(li.getChrom()) + ":" +
                (li.getPosition()-1) + " " +   // arachne output is 0-based
                ref + " " +
                gt + " " +
                String.format("%f %f", btnb,btr) + " " +
                type + " " +
                "A:" + a + " " +
                "C:" + c + " " +
                "G:" + g + " " +
                "T:" + t + " " +
                bases.length() + " " +
                "0 1 1 " + // used prior, is alignable, bait present
                theoryString
        );


        writer.write("\n");
    }

    abstract protected SortedSet<GenotypeTheory> call(char ref, String bases, List<Byte> quals);
    

    /**
     * Apply a general population-based prior to the likelihood:
     * <ul>
     * <li>ref is .999</li>
     * <li>het is 10^-3</li>
     * <li>homozygous, non-reference is 10^-5</li>
     *
     * @param ref reference allele
     * @param allele1 first allele of the genotype
     * @param allele2 second allele of the genotype
     * @return prior, given the reference and genotype alleles
     */
    protected double getPrior(final char ref, final DiploidGenotype gt) {
        final double prior;
        if (gt.isHom() && gt.getAllele1() == ref) {
            prior = 0.999; // reference
        } else {
            if (gt.getAllele1() != ref && gt.getAllele2() != ref) {
                prior = 0.00001; // neither base is reference
            } else {
                prior = 0.001; // het, one base is reference
            }
        }
        return prior;
    }

    // --------------------------------------------------------------------------------------------
    // Helper methods below this point...
    // --------------------------------------------------------------------------------------------


    private final String[] chroms = new String[]{"chrM","chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY","chr1_random","chr2_random","chr3_random","chr4_random","chr5_random","chr6_random","chr7_random","chr8_random","chr9_random","chr10_random","chr11_random","chr13_random","chr15_random","chr16_random","chr17_random","chr18_random","chr19_random","chr21_random","chr22_random","chrX_random"};
    private int translateChromToContig(final String chrom) {
        for(int i=0; i<chroms.length; i++) {
            if (chrom.equals(chroms[i])) { return i; }
        }
        return -1;
    }

    public boolean isHet(final String alleles) {
        return (alleles.charAt(0) != (alleles.charAt(1)));
    }


}

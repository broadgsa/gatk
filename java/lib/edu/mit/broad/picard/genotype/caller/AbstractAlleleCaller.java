package edu.mit.broad.picard.genotype.caller;

import edu.mit.broad.picard.sam.SamLocusIterator;
import edu.mit.broad.sam.SAMFileHeader;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequence;
import edu.mit.broad.picard.PicardException;

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
    private final ReferenceSequenceFile referenceSequenceFile;
    private final SAMFileHeader samHeader;
    private ReferenceSequence referenceSequence;

    public AbstractAlleleCaller(final File reference, final SAMFileHeader samHeader, final BufferedWriter writer) {
        this.writer = writer;
        this.referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(reference);
        this.samHeader = samHeader;
    }


    /**
     * emit allele calls to the writer specified in the constructor
     * 
     * @param li Locus to call
     */
    public void callAlleles(final SamLocusIterator.LocusInfo li) throws IOException {


        cacheReferenceSequence(li.getSequenceIndex());

        final char ref = Character.toUpperCase((char)(referenceSequence.getBases()[li.getPosition() - 1] & 0xff));


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
                li.getSequenceIndex() + ":" +
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

    /**
     * Ensure that the referenceSequence member points to the sequenceIndex-th sequence.  Note that
     * this is not random access.  It is required that current sequenceIndex is >= the arg in the previous
     * call to this method.
     */
    private void cacheReferenceSequence(int sequenceIndex) {
        if (referenceSequence != null && referenceSequence.getContigIndex() == sequenceIndex) {
            return;
        }
        referenceSequence = null;
        for(referenceSequence = referenceSequenceFile.nextSequence();
                referenceSequence != null;
                referenceSequence = referenceSequenceFile.nextSequence()) {
            // Sanity check the sequence names against the sequence dictionary while scanning through.
            if (!referenceSequence.getName().equals(samHeader.getSequence(referenceSequence.getContigIndex()).getSequenceName())) {
                throw new PicardException("Sequence name mismatch at sequence index " + referenceSequence.getContigIndex() +
                ": " + referenceSequence.getName() + " != " +
                        samHeader.getSequence(referenceSequence.getContigIndex()).getSequenceName());
            }
            if (referenceSequence.getContigIndex() == sequenceIndex) {
                break;
            }
            if (referenceSequence.getContigIndex() > sequenceIndex) {
                throw new PicardException("Never found reference sequence with index " + sequenceIndex);
            }
        }
        if (referenceSequence == null) {
            throw new PicardException("Reference sequence with index " + sequenceIndex + " was not found");
        }
    }

    /**
     * Override this to implement a concrete genotype caller
     * @param ref the reference base
     * @param bases each element in the String is the base at current locus for a given read
     * @param quals same length as bases. the ith element corresponds to the ith element of bases.
     * @return
     */
    abstract protected SortedSet<GenotypeTheory> call(char ref, String bases, List<Byte> quals);
    

    /**
     * Apply a general population-based prior to the likelihood:
     * <ul>
     * <li>ref is .999</li>
     * <li>het is 10^-3</li>
     * <li>homozygous, non-reference is 10^-5</li>
     *
     * @param ref reference allele
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


    public boolean isHet(final String alleles) {
        return (alleles.charAt(0) != (alleles.charAt(1)));
    }


}

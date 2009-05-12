package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.HapMapAlleleFrequenciesROD;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import net.sf.samtools.SAMRecord;

import java.util.List;
import java.util.HashMap;
import java.util.Map;

public class SomaticMutationWalker extends LocusWalker<Integer, Integer> {
    protected static class QualitySums {
        private int a = 0;
        private int c = 0;
        private int g = 0;
        private int t = 0;

        public int get(final char base) {
            if (base == 'a' || base == 'A') { return a; }
            if (base == 'c' || base == 'C') { return c; }
            if (base == 'g' || base == 'G') { return g; }
            if (base == 't' || base == 'T') { return t; }
            throw new RuntimeException("Unknown base: " + base);
        }

        public void incrementSum(final char base, final byte qual) {
            if (base == 'a' || base == 'A') { a += qual; }
            else if (base == 'c' || base == 'C') { c += qual; }
            else if (base == 'g' || base == 'G') { g += qual; }
            else if (base == 't' || base == 'T') { t += qual; }
            else throw new RuntimeException("Unknown base: " + base);
        }

        public int getOtherQualities(final char base) {
            int total = a + c + g + t;
            if (base == 'a' || base == 'A') { return total-a; }
            else if (base == 'c' || base == 'C') { return total-c; }
            else if (base == 'g' || base == 'G') { return total-g; }
            else if (base == 't' || base == 'T') { return total-t; }
            else throw new RuntimeException("Unknown base: " + base);
        }

        public void reset() {
            a = 0; c = 0; g = 0; t = 0;
        }
    }

    @Argument(fullName = "tumor_sample_name", shortName = "s1", required = true, doc="Name of the tumor sample")
    public String tumorSampleName;

    @Argument(fullName = "normal_sample_name", shortName = "s2", required = true, doc="Name of the normal sample")
    public String normalSampleName;

    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public static int MAX_INSERT_SIZE = 10000;
    public static int MIN_MUTANT_SUM_PRETEST = 60;
    public static int MIN_MUTANT_SUM = 100;

    public static double TUMOR_LOD = 6.3d;
//    public static double TUMOR_LOD = 1.0d;
    public static double NORMAL_LOD = 2.3d;


    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        Map<Character, Integer> midp = new HashMap<Character, Integer>();
        midp.put('A', Integer.MAX_VALUE);
        midp.put('C', Integer.MAX_VALUE);
        midp.put('G', Integer.MAX_VALUE);
        midp.put('T', Integer.MAX_VALUE);

        // First, exclude DBSNP Sites
        for ( ReferenceOrderedDatum datum : tracker.getAllRods() )
        {
            if ( datum != null )
            {
                if ( datum instanceof rodDbSNP) {
                    // we're at a dbsnp site... move along... nothing to see here... 
                    return -1;
                }
            }
        }

        char upRef = Character.toUpperCase(ref);

        // TODO: should we be able to call mutations at bases where ref is N?
        if (upRef == 'N') { return -1; }
        
        List<SAMRecord> reads = context.getReads();

        GenotypeLikelihoods tumorGL = new GenotypeLikelihoods();
        GenotypeLikelihoods normalGL = new GenotypeLikelihoods();

        QualitySums tumorQualitySums = new QualitySums();
        QualitySums normalQualitySums = new QualitySums();

        StringBuilder tumorBases = new StringBuilder();
        StringBuilder normalBases = new StringBuilder();
        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() <= 0

//               || (read.getReadPairedFlag() && (!read.getProperPairFlag() || read.getInferredInsertSize() >= MAX_INSERT_SIZE))
                    ) {
                continue;
            }

            String rg = (String) read.getAttribute("RG");
            String sample = read.getHeader().getReadGroup(rg).getSample();

            int offset = context.getOffsets().get(i);
            char base = read.getReadString().charAt(offset);
            byte qual = read.getBaseQualities()[offset];


            if (base == 'N' || base == 'n') { continue; }
            
            if (normalSampleName.equals(sample)) {
                normalGL.add(ref, base, qual);
                normalQualitySums.incrementSum(base, qual);
                normalBases.append(base);
            } else if (tumorSampleName.equals(sample)) {
                tumorGL.add(ref, base, qual);
                tumorQualitySums.incrementSum(base, qual);
                tumorBases.append(base);

                int midDist = Math.abs((int)(read.getReadLength() / 2) - offset);
                if (midDist < midp.get(base)) { midp.put(base, midDist); }

            } else {
                throw new RuntimeException("Unknown Sample Name: " + sample);
            }
        }

        // pretest: if the sum of the quality scores for all non-ref alleles < 60, just quit looking now
        if (tumorQualitySums.getOtherQualities(ref) < MIN_MUTANT_SUM_PRETEST) {
            return -1;
        }


        // Test each of the poosible alternate alleles

        for (char altAllele : new char[]{'A','C','G','T'}) {
            if (altAllele == upRef) { continue; }


            double[] tumorRefAlt = extractRefAlt(tumorGL, ref, altAllele);
            double tumorLod = Math.log10(tumorRefAlt[1] + tumorRefAlt[2]) - tumorRefAlt[0];

            // (i) either an adjusted quality score sum in the tumor for the mutant base must be
            //     at least 100 or the LOD score for mutant:ref + mutant:mutant vs ref:ref must
            //     be at least 6.3;
            if (tumorQualitySums.get(altAllele) < MIN_MUTANT_SUM && tumorLod < TUMOR_LOD) {
                continue;
            }

            // make sure we've seen at least 1 obs of the alternate allele within 20bp of the read-middle
            if (midp.get(altAllele) > 20) {
                out.println("Rejecting due to midpoint check!");
                return 0;
            }

            double[] refAlt = extractRefAlt(normalGL, ref, altAllele);
            double normalLod = refAlt[0] - Math.log10(refAlt[1] + refAlt[2]);

            // (ii) the quality score sum for the mutant base in the normal must be < 50 and the
            //      LOD score for ref:ref vs mutant:ref + mutant:mutant must be at least 2.3.
            if ( normalQualitySums.get(altAllele) > 50 || normalLod < NORMAL_LOD ) {
                continue;
            }

            // if we're still here... we've got a somatic mutation!  Output the results
            // and stop looking for mutants!
            out.println(context.getLocation() + " " + upRef + " " + altAllele +
            " TScore:" + tumorLod +
            " TRefSum: " + tumorQualitySums.get(ref) +
            " TAltSum: " + tumorQualitySums.get(altAllele) +
            " NScore:" + normalLod +
            " NRefSum: " + normalQualitySums.get(ref) +
            " NAltSum: " + normalQualitySums.get(altAllele) + " " +
                    tumorBases.toString() + " " +
                    normalBases.toString()
                    );


            return 1;


        }

        return -1;
    }

    /**
     * Extract the LOD comparing ref:ref to ref:alt and alt:alt
     */
    private double[] extractRefAlt(GenotypeLikelihoods gl, char ref, char altAllele) {
        double refRef = 0;
        double altRef = 0;
        double altAlt = 0;
        for(int j=0; j<10; j++) {
            String gt = gl.genotypes[j];
            double likelihood = gl.likelihoods[j];

            // the ref:mutant theory
            if ( (gt.charAt(0) == ref && gt.charAt(1) == altAllele) ||
                 (gt.charAt(0) == altAllele && gt.charAt(1) == ref) ) {
                altRef += Math.pow(10, likelihood);
            }

            if ( gt.charAt(0) == altAllele && gt.charAt(1) == altAllele) {
                altAlt += Math.pow(10, likelihood);
            }

            if ( gt.charAt(0) == ref && gt.charAt(1) == ref) {
                refRef = likelihood;
            }

        }
        return new double[]{refRef, altRef, altAlt};
    }

    // Given result of map function
    public Integer reduceInit() {
        return 0;
    }
    public Integer reduce(Integer value, Integer sum) {
        return value + sum;
    }

    @Override
    public void onTraversalDone(Integer result) {
//        out.println(String.format("FINAL - %d %d %d %d", totalSites, tumorCovered, normalCovered, somaticCovered));
    }



}
package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.StringUtil;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.playground.indels.SWPairwiseAlignment;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.io.IOException;
import java.util.*;

@By(DataSource.REFERENCE)
public class  SomaticMutationWalker extends LocusWalker<Integer, Integer> {
    protected enum MutantFailureReason {
        StrandImbalance,
        Misalignment
    }
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

    @Argument(fullName = "tumor_lod", required = false, doc = "LOD threshold for calling tumor variant")
    public float TUMOR_LOD_THRESHOLD = 6.3f;

    @Argument(fullName = "normal_lod", required = false, doc = "LOD threshold for calling normal non-variant")
    public float NORMAL_LOD_THRESHOLD = 2.3f;

    @Argument(fullName = "min_mutant_sum", required = false, doc = "threshold for sum of mutant allele quality scores")
    public int MIN_MUTANT_SUM = 100;

    @Argument(fullName = "mode", required = false, doc="Mode of operation (detect, full)")
    public String mode = "full";

    public float SKEW_LOD_THRESHOLD = 1.5f;

//    @Argument(fullName = "output_failures", required = false, doc="produce output for failed sites")
    public boolean OUTPUT_FAILURES = true;


    public void initialize() {
    }

    public String walkerType() { return "ByLocus"; }

    // Do we actually want to operate on the context?
    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return true;    // We are keeping all the reads
    }

    public static int MAX_INSERT_SIZE = 10000;
    public static int MIN_MUTANT_SUM_PRETEST = 60;

    public static int MIN_QSCORE = 13;


    private static class LocusReadPile {
        private char refBase;
        public List<SAMRecord> reads = new ArrayList<SAMRecord>();
        public List<Integer> offsets = new ArrayList<Integer>();
        public GenotypeLikelihoods likelihoods = new GenotypeLikelihoods();
        public QualitySums qualitySums = new QualitySums();

        public LocusReadPile(char refBase) {
            this.refBase = refBase;
        }

        public void add(SAMRecord read, int offset) {
            char base = read.getReadString().charAt(offset);
            byte qual = read.getBaseQualities()[offset];

            if (base == 'N' || base == 'n') { return; }

            reads.add(read);
            offsets.add(offset);


            likelihoods.add(refBase, base, qual);
            if (qual > MIN_QSCORE) qualitySums.incrementSum(base, qual);

        }

        public String getLocusBases() {
            return getLocusBases(0);
        }

        public String getLocusBases(int locusOffset) {
            StringBuilder sb = new StringBuilder();
            for(int i=0; i<reads.size(); i++) {
                SAMRecord read = reads.get(i);
                int readOffset = offsets.get(i);

                int offset = readOffset + locusOffset;
                if (offset >= 0 && offset < read.getReadString().length()) {
                    char base = read.getReadString().charAt(offset);
                    sb.append(base);
                }
            }
            return sb.toString();
        }

        public double getAltVsRef(char altAllele) {
            double[] tumorRefAlt = extractRefAlt(this.likelihoods, this.refBase, altAllele);
            double tumorLod = Math.log10(tumorRefAlt[1] + tumorRefAlt[2]) - tumorRefAlt[0];
            return tumorLod;
        }

        public double getRefVsAlt(char altAllele) {
            double[] refAlt = extractRefAlt(this.likelihoods, this.refBase, altAllele);
            double normalLod = refAlt[0] - Math.log10(refAlt[1] + refAlt[2]);
            return normalLod;
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

        private GenotypeLikelihoods getLikelihood(int locusOffset) {
            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods();


            for(int i=0; i<reads.size(); i++) {
                SAMRecord read = reads.get(i);
                int readOffset = offsets.get(i);

                int offset = readOffset + locusOffset;
                if (offset >= 0 && offset < read.getReadString().length()) {

                    char base = read.getReadString().charAt(offset);
                    byte qual = read.getBaseQualities()[offset];

                    likelihoods.add(refBase, base, qual);
                }

            }
            return likelihoods;
        }

        public double[] getNormalizedProbs(int locusOffset, int qThreshold) {
            GenotypeLikelihoods likelihoods = new GenotypeLikelihoods();


            // FIXME: gaddy suggested this "correction" which evidently has a name...
            //likelihoods.add(refBase, refBase, (byte) 20);
            
            for(int i=0; i<reads.size(); i++) {
                SAMRecord read = reads.get(i);
                int readOffset = offsets.get(i);

                int offset = readOffset + locusOffset;
                if (offset >= 0 && offset < read.getReadString().length()) {

                    char base = read.getReadString().charAt(offset);
                    byte qual = read.getBaseQualities()[offset];

                    if (qual >= qThreshold) {
                        likelihoods.add(refBase, base, qual);
//                        if (locusOffset == -43) System.out.println("\tUSING " + base + " " + qual);
                    } else {
//                        if (locusOffset == -43) System.out.println("\tDropping " + base + " " + qual);
                    }
                }

            }


            double[] logLikelihood = likelihoods.likelihoods;
            double[] nonLogLikelihood = new double[10];
            double sum = 0;
            for(int i=0; i<10; i++) {
                nonLogLikelihood[i] = Math.pow(10, logLikelihood[i]);
                sum += nonLogLikelihood[i];
            }

            double[] normalizedProbs = new double[10];
            for(int i=0; i<10; i++) {
                normalizedProbs[i] = nonLogLikelihood[i] / sum;
            }


            //quick sanity check
//            sum=0;
//            for(int i=0; i<10; i++) {
//                sum += normalizedProbs[i];
//            }
//            System.out.println("normalized probs = " + sum);


            return normalizedProbs;
        }
        
    }

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

        LocusReadPile tumorReadPile = new LocusReadPile(upRef);
        LocusReadPile normalReadPile = new LocusReadPile(upRef);

        for ( int i = 0; i < reads.size(); i++ )
        {
            SAMRecord read = reads.get(i);

            if (read.getNotPrimaryAlignmentFlag() ||
                read.getDuplicateReadFlag() ||
                read.getReadUnmappedFlag() ||
                read.getMappingQuality() <= 0

               || (read.getReadPairedFlag() && (!read.getProperPairFlag() || read.getInferredInsertSize() >= MAX_INSERT_SIZE))
                    ) {
                continue;
            }

            String rg = (String) read.getAttribute("RG");
            String sample = read.getHeader().getReadGroup(rg).getSample();

            int offset = context.getOffsets().get(i);

            char base = read.getReadString().charAt(offset);
            if (base == 'N' || base == 'n') { continue; }


            // TODO: build a pile of reads and offsets, then pass that into a
            // constructor for the normalGL class
            // that way, we can build a different pile of reads later on and extract the genotype
            if (normalSampleName.equals(sample)) {
                normalReadPile.add(read, offset);

            } else if (tumorSampleName.equals(sample)) {
                tumorReadPile.add(read, offset);


                int midDist = Math.abs((int)(read.getReadLength() / 2) - offset);
                if (midDist < midp.get(base)) { midp.put(base, midDist); }

            } else {
                throw new RuntimeException("Unknown Sample Name: " + sample);
            }
        }

        // pretest: if the sum of the quality scores for all non-ref alleles < 60, just quit looking now
        if (tumorReadPile.qualitySums.getOtherQualities(upRef) < MIN_MUTANT_SUM_PRETEST) {
            return -1;
        }


        // Test each of the poosible alternate alleles

        for (char altAllele : new char[]{'A','C','G','T'}) {
            if (altAllele == upRef) { continue; }


            // (i) either an adjusted quality score sum in the tumor for the mutant base must be
            //     at least 100 or the LOD score for mutant:ref + mutant:mutant vs ref:ref must
            //     be at least 6.3;
            int mutantSum = tumorReadPile.qualitySums.get(altAllele);
            int refSum = tumorReadPile.qualitySums.get(upRef);

            if (tumorReadPile.getAltVsRef(altAllele) >= TUMOR_LOD_THRESHOLD
//                    ||
//                (mutantSum >= MIN_MUTANT_SUM && (float)mutantSum / (float) refSum >= 0.05f)
                ) {
                // yep -- just fall through... easier to think about this way!
            } else {
                continue;
            }

            // (ii) the quality score sum for the mutant base in the normal must be < 50 and the
            //      LOD score for ref:ref vs mutant:ref + mutant:mutant must be at least 2.3.
            double normalLod = normalReadPile.getRefVsAlt(altAllele);
            if ( normalReadPile.qualitySums.get(altAllele) > 50 || normalLod < NORMAL_LOD_THRESHOLD) {
                continue;
            }

            // make sure we've seen at least 1 obs of the alternate allele within 20bp of the read-middle
            boolean failedMidpointCheck = midp.get(altAllele) > 20;
//            if (failedMidpointCheck) {
//                out.println("Rejecting due to midpoint check!");
//                return 0;
//            }


            // do a MSA to figure out if the alternate reads comes from a cluster of reads with more
            // than one alternate allele (hints at being an alignment artifact)
            ReferenceSequence refSeq;
            // TODO: don't hardcode.  Make this the max read length in the pile 
            long refStart = context.getLocation().getStart() - 150;
            //tumorReadPile.offsets.get(0);
            long refStop = context.getLocation().getStart() + 150;
            try {
                IndexedFastaSequenceFile seqFile = new IndexedFastaSequenceFile(getToolkit().getArguments().referenceFile);
                refSeq = seqFile.getSubsequenceAt(context.getContig(),refStart, refStop);

            } catch (IOException ioe) {
                throw new RuntimeException(ioe);
            }

//            System.out.println("TESTING " + context.getContig() + ":" + context.getPosition()); 
            LocusReadPile t2 = filterHighMismatchScoreReads(tumorReadPile, StringUtil.bytesToString(refSeq.getBases()), refStart);

            // TEST the LOD score again!
            // TODO: extract this since we'll do it multiple times...

            // (i) either an adjusted quality score sum in the tumor for the mutant base must be
            //     at least 100 or the LOD score for mutant:ref + mutant:mutant vs ref:ref must
            //     be at least 6.3;
            double tumorLod = t2.getAltVsRef(altAllele);
            if (mode.equals("full") && t2.qualitySums.get(altAllele) < MIN_MUTANT_SUM && tumorLod < TUMOR_LOD_THRESHOLD) {
                if (OUTPUT_FAILURES) {

                    String msg = "FAILED  due to MAX MM QSCORE TEST." +
                            " LOD was " + tumorReadPile.getAltVsRef(altAllele) +
                            " LOD is now " + t2.getAltVsRef(altAllele) +
                            " QSUM was " + tumorReadPile.qualitySums.get(altAllele) +
                            " QSUM is now " + t2.qualitySums.get(altAllele);

                    out.println(
                            context.getContig() + "\t" +
                                    context.getPosition() + "\t" +
                                    context.getPosition() + "\t"
                            + msg.replace(' ','_')
                            );
                }
                continue;
            }

//            // TODO: using the original pile here since often these artifacts will be supported
//            // by those reads that get thrown out!  Maybe that means we don't need the noise filter...
//            boolean shouldDisalign =
//                    disaligner(context.getPosition(), tumorReadPile, StringUtil.bytesToString(refSeq.getBases()), refStart);
            Pair<MutantFailureReason,Double> failureReason =
                    readPileSkew(t2, upRef, altAllele, StringUtil.bytesToString(refSeq.getBases()), refStart);
                    

            if (mode.equals("full") && failureReason.first != null) {
                if (OUTPUT_FAILURES) {
                    String msg = "FAILED due to " + failureReason.first.name() +
                            " mutAllele " + altAllele + "_" +
                            " maxSkewLod " + failureReason.second
                            ;

                    out.println(
                            context.getContig() + "\t" +
                                    context.getPosition() + "\t" +
                                    context.getPosition() + "\t" +
                                    msg.replace(' ','_')
                            );

                }
                continue;
            }







            // if we're still here... we've got a somatic mutation!  Output the results
            // and stop looking for mutants!
            String msg =
                        (failedMidpointCheck?"__FAILED-MPCHECK":"") +
                        "TScore:" + tumorLod +
                        "__TRefSum: " + tumorReadPile.qualitySums.get(upRef) +
                        "__TAltSum: " + tumorReadPile.qualitySums.get(altAllele) +
                        "__NScore:" + normalLod +
                        "__NRefSum: " + normalReadPile.qualitySums.get(upRef) +
                        "__NAltSum: " + normalReadPile.qualitySums.get(altAllele) +
                        "__maxSkewLod_" + failureReason.second + "_" +
                        "__MIDP: " + midp.get(altAllele);

            out.println(
                    context.getContig() + "\t" +
                            context.getPosition() + "\t" +
                            context.getPosition() + "\t"
                    + msg.replace(' ','_')
                    );




            return 1;


        }

        return -1;
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



    private Pair<MutantFailureReason, Double> readPileSkew(LocusReadPile pile, char refAllele, char mutantAllele, String reference, long leftmostIndex) {
        // first split into two piles, those supporting the mutant and those not
        LocusReadPile mutantPile = new LocusReadPile(mutantAllele);
        LocusReadPile refPile =  new LocusReadPile(mutantAllele);


        for (int i=0; i<pile.reads.size(); i++) {
            SAMRecord read = pile.reads.get(i);
            int offset = pile.offsets.get(i);

            if (read.getReadString().charAt(offset) == mutantAllele) {
                mutantPile.add(read, offset);
            } else if (read.getReadString().charAt(offset) == refAllele) {
                refPile.add(read, offset);
            } else {
                // just drop the read...
            }
        }

        // now we can ask questons abot the two piles.

        // e.g. is the mutant allele seen in both strands?
        boolean seenOnPositive = false;
        boolean seenOnNegative = false;
        for(SAMRecord read : mutantPile.reads) {
            if (read.getReadNegativeStrandFlag()) { seenOnNegative = true; } else { seenOnPositive = true; }
        }
        if (!seenOnPositive || !seenOnNegative) {
//            return MutantFailureReason.StrandImbalance;
        }


        //chr17:4979257
        // are the two alleles distributed differently
        //fixme: this seems to degrade as you lose reads?
        //fixme: what should the threshold be here?
        SortedMap<Integer, Double> skewLodOffsets = new TreeMap<Integer, Double>();
        double maxSkewLod = 0;

        //fixme: calculate this range properly
        int MAX_OFFSET_DISTANCE = 60;
        for(int offset=-1 * MAX_OFFSET_DISTANCE; offset<MAX_OFFSET_DISTANCE; offset++) {
            // allow for doubletons
            if (offset >= -1 && offset <= 1 ) { continue; }

            int SKEW_QSCORE_THRESHOLD = 10;
            double[] mutantNormProbs = mutantPile.getNormalizedProbs(offset, SKEW_QSCORE_THRESHOLD);
            double[] otherNormProbs = refPile.getNormalizedProbs(offset, SKEW_QSCORE_THRESHOLD);

            double J = 0;
            for(int i=0; i<10; i++) {
                J += mutantNormProbs[i] * otherNormProbs[i];
            }
            double skewLod = Math.log10( (1-J) / J);


            // fixme: is it legit to require that we see it in at least 2 reads each?
            int mutantReadCounts = mutantPile.getLocusBases(offset).length();
            int otherReadCounts = refPile.getLocusBases(offset).length();
            if (mutantReadCounts >= 2 && otherReadCounts >= 2) {

                if (skewLod > maxSkewLod) { maxSkewLod = skewLod; }
                if (skewLod > SKEW_LOD_THRESHOLD) {

//                    System.out.println( "Offset: " + offset +
//                            " mutant_reads: " + mutantReadCounts +
//                            " mutant bases: " + mutantPile.getLocusBases(offset) +
//                            " other_reads: " + otherReadCounts +
//                            " other bases: " + refPile.getLocusBases(offset) +
//                            " skewLod:" + skewLod );

                    skewLodOffsets.put(offset, skewLod);
                }
            }
        }

        if (skewLodOffsets.size() > 0) {
            return new Pair<MutantFailureReason, Double>(MutantFailureReason.Misalignment, maxSkewLod);
        } else {
            return new Pair<MutantFailureReason, Double>(null, maxSkewLod);
//            System.out.println("MAX SKEW: " + maxSkewLod);
        }

    }

    int MAX_READ_MISMATCH_QUALITY_SCORE_SUM = 100;
    private LocusReadPile filterHighMismatchScoreReads(LocusReadPile pile, String reference, long leftmostIndex) {
        LocusReadPile newPile = new LocusReadPile(pile.refBase);

        for ( int i=0; i< pile.reads.size(); i++) {
            SAMRecord read = pile.reads.get(i);
            int offset = pile.offsets.get(i);

            AlignedRead aRead = new AlignedRead(read);
            int mismatchScore = mismatchQualitySum(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex);

            if (mismatchScore <= MAX_READ_MISMATCH_QUALITY_SCORE_SUM) {
                newPile.add(read, offset);
            } else {
                char base = read.getReadString().charAt(offset);

                // TODO: use a logger to output this information for debugging/Picard folks
                //out.println("MMQSUM: Filtering out " + read.getReadName() + " with locus base " + base + " and ref of " + pile.refBase);
            }
        }
        return newPile;
    }


    // also pass in the variant position
    private boolean disaligner(long mutationLocus, LocusReadPile pile, String reference, long leftmostIndex) {
        int VARIANT_SITE_MIN_QSCORE = 40; // ie 2x20
        // build a 2d array of read # vs position
//        Map<Long, Map<Character, List<Mismatch>>> matrix = new TreeMap<Long, Map<Character, List<Mismatch>>>();

        int mutationLocusOffset = (int)(mutationLocus - leftmostIndex);
        int[][] mmm = new int[pile.reads.size()][reference.length()];

        // for every read which has an offset
        

        for ( int i=0; i< pile.reads.size(); i++) {
            SAMRecord read = pile.reads.get(i);
            int offset = pile.offsets.get(i);

            AlignedRead aRead = new AlignedRead(read);
            List<Mismatch> mismatches = getMismatches(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex);

            // now fill in with -1 where we had reference
            // TODO: deal with other CIGAR strings?
            for(long j=read.getAlignmentStart(); j <= read.getAlignmentEnd(); j++) {
                mmm[i][(int)( j - leftmostIndex )] = -1;
            }


            for(Mismatch mm : mismatches) {
                if (mm.mismatchBase == 'N') { continue; }

                mmm[i][(int)mm.position] = mm.qualityScore;
            }
        }

        // traverse in column-order
        int[] altSums = new int[mmm[0].length];
        for(int j=0; j<altSums.length; j++) {

            for(int i=0; i<mmm.length; i++) {
                if (mmm[i][j] > 0) altSums[j] += mmm[i][j];
            }
        }

        // now check that we have at least 3 alignment artifacts, and that
        // they are at least 5bp from eachother
        // FIXME: what should these numbers be?
        int varSites = 0;
        int lastVarSite = -1;
        for(int j=0; j<altSums.length; j++) {
            if (altSums[j] >= VARIANT_SITE_MIN_QSCORE) {
                if (lastVarSite > -1) {
                    if (j - lastVarSite > 5) {
                        varSites++;
                        lastVarSite = j;
                    }
                }
            }
        }
        if (varSites < 3) { return false; }
        
        // now go through each read which has the alternate allele at the
        // current site (the potential somatic mutation) and see what fraction of
        // the "covered" variant bases were called as variant
        Map<Integer, Float> fractionVariantSites = new HashMap<Integer, Float>();
        for(int i=0; i<mmm.length; i++) {

            // does this read have the variant allele?
            if (mmm[i][mutationLocusOffset] > 0) {
                float calledSites = 0;
                float variantSites = 0;
                for(int j=0; j<altSums.length; j++) {
                    if (altSums[j] > VARIANT_SITE_MIN_QSCORE && j != mutationLocusOffset) {
                        if (mmm[i][j] != 0) { calledSites++; }
                        if (mmm[i][j] > 0) { variantSites++; }
                    }
                }
                fractionVariantSites.put(i,  variantSites / calledSites);
            }
        }

        // FIXME: what's the right statistic to use here?
        if (getMedian(fractionVariantSites.values()) >= 0.75) {
            return true;
        }

        return false;

    }

    protected float getMedian(Collection<Float> in) {
        float[] data = new float[in.size()];
        int i=0;
        for(float f : in) { data[i++] = f; }

        Arrays.sort(data);

        int middle = data.length/2;  // subscript of middle element
        if (data.length%2 == 1) {
            return data[middle];
        } else {
           // Even number -- return average of middle two
           // Must cast the numbers to double before dividing.
           return (data[middle-1] + data[middle]) / 2.0f;
        }
    }


    private void clean(LocusReadPile pile, String reference, long leftmostIndex) {
//        LocusReadPile newPile = new LocusReadPile(pile.refBase);

        ArrayList<SAMRecord> refReads = new ArrayList<SAMRecord>();
        ArrayList<AlignedRead> altReads = new ArrayList<AlignedRead>();
        ArrayList<Boolean> altAlignmentsToTest = new ArrayList<Boolean>();
        int totalMismatchSum = 0;

        int MIN_MISMATCH_QSCORE = 10;
        // decide which reads potentially need to be cleaned
        for ( int i=0; i< pile.reads.size(); i++) {
            SAMRecord read = pile.reads.get(i);
            int offset = pile.offsets.get(i);

            AlignedRead aRead = new AlignedRead(read);
            int mismatchScore = mismatchQualitySum(aRead, reference, read.getAlignmentStart()-(int)leftmostIndex, MIN_MISMATCH_QSCORE);

            // if this doesn't match perfectly to the reference, let's try to clean it
            if ( mismatchScore > 0 ) {

                if ( mismatchScore > 100) {
                    out.println(read.getReadName() + " has a sum of mismatch quality scores of " + mismatchScore);

                } else {
                    altReads.add(aRead);
                    altAlignmentsToTest.add(true);
                    totalMismatchSum += mismatchScore;
                    aRead.setMismatchScoreToReference(mismatchScore);
                }
            }
            // otherwise, we can emit it as is
            else {
                refReads.add(read);
            }
        }

        // build a map of position ->
        Consensus bestConsensus = null;

        // for each alternative consensus to test, align it to the reference and create an alternative consensus
        for ( int index = 0; index < altAlignmentsToTest.size(); index++ ) {
            if ( altAlignmentsToTest.get(index) ) {

                // do a pairwise alignment against the reference
                AlignedRead aRead = altReads.get(index);
                SWPairwiseAlignment swConsensus = new SWPairwiseAlignment(reference, aRead.getReadString());
                int refIdx = swConsensus.getAlignmentStart2wrt1();
                if ( refIdx < 0 )
                    continue;

                // create the new consensus
                StringBuffer sb = new StringBuffer();
                sb.append(reference.substring(0, refIdx));
                Cigar c = swConsensus.getCigar();
                logger.debug("CIGAR = " + cigarToString(c));

                int indelCount = 0;
                int altIdx = 0;
                boolean ok_flag = true;
                for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
                    CigarElement ce = c.getCigarElement(i);
                    switch( ce.getOperator() ) {
                        case D:
                            indelCount++;
                            refIdx += ce.getLength();
                            break;
                        case M:
                            if ( reference.length() < refIdx+ce.getLength() )
                                ok_flag = false;
                            else
                                sb.append(reference.substring(refIdx, refIdx+ce.getLength()));
                            refIdx += ce.getLength();
                            altIdx += ce.getLength();
                            break;
                        case I:
                            sb.append(aRead.getReadString().substring(altIdx, altIdx+ce.getLength()));
                            altIdx += ce.getLength();
                            indelCount++;
                            break;
                    }
                }
                // make sure that there is at most only a single indel and it aligns appropriately!
                if ( !ok_flag || indelCount > 1 || reference.length() < refIdx )
                    continue;

                sb.append(reference.substring(refIdx));
                String altConsensus =  sb.toString();

                // for each imperfect match to the reference, score it against this alternative
                Consensus consensus = new Consensus(altConsensus, c, swConsensus.getAlignmentStart2wrt1());
                for ( int j = 0; j < altReads.size(); j++ ) {
                    AlignedRead toTest = altReads.get(j);
                    Pair<Integer, Integer> altAlignment = findBestOffset(altConsensus, toTest);

                    // the mismatch score is the min of its alignment vs. the reference and vs. the alternate
                    int myScore = altAlignment.getSecond();
                    if ( myScore >= toTest.getMismatchScoreToReference() )
                        myScore = toTest.getMismatchScoreToReference();
                    // keep track of reads that align better to the alternate consensus
                    else
                        consensus.readIndexes.add(new Pair<Integer, Integer>(j, altAlignment.getFirst()));

                    logger.debug(aRead.getReadString() +  " vs. " + toTest.getReadString() + " => " + myScore + " - " + altAlignment.getFirst());
                    consensus.mismatchSum += myScore;
                    if ( myScore == 0 )
                        // we already know that this is its consensus, so don't bother testing it later
                        altAlignmentsToTest.set(j, false);
                }
                logger.debug(aRead.getReadString() +  " " + consensus.mismatchSum);
                if ( bestConsensus == null || bestConsensus.mismatchSum > consensus.mismatchSum) {
                    bestConsensus = consensus;
                    logger.debug(aRead.getReadString() +  " " + consensus.mismatchSum);
                }
            }
        }
        System.out.println("Done with MSA");
//        // if the best alternate consensus has a smaller sum of quality score mismatches (more than the LOD threshold), then clean!
//        if ( bestConsensus != null && ((double)(totalMismatchSum - bestConsensus.mismatchSum))/10.0 >= LOD_THRESHOLD ) {
//            logger.debug("CLEAN: " + bestConsensus.str );
//            if ( indelOutput != null && bestConsensus.cigar.numCigarElements() > 1 ) {
//                StringBuffer str = new StringBuffer();
//                str.append(reads.get(0).getReferenceName());
//                int position = bestConsensus.positionOnReference + bestConsensus.cigar.getCigarElement(0).getLength();
//                str.append(":" + (leftmostIndex + position));
//                CigarElement ce = bestConsensus.cigar.getCigarElement(1);
//                str.append("\t" + ce.getLength() + ce.getOperator());
//                str.append("\t" + (((double)(totalMismatchSum - bestConsensus.mismatchSum))/10.0) + "\n");
//                try {
//                    indelOutput.write(str.toString());
//                    indelOutput.flush();
//                } catch (Exception e) {}
//            }
//
//            // We need to update the mapping quality score of the cleaned reads;
//            // however we don't have enough info to use the proper MAQ scoring system.
//            // For now, we'll use a heuristic:
//            // the mapping quality score is improved by the LOD difference in mismatching
//            // bases between the reference and alternate consensus
//            int improvement = (totalMismatchSum - bestConsensus.mismatchSum) / 10;
//
//            // clean the appropriate reads
//            for ( Pair<Integer, Integer> indexPair : bestConsensus.readIndexes ) {
//                AlignedRead aRead = altReads.get(indexPair.getFirst());
//                updateRead(bestConsensus.cigar, bestConsensus.positionOnReference, indexPair.getSecond(), aRead, (int)leftmostIndex);
//                aRead.getRead().setMappingQuality(Math.min(aRead.getRead().getMappingQuality() + improvement, 255));
//                aRead.getRead().setAttribute("NM", numMismatches(aRead.getRead(), reference, aRead.getRead().getAlignmentStart()-(int)leftmostIndex));
//            }
//        }
//
//        // write them out
//        for ( SAMRecord rec : refReads )
//            readsToWrite.add(new ComparableSAMRecord(rec));
//        for ( AlignedRead aRec : altReads )
//            readsToWrite.add(new ComparableSAMRecord(aRec.getRead()));
    }



    // ----------------------------------- PRIVATE IN IntervalCleanerWalker
    public static final int MAX_QUAL = 99;

    private static class Mismatch {
        AlignedRead aRead;
        int offset;
        long position;
        char mismatchBase;
        int qualityScore;

        private Mismatch(AlignedRead aRead, int offset, long position, char mismatchBase, int qualityScore) {
            this.aRead = aRead;
            this.offset = offset;
            this.position = position;
            this.mismatchBase = mismatchBase;
            this.qualityScore = qualityScore;
        }
    }

    private static int mismatchQualitySum(AlignedRead aRead, String refSeq, int refIndex) {
        return mismatchQualitySum(aRead, refSeq, refIndex, 0);
    }

    private static int mismatchQualitySum(AlignedRead aRead, String refSeq, int refIndex, int minMismatchQualityScore) {
        List<Mismatch> mismatches = getMismatches(aRead, refSeq, refIndex);
        int sum = 0;
        for(Mismatch mm : mismatches) {
            if (mm.qualityScore >= minMismatchQualityScore) {
                sum += mm.qualityScore;
            }
        }
        return sum;

    }




    /**
     * Returns a list of position
     * @param aRead
     * @param refSeq
     * @param refIndex
     * @param minMismatchQualityScore
     * @return
     */
    private static List<Mismatch> getMismatches(AlignedRead aRead, String refSeq, int refIndex) {
        List<Mismatch> mismatches = new ArrayList<Mismatch>();

        String readSeq = aRead.getReadString();
        String quals = aRead.getBaseQualityString();
        int readIndex = 0;
        int sum = 0;
        Cigar c = aRead.getCigar();
        for (int i = 0 ; i < c.numCigarElements() ; i++) {
            CigarElement ce = c.getCigarElement(i);
            switch ( ce.getOperator() ) {
                case M:
                    for (int j = 0 ; j < ce.getLength() ; j++, refIndex++, readIndex++ ) {
                        // FIXME: what is this case????
                        if ( refIndex >= refSeq.length() ) {
                            sum += MAX_QUAL;
                        } else {
                            char readBase = Character.toUpperCase(readSeq.charAt(readIndex));
                            char refBase = Character.toUpperCase(refSeq.charAt(refIndex));

                            if ( readBase != refBase ) {
                                int qual = (int)quals.charAt(readIndex) - 33;
                                mismatches.add(new Mismatch(aRead, readIndex, refIndex, readBase, qual));
                            }
                        }
                    }
                    break;
                case I:
                    readIndex += ce.getLength();
                    break;
                case D:
                    refIndex += ce.getLength();
                    break;
            }

        }
        return mismatches;
    }

    private Pair<Integer, Integer> findBestOffset(String ref, AlignedRead read) {
        int attempts = ref.length() - read.getReadLength() + 1;
        int bestScore = mismatchQualitySum(read, ref, 0);
        int bestIndex = 0;
        for ( int i = 1; i < attempts; i++ ) {
            // we can't get better than 0!
            if ( bestScore == 0 )
                return new Pair<Integer, Integer>(bestIndex, 0);
            int score = mismatchQualitySum(read, ref, i);
            if ( score < bestScore ) {
                bestScore = score;
                bestIndex = i;
            }
        }
        return new Pair<Integer, Integer>(bestIndex, bestScore);
    }
    

    private class AlignedRead {
        SAMRecord read;
        int mismatchScoreToReference;

        public AlignedRead(SAMRecord read) {
            this.read = read;
            mismatchScoreToReference = 0;
        }

        public SAMRecord getRead() {
               return read;
        }

        public String getReadString() {
               return read.getReadString();
        }

        public int getReadLength() {
               return read.getReadLength();
        }

        public Cigar getCigar() {
            return read.getCigar();
        }

        public void setCigar(Cigar cigar) {
            read.setCigar(cigar);
        }

        public String getBaseQualityString() {
            return read.getBaseQualityString();
        }

        public void setMismatchScoreToReference(int score) {
            mismatchScoreToReference = score;
        }

        public int getMismatchScoreToReference() {
            return mismatchScoreToReference;
        }
    }

    private class Consensus {
        public String str;
        public int mismatchSum;
        public int positionOnReference;
        public Cigar cigar;
        public ArrayList<Pair<Integer, Integer>> readIndexes;

        public Consensus(String str, Cigar cigar, int positionOnReference) {
            this.str = str;
            this.cigar = cigar;
            this.positionOnReference = positionOnReference;
            mismatchSum = 0;
            readIndexes = new ArrayList<Pair<Integer, Integer>>();
        }
    }

    public static String cigarToString(Cigar cig) {
        StringBuilder b = new StringBuilder();

        for ( int i = 0 ; i < cig.numCigarElements() ; i++ ) {
            char c='?';
            switch ( cig.getCigarElement(i).getOperator() ) {
                case M : c = 'M'; break;
                case D : c = 'D'; break;
                case I : c = 'I'; break;
            }
            b.append(cig.getCigarElement(i).getLength());
            b.append(c);
        }
        return b.toString();
    }

}
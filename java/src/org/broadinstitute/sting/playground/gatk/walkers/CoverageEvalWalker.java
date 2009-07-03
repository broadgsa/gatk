package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.utils.IndelLikelihood;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;

import java.util.List;
import java.util.ArrayList;
import java.io.PrintStream;
import java.io.FileNotFoundException;
import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Jun 30, 2009
 * Time: 12:38:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class CoverageEvalWalker extends LocusWalker<List<String>, String> {

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName="lod_threshold", shortName="lod", doc="The lod threshold on which variants should be filtered", required=false) public Double LOD_THRESHOLD = 5.0;
    @Argument(fullName="format_geli", shortName="geli", doc="Output variant calls in Geli/Picard format", required=false) public boolean GELI_OUTPUT_FORMAT = false;

    @Argument(fullName="variants_out", shortName="varout", doc="File to which variants should be written", required=true) public File VARIANTS_FILE;

    public PrintStream variantsOut;

    public void initialize() {
        try {
            variantsOut = new PrintStream(VARIANTS_FILE);
        } catch (FileNotFoundException e) {
            err.format("Unable to open file '%s'. Perhaps the parent directory does not exist or is read-only.\n", VARIANTS_FILE.getAbsolutePath());
            System.exit(-1);
        }
                
        String header = GELI_OUTPUT_FORMAT ? AlleleFrequencyEstimate.geliHeaderString() : AlleleFrequencyEstimate.asTabularStringHeader();
        variantsOut.println("#DownsampledCoverage\tAvailableCoveragt \t"+header);
    }

    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref) != -1 && context.getReads().size() != 0);
    }

    public List<String> map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        rodGFF hapmap_chip = (rodGFF)tracker.lookup("hapmap-chip", null);
        String hc_genotype;
        if (hapmap_chip != null) {
            hc_genotype = hapmap_chip.getFeature();
        }else{
            hc_genotype = new String(new char[] {ref, ref});
        }

        //if (tracker.hasROD("hapmap-chip")) {
        ArrayList<String> Gs = new ArrayList<String>();

        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        // Iterate over coverage levels
        int coverage_available = reads.size();
        int coverage_levels[] = {4, 10, 20, Integer.MAX_VALUE};
        int downsampling_repeats = 10; // number of times to random re-sample each coverage_level
        for (int coverage : coverage_levels) {
            coverage = Math.min(coverage_available, coverage); // don't exceed max available coverage
            for (int r=0; r<downsampling_repeats; r++) {
                List<Integer> subset_indices = ListUtils.randomSubsetIndices(coverage, coverage_available);
                List<SAMRecord> sub_reads = ListUtils.subsetListByIndices(subset_indices, reads);
                List<Integer> sub_offsets = ListUtils.subsetListByIndices(subset_indices, offsets);

                // Call genotypes on subset of reads and offsets
                GenotypeLikelihoods G = callGenotype(tracker, ref, pileup, sub_reads, sub_offsets);
                String geliString = G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods, "sample").asGeliString();

                Gs.add(hc_genotype+"\t"+coverage+"\t"+coverage_available+"\t"+geliString);
            }
        }

        return Gs;
    }

    /**
     * Calls the underlying, single locus genotype of the sample
     *
     * @param tracker  the meta data tracker
     * @param ref      the reference base
     * @param pileup   the pileup object for the given locus
     * @param reads    the reads that overlap this locus
     * @param offsets  the offsets per read that identify the base at this locus
     * @return the likelihoods per genotype
     */
    private GenotypeLikelihoods callGenotype(RefMetaDataTracker tracker, char ref, ReadBackedPileup pileup, List<SAMRecord> reads, List<Integer> offsets) {
        GenotypeLikelihoods G;

        G = new GenotypeLikelihoods(); 

        for ( int i = 0; i < reads.size(); i++ ) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);

            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }

        G.ApplyPrior(ref, 'N', -1);

        return G;
    }

    public String reduceInit() {
        return "";
    }

    public String reduce(List<String> alleleFreqLines, String sum) {
        //GenomeLoc a =  GenomeLocParser.parseGenomeLoc("chr1:42971309");
        //if ((alleleFreq != null && alleleFreq.lodVsRef >= LOD_THRESHOLD)) { // || (alleleFreq.location == a) )   {
        for (String line : alleleFreqLines) {
            variantsOut.println(line);
        }

        return "";
	}
}






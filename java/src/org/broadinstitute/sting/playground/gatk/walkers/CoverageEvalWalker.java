package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ListUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.GenotypeCall;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

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
    @Argument(fullName="min_coverage", shortName="mincov", doc="Mininum coverage to downsample to", required=false) public int min_coverage=1;
    @Argument(fullName="max_coverage", shortName="maxcov", doc="Maximum coverage to downsample to", required=false) public int max_coverage=20;
    @Argument(fullName="downsampling_repeats", shortName="repeat", doc="Number of times to repeat downsampling at each coverage level", required=false) public int downsampling_repeats=20;

    public PrintStream variantsOut;

    SingleSampleGenotyper SSG;

    public void initialize() {
        SSG = new SingleSampleGenotyper();
        SSG.VARIANTS_FILE = VARIANTS_FILE;
        SSG.initialize();

        try {
            variantsOut = new PrintStream(VARIANTS_FILE);
        } catch (FileNotFoundException e) {
            err.format("Unable to open file '%s'. Perhaps the parent directory does not exist or is read-only.\n", VARIANTS_FILE.getAbsolutePath());
            System.exit(-1);
        }
                
        String header = "#Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod dbSNP   AA      AC      AG      AT      CC      CG      CT      GG      GT      TT";
        variantsOut.println("DownsampledCoverage\tAvailableCoverage\tHapmapChipGenotype\tGenotypeCallType\t"+header.substring(1));
    }

    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref) != -1 && context.getReads().size() != 0);
    }

    public List<String> map(RefMetaDataTracker tracker, char ref, LocusContext context) {

        rodGFF hapmap_chip = (rodGFF)tracker.lookup("hapmap-chip", null);
        String hc_genotype;

        if (hapmap_chip != null) {
            hc_genotype = hapmap_chip.getFeature();

            ArrayList<String> GenotypeCalls = new ArrayList<String>();

            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();

            int coverage_available = reads.size();
            List<Integer> coverage_levels = new ArrayList<Integer>();// = {4, 7, 10, 20, Integer.MAX_VALUE};
            for (int coverage = min_coverage; coverage <= max_coverage; coverage++) {
                coverage_levels.add(coverage);
            }
            coverage_levels.add(coverage_available); // Run on all available reads

            // Iterate over coverage levels
            for (int coverage : coverage_levels) {
                coverage = Math.min(coverage_available, coverage); // don't exceed max available coverage
                for (int r=0; r<downsampling_repeats; r++) {
                    List<Integer> subset_indices = ListUtils.sampleIndicesWithReplacement(coverage, coverage_available);
                    List<SAMRecord> sub_reads = ListUtils.sliceListByIndices(subset_indices, reads);
                    List<Integer> sub_offsets = ListUtils.sliceListByIndices(subset_indices, offsets);

                    LocusContext subContext = new LocusContext(context.getLocation(), sub_reads, sub_offsets);
                    GenotypeCall call = SSG.map(tracker, ref, subContext);

                    String callType = (call.isVariation()) ? ((call.getBestVrsRef().first.isHom()) ? "HomozygousSNP" : "HeterozygousSNP") : "HomozygousReference";
                    if (call != null) {
                        GenotypeCalls.add(coverage+"\t"+coverage_available+"\t"+hc_genotype+"\t"+callType+"\t"+toGeliString(call));
                    }
                }
            }
            return GenotypeCalls;
        }else{
            return new ArrayList<String>();
        }
    }

    public String reduceInit() {
        return "";
    }

    public String reduce(List<String> alleleFreqLines, String sum) {
        for (String line : alleleFreqLines) {
            variantsOut.println(line);
        }

        return "";
	}

    // a method to support getting the geli string, since the AlleleFrequencyEstimate is going away
    public String toGeliString (GenotypeCall locus) {
        if (locus.getPosteriors().size() != 10) throw new IllegalArgumentException("Geli text only supports SNP calls, with a diploid organism (i.e. posterior array size of 10)");


        // this is to perserve the format string that we used to use
        double[] likelihoods = new double[10];
        int index = 0;
        List<Genotype> lt = locus.getLexigraphicallySortedGenotypes();
        for (Genotype G: lt) {
            likelihoods[index] = G.getLikelihood();
            index++;
        }

        return String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
	                                        locus.getLocation().getContig(),
                                            locus.getLocation().getStart(),
											locus.getReferencebase(),
                                            locus.getReadDepth(),
                                            -1,
	                                        locus.getGenotypes().get(0).getBases(),
	                                        locus.getBestVrsRef().second.getScore(),
	                                        locus.getBestVrsNext().second.getScore(),
                                            likelihoods[0],
                                            likelihoods[1],
                                            likelihoods[2],
                                            likelihoods[3],
                                            likelihoods[4],
                                            likelihoods[5],
                                            likelihoods[6],
                                            likelihoods[7],
                                            likelihoods[8],
                                            likelihoods[9]);
    }
}






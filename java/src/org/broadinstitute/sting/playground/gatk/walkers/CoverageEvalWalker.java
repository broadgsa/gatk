package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodGenotypeChipAsGFF;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.genotyper.GenotypeCall;
import org.broadinstitute.sting.gatk.walkers.genotyper.UnifiedGenotyper;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.ListUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.ArrayList;
import java.util.Arrays;
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

    @Argument(fullName="min_coverage", shortName="mincov", doc="Mininum coverage to downsample to", required=false) public int min_coverage=1;
    @Argument(fullName="max_coverage", shortName="maxcov", doc="Maximum coverage to downsample to", required=false) public int max_coverage=Integer.MAX_VALUE;
    @Argument(fullName="downsampling_repeats", shortName="repeat", doc="Number of times to repeat downsampling at each coverage level", required=false) public int downsampling_repeats=1;

    UnifiedGenotyper UG;

    public void initialize() {
        UG = new UnifiedGenotyper();
        UG.initializePublic();

        String header = "#Sequence       Position        ReferenceBase   NumberOfReads   MaxMappingQuality       BestGenotype    BtrLod  BtnbLod AA      AC      AG      AT      CC      CG      CT      GG      GT      TT";
        out.println("DownsampledCoverage\tAvailableCoverage\tHapmapChipGenotype\tGenotypeCallType\t"+header.substring(1));
    }

    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 && context.getReads().size() != 0);
    }

    public List<String> map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {

        RodGenotypeChipAsGFF hapmap_chip = (RodGenotypeChipAsGFF)tracker.lookup("hapmap-chip", null);
        String hc_genotype;

        if (hapmap_chip != null) {
            hc_genotype = hapmap_chip.getFeature();

            ArrayList<String> GenotypeCalls = new ArrayList<String>();

            List<SAMRecord> reads = context.getReads();
            List<Integer> offsets = context.getOffsets();

            int coverage_available = reads.size();
            List<Integer> coverage_levels = new ArrayList<Integer>();// = {4, 7, 10, 20, Integer.MAX_VALUE};
            Integer this_max_coverage = Math.min(max_coverage, coverage_available); 
            for (int coverage = min_coverage; coverage <= this_max_coverage; coverage++) {
                coverage_levels.add(coverage);
            }
            //coverage_levels.add(coverage_available); // Run on all available reads

            // Iterate over coverage levels
            for (int coverage : coverage_levels) {
                coverage = Math.min(coverage_available, coverage); // don't exceed max available coverage
                for (int r=0; r<downsampling_repeats; r++) {
                    List<Integer> subset_indices = ListUtils.sampleIndicesWithReplacement(coverage_available, coverage);
                    List<SAMRecord> sub_reads = ListUtils.sliceListByIndices(subset_indices, reads);
                    List<Integer> sub_offsets = ListUtils.sliceListByIndices(subset_indices, offsets);

                    AlignmentContext subContext = new AlignmentContext(context.getLocation(), sub_reads, sub_offsets);
                    Pair<List<GenotypeCall>, GenotypeMetaData> calls = UG.map(tracker, ref, subContext);
                    if (calls != null && calls.first != null) {
                        GenotypeCall call = calls.first.get(0);
                        String callType = (call.isVariant(call.getReference())) ? ((call.isHom()) ? "HomozygousSNP" : "HeterozygousSNP") : "HomozygousReference";
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

    public void onTraversalDone(String result) {} // Don't print the reduce result

    public String reduce(List<String> alleleFreqLines, String sum) {
        for (String line : alleleFreqLines) {
            out.println(line);
        }

        return "";
	}

    // a method to support getting the geli string, since the AlleleFrequencyEstimate is going away
    public String toGeliString (Genotype locus) {
        double posteriors[];
        int readDepth = -1;
        double nextVrsBest = 0;
        double nextVrsRef = 0;

        char ref = locus.getReference();

        if (locus instanceof ReadBacked) {
            readDepth = ((ReadBacked)locus).getReadCount();
        }
        if (!(locus instanceof GenotypesBacked)) {
            posteriors = new double[10];
            Arrays.fill(posteriors, 0.0);
        } else {
            posteriors = ((PosteriorsBacked) locus).getPosteriors();
            double[] lks;
            lks = Arrays.copyOf(posteriors,posteriors.length);
            Arrays.sort(lks);
            nextVrsBest = lks[9] - lks[8];
            if (ref != 'X')  {
                int index = (DiploidGenotype.valueOf(Utils.dupString(ref,2)).ordinal());
                nextVrsRef = lks[9] - posteriors[index];
            }
        }
        // we have to calcuate our own

       return new String(String.format("%s    %16d  %c  %8d  %d  %s %.6f %.6f    %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f %6.6f",
                                      locus.getLocation().getContig(),
                                      locus.getLocation().getStart(),
                                      ref,
                                      readDepth,
                                      -1,
                                      locus.getBases(),
                                      nextVrsRef,
                                      nextVrsBest,
                                      posteriors[0],
                                      posteriors[1],
                                      posteriors[2],
                                      posteriors[3],
                                      posteriors[4],
                                      posteriors[5],
                                      posteriors[6],
                                      posteriors[7],
                                      posteriors[8],
                                      posteriors[9]));
    }
}






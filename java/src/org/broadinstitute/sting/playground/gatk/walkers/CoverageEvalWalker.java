package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.playground.utils.IndelLikelihood;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.playground.utils.AlleleFrequencyEstimate;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.*;

import java.util.List;
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
public class CoverageEvalWalker extends LocusWalker<AlleleFrequencyEstimate, String> {

    // Control how the genotype hypotheses are weighed
    @Argument(fullName="priors_any_locus", shortName="plocus", doc="Comma-separated prior likelihoods for any locus (homref,het,homvar)", required=false) public String PRIORS_ANY_LOCUS = "0.999,1e-3,1e-5";
    @Argument(fullName="priors_hapmap", shortName="phapmap", doc="Comma-separated prior likelihoods for Hapmap loci (homref,het,homvar)", required=false) public String PRIORS_HAPMAP = "0.999,1e-3,1e-5";
    @Argument(fullName="priors_dbsnp", shortName="pdbsnp", doc="Comma-separated prior likelihoods for dbSNP loci (homref,het,homvar)", required=false) public String PRIORS_DBSNP = "0.999,1e-3,1e-5";
    @Argument(fullName="priors_2nd_on", shortName="p2ndon", doc="Comma-separated prior likelihoods for the secondary bases of primary on-genotype bases (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT)", required=false) public String PRIORS_2ND_ON = "0.000,0.302,0.366,0.142,0.000,0.548,0.370,0.000,0.319,0.000";
    @Argument(fullName="priors_2nd_off", shortName="p2ndoff", doc="Comma-separated prior likelihoods for the secondary bases of primary off-genotype bases (AA,AC,AG,AT,CC,CG,CT,GG,GT,TT)", required=false) public String PRIORS_2ND_OFF = "0.480,0.769,0.744,0.538,0.575,0.727,0.768,0.589,0.762,0.505";

    // Control what goes into the variants file and what format that file should have
    @Argument(fullName="lod_threshold", shortName="lod", doc="The lod threshold on which variants should be filtered", required=false) public Double LOD_THRESHOLD = 5.0;
    @Argument(fullName="format_geli", shortName="geli", doc="Output variant calls in Geli/Picard format", required=false) public boolean GELI_OUTPUT_FORMAT = false;

    @Argument(fullName="variants_out", shortName="varout", doc="File to which variants should be written", required=true) public File VARIANTS_FILE;

    public double[] plocus;
    public double[] phapmap;
    public double[] pdbsnp;
    public double[] p2ndon;
    public double[] p2ndoff;

    public PrintStream variantsOut;

    public void initialize() {
        try {
            variantsOut = new PrintStream(VARIANTS_FILE);
        } catch (FileNotFoundException e) {
            err.format("Unable to open file '%s'. Perhaps the parent directory does not exist or is read-only.\n", VARIANTS_FILE.getAbsolutePath());
            System.exit(-1);
        }
                
        String header = GELI_OUTPUT_FORMAT ? AlleleFrequencyEstimate.geliHeaderString() : AlleleFrequencyEstimate.asTabularStringHeader();
        variantsOut.println(header);

        plocus = priorsArray(PRIORS_ANY_LOCUS);
        phapmap = priorsArray(PRIORS_HAPMAP);
        pdbsnp = priorsArray(PRIORS_DBSNP);
        p2ndon = priorsArray(PRIORS_2ND_ON);
        p2ndoff = priorsArray(PRIORS_2ND_OFF);

    }    

    public boolean filter(RefMetaDataTracker tracker, char ref, LocusContext context) {
        return (BaseUtils.simpleBaseToBaseIndex(ref) != -1 && context.getReads().size() != 0);
    }

    private double[] priorsArray(String priorsString) {
        String[] pstrs = priorsString.split(",");
        double[] pdbls = new double[pstrs.length];

        for (int i = 0; i < pstrs.length; i++) {
            pdbls[i] = Double.valueOf(pstrs[i]);
        }

        return pdbls;
    }

    public AlleleFrequencyEstimate map(RefMetaDataTracker tracker, char ref, LocusContext context) {
    

        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();

        // Handle single-base polymorphisms.
        GenotypeLikelihoods G = callGenotype(tracker, ref, pileup, reads, offsets);

        return G.toAlleleFrequencyEstimate(context.getLocation(), ref, bases.length(), bases, G.likelihoods, "sample");
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

        G = new GenotypeLikelihoods(plocus[0], plocus[1], plocus[2], p2ndon, p2ndoff);

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

    public String reduce(AlleleFrequencyEstimate alleleFreq, String sum) {
        GenomeLoc a =  GenomeLocParser.parseGenomeLoc("chr1:42971309");
        if ((alleleFreq != null && alleleFreq.lodVsRef >= LOD_THRESHOLD) || (alleleFreq.location == a) )   {
            String line = GELI_OUTPUT_FORMAT ? alleleFreq.asGeliString() : alleleFreq.asTabularString();
		    variantsOut.println(line);
        }

        return "";
	}
}






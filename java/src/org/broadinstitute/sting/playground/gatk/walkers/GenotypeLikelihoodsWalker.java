package org.broadinstitute.sting.playground.gatk.walkers;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.filters.ZeroMappingQualityReadFilter;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.ReadFilters;
import org.broadinstitute.sting.playground.utils.GenotypeLikelihoods;
import org.broadinstitute.sting.utils.ReadBackedPileup;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.GenotypeWriter;
import org.broadinstitute.sting.utils.genotype.GenotypeWriterFactory;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;

import java.io.File;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class LikelihoodWalker
 *
 * a simple walker to calculate the genotype likelihoods
 */
@ReadFilters(ZeroMappingQualityReadFilter.class)
public class GenotypeLikelihoodsWalker extends LocusWalker<LikelihoodWrapper, GenotypeWriter> {
    @Argument(fullName = "variants_out", shortName = "varout", doc = "File to which variants should be written", required = true) public File VARIANTS_FILE;
    @Argument(fullName = "metrics_out", shortName = "metout", doc = "File to which metrics should be written", required = false) public File METRICS_FILE = new File("/dev/stderr");
    @Argument(fullName = "variant_output_format", shortName = "vf", doc = "File to which metrics should be written", required = false) public GenotypeWriterFactory.GENOTYPE_FORMAT VAR_FORMAT = GenotypeWriterFactory.GENOTYPE_FORMAT.GLF;


    @Override
    public LikelihoodWrapper map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref, context);
        String bases = pileup.getBases();
        List<SAMRecord> reads = context.getReads();
        List<Integer> offsets = context.getOffsets();
        double rmsSum = 0.0;
        // Handle single-base polymorphisms.
        GenotypeLikelihoods G = new GenotypeLikelihoods();
        for (int i = 0; i < reads.size(); i++) {
            SAMRecord read = reads.get(i);
            int offset = offsets.get(i);
            rmsSum += (read.getMappingQuality() * read.getMappingQuality());
            G.add(ref, read.getReadString().charAt(offset), read.getBaseQualities()[offset]);
        }

        // our return
        LikelihoodWrapper wrap = new LikelihoodWrapper();

        return wrapLikelihoods(ref, context, reads, rmsSum, G, wrap);
    }

    /**
     * wrap the likelihood values in with the other data we'll need
     * @param ref the ref base
     * @param context the locus context
     * @param reads the reads
     * @param rmsSum the rms square total (not the actual rms yet)
     * @param g the genotypeLikelihoods
     * @param wrap the object to place the values into
     * @return a likelihood wrapper
     */
    private LikelihoodWrapper wrapLikelihoods(char ref, LocusContext context, List<SAMRecord> reads, double rmsSum, GenotypeLikelihoods g, LikelihoodWrapper wrap) {
        wrap.obj = new LikelihoodObject();
        wrap.obj.setLikelihoodType(LikelihoodObject.LIKELIHOOD_TYPE.LOG);
        for (int x = 0; x < GenotypeLikelihoods.genotypes.length; x++) {
            wrap.obj.setLikelihood(LikelihoodObject.GENOTYPE.valueOf(GenotypeLikelihoods.genotypes[x]), g.likelihoods[x]*10.0);
        }
        wrap.obj.setLikelihoodType(LikelihoodObject.LIKELIHOOD_TYPE.NEGITIVE_LOG);
        wrap.loc = GenomeLocParser.getContigInfo(context.getContig());
        wrap.readDepth = reads.size();
        float rms = (float)(Math.sqrt(rmsSum/reads.size()));
        wrap.rms = (rms > 255) ? 255 : rms;
        wrap.ref = ref;
        wrap.position = context.getLocation().getStart();
        return wrap;
    }

    /**
     * Provide an initial value for reduce computations.
     *
     * @return Initial value of reduce.
     */
    @Override
    public GenotypeWriter reduceInit() {
        return GenotypeWriterFactory.create(VAR_FORMAT, GenomeAnalysisEngine.instance.getEngine().getSAMHeader(),VARIANTS_FILE);    
    }

    /**
     * Reduces a single map with the accumulator provided as the ReduceType.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     *
     * @return accumulator with result of the map taken into account.
     */
    @Override
    public GenotypeWriter reduce(LikelihoodWrapper value, GenotypeWriter sum) {

        sum.addGenotypeCall(value.loc,(int)value.position,value.rms,value.ref,value.readDepth,value.obj);
        return sum;
    }

    /** Close the variant file. */
    public void onTraversalDone(GenotypeWriter result) {
        result.close();
    }
}


class LikelihoodWrapper {
    public LikelihoodObject obj;
    public SAMSequenceRecord loc;
    public long position;
    public float rms;
    public char ref;
    public int readDepth;
}
package org.broadinstitute.sting.playground.gatk.walkers.variants;

import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.*;

@Requires(value={DataSource.READS, DataSource.REFERENCE},referenceMetaData={@RMD(name="hapmap",type=HapMapGenotypeROD.class),@RMD(name="dbsnp",type=rodDbSNP.class),@RMD(name="variant",type=rodVariants.class)})
public class StrandTestPrinter extends LocusWalker<Integer, Integer> {
    @Argument(fullName="concordance_out", shortName="CO", doc="File to which concordant variant stats should be written", required=true)
    File CONCORDANCE_OUT = null;
    @Argument(fullName="sample_name", shortName="sample", doc="Sample name (e.g. NA12878)", required=true)
    String SAMPLE_NAME = null;

    private PrintWriter writer = null;

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        try {
            writer = new PrintWriter(CONCORDANCE_OUT);
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file for writing"));
        }
        writer.println("location\tvariant\tlodBestVsRef\tStrandScore\tForwardCount\tReverseCount\tInHapmap\tIndbSNP\tOurCall\tHapMapCall");
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    /**
     * For each site of interest, rescore the genotype likelihoods by applying the specified feature set.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus is a het site, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, char ref, LocusContext context) {
        HapMapGenotypeROD hapmap = (HapMapGenotypeROD)tracker.lookup("hapmap", null);
        rodDbSNP dbsnp = (rodDbSNP)tracker.lookup("dbsnp", null);
        rodVariants variant = (rodVariants)tracker.lookup("variant", null);

        StringBuffer sb = new StringBuffer();
        char refUpper = Character.toUpperCase(ref);

        boolean realHetVariant = false, realVariant = false, calledHetVariant = false;
        if ( hapmap != null ) {
            String genotype = hapmap.get(SAMPLE_NAME).toUpperCase();
            if ( genotype.charAt(0) != refUpper || genotype.charAt(1) != refUpper ) {
                realVariant = true;
                if ( genotype.charAt(0) != genotype.charAt(1) ) {
                    realHetVariant = true;
                    sb.append(context.getLocation() + "\t" + genotype + "\t");
                }
            }
        }
        if ( variant != null ) {
            String genotype = variant.getBestGenotype().toUpperCase();
            if ( (genotype.charAt(0) != refUpper || genotype.charAt(1) != refUpper) &&
                 (genotype.charAt(0) != genotype.charAt(1)) ) {
                calledHetVariant = true;
                if ( !realHetVariant )
                    sb.append(context.getLocation() + "\t" + genotype + "\t");
            }
        }

        if ( !realHetVariant && !calledHetVariant )
            return 0;

        if ( variant != null )
            sb.append(variant.getLodBtr() + "\t");
        else
            sb.append("-1\t");

        int allele1, allele2;
        if ( hapmap != null ) {
            allele1 = BaseUtils.simpleBaseToBaseIndex(hapmap.get(SAMPLE_NAME).charAt(0));
            allele2 = BaseUtils.simpleBaseToBaseIndex(hapmap.get(SAMPLE_NAME).charAt(1));
        } else {
            allele1 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(0));
            allele2 = BaseUtils.simpleBaseToBaseIndex(variant.getBestGenotype().charAt(1));
        }
        VECFisherStrand.strandTest(ref, context, allele1, allele2, -1, sb);

        sb.append((hapmap != null ? "1" : "0") + "\t" + (dbsnp != null ? "1" : "0") + "\t");
        if ( variant == null)
            sb.append("REF\t");
        else if ( calledHetVariant )
            sb.append("HET\t");
        else
            sb.append("HOMNONREF\t");

        if ( realHetVariant )
            sb.append("HET\n");
        else if ( realVariant )
            sb.append("HOMNONREF\n");
        else
            sb.append("REF\n");

        writer.print(sb.toString());

        return 0;
    }

    /**
     * Increment the number of loci processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of loci processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of hets processed and close out the files.
     *
     * @param result  the number of variants seen.
     */
    public void onTraversalDone(Integer result) {
        out.printf("Processed %d hets.\n", result);

        writer.close();
    }
}
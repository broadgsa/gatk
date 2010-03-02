package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.HapMapGenotypeROD;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import java.io.File;

/**
 * Converts HapMap genotype information to VCF format
 */
public class HapMap2VCF extends RodWalker<Integer, Integer> {

    @Argument(fullName="vcfOutput", shortName="vcf", doc="VCF file to which all variants should be written with annotations", required=true)
    protected File VCF_OUT;

    private VCFWriter vcfWriter;
    String[] sample_ids = null;

    public void initialize() {
        vcfWriter = new VCFWriter(VCF_OUT);
    }

    /**
     * For each HapMap record, generate and print the VCF record string
     * Output the VCF header if this is the first record
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        for ( ReferenceOrderedDatum rod : tracker.getAllRods() ) {
            if ( rod instanceof HapMapGenotypeROD ) {
                HapMapGenotypeROD hapmap_rod = (HapMapGenotypeROD) rod;

                // If this is the first time map is called, we need to fill out the sample_ids from the rod
                if (sample_ids == null) {
                    // ensure that there are no duplicate sample IDs
                    sample_ids = hapmap_rod.getSampleIDs();
                    Set<String> sample_id_set = new LinkedHashSet<String>(Arrays.asList(sample_ids));
                    if (sample_id_set.size() != sample_ids.length)
                        throw new IllegalStateException("Sample set passed into HapMap2VCF has repeated sample IDs");

                    // setup the header fields
                    Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
                    hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
                    hInfo.add(new VCFHeaderLine("source", "HapMap2VCF"));
                    hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

                    // write out the header once
                    vcfWriter.writeHeader(new VCFHeader(hInfo, sample_id_set));
                }

                // Get reference base
                Character ref_allele = ref.getBase();
                VCFGenotypeEncoding refAllele = new VCFGenotypeEncoding(ref_allele.toString());

                // Create a new record
                VCFRecord record = new VCFRecord(Character.toString(ref_allele), context.getLocation(), "GT");

                // Record each sample's genotype info
                String[] hapmap_rod_genotypes = hapmap_rod.getGenotypes();
                for (int i = 0; i < hapmap_rod_genotypes.length; i++) {
                    String sample_id = sample_ids[i];
                    String hapmap_rod_genotype = hapmap_rod_genotypes[i];

                    // for each sample, set the genotype if it exists
                    if (!hapmap_rod_genotype.contains("N")) {
                        List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>();
                        VCFGenotypeEncoding allele1 = new VCFGenotypeEncoding(hapmap_rod_genotype.substring(0,1));
                        VCFGenotypeEncoding allele2 = new VCFGenotypeEncoding(hapmap_rod_genotype.substring(1));
                        alleles.add(allele1);
                        alleles.add(allele2);

                        VCFGenotypeRecord genotype = new VCFGenotypeRecord(sample_id, alleles, VCFGenotypeRecord.PHASE.UNPHASED);
                        record.addGenotypeRecord(genotype);
                        if ( !allele1.equals(refAllele) )
                            record.addAlternateBase(allele1);
                        if ( !allele2.equals(refAllele) )
                            record.addAlternateBase(allele2);
                    }
                }

                // Add dbsnp ID
                rodDbSNP dbsnp = rodDbSNP.getFirstRealSNP(tracker.getTrackData("dbsnp", null));
                if (dbsnp != null)
                    record.setID(dbsnp.getRS_ID());

                // Write the VCF record
                vcfWriter.addRecord(record);
            }
        }

        return 1;
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() {
        return 0;
    }

    /**
     * Increment the number of rods processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of rods processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    public void onTraversalDone(Integer value) {}

}

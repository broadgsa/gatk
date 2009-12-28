package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.HapMapGenotypeROD;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.RodWalker;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Converts HapMap genotype information to VCF format
 */
public class HapMap2VCF extends RodWalker<Integer, Integer> {

    String[] sample_ids;

    /**
     * For each HapMap record, generate and print the VCF record string
     * Output the VCF header if this is the first record
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Iterator<ReferenceOrderedDatum> rods = tracker.getAllRods().iterator();
        boolean first_hapmap_rod = true;
        while ( rods.hasNext() ) {
            ReferenceOrderedDatum rod = rods.next();
            if ( rod instanceof HapMapGenotypeROD ) {
                HapMapGenotypeROD hapmap_rod = (HapMapGenotypeROD) rod;
                if (first_hapmap_rod) {
                    out.println(VCFHeaderString(hapmap_rod.getSampleIDs()));
                    first_hapmap_rod = false;
                }

                // Get reference and alternate bases

                Character alt_allele, ref_allele;

                Set<Character> observed_alleles = getObservedAlleles (hapmap_rod.getGenotypes());
                if (observed_alleles.contains('N'))
                    observed_alleles.remove('N');
                ref_allele = ref.getBase();
                if (observed_alleles.contains(ref_allele))
                    observed_alleles.remove(ref_allele);
                if (observed_alleles.isEmpty()) {
                    alt_allele = ref_allele; // ## todo: Confirm that alt allele becomes ref allele 
                }else{
                    if (observed_alleles.size() != 1) {
                        out.println("Error: more than 2 alleles found in hapmap chip position");
                        alt_allele = 'N';
                        System.exit(1);
                    }else{
                        alt_allele = observed_alleles.iterator().next();
                    }
                }

                // Print all position specific info
                String vcf = hapmap_rod.get("chrom")+"\t"+
                             hapmap_rod.get("pos")+"\t"+
                             hapmap_rod.get("rs#")+"\t"+
                             ref_allele+"\t"+
                             alt_allele+"\t"+
                             "99\t0\t.\tGT:GQ";
                out.print(vcf);

                // Print each sample's genotype info

                String allele_strings = "";
                for (String genotype : hapmap_rod.getGenotypes()) {
                    String allele_str = ""; // one allele string
                    Integer GQ = 99;
                    for (Character allele : genotype.toCharArray()) {
                        if (allele == ref_allele) {
                            allele_str += "0";
                        }else{
                            if (allele == alt_allele) {
                                allele_str += "1";
                            }else{
                                if (allele == 'N') {
                                    allele_str += "0";
                                    GQ = 0;
                                }else{
                                    out.println("ERROR: Unexpected tri-allelic site detected");
                                    System.exit(1);
                                }
                            }
                        }
                        if (allele_str.length() == 1) {
                            allele_str += "/";
                        }
                    }
                    allele_strings += "\t" + allele_str+":"+GQ;
                }

                out.println(allele_strings);
            }
        }

        return 1;
    }

    public static Set<Character> getObservedAlleles (String[] genotypes) {
        Set<Character> observed_alleles = new HashSet<Character>();
        for (String genotype : genotypes)
            for (Character allele : genotype.toCharArray())
                if (!observed_alleles.contains(allele))
                    observed_alleles.add(allele);

        return observed_alleles;
    }

    public String VCFHeaderString (String[] sample_ids) {
        String header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (String sample_id : sample_ids)
            header += "\t" + sample_id;

        return header;
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() {
        out.println("#fileformat=VCFv3.3");
        out.println("##reference=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta");
        out.println("##source=VariantsToVCF");

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

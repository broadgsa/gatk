package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeEncoding;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;

import java.util.*;


public class VariantContextAdaptors {
    public static boolean canBeConvertedToVariantContext(String name, Object variantContainingObject) {
        return convertToVariantContext(name, variantContainingObject) != null;
    }

    public static VariantContext convertToVariantContext(String name, Object variantContainingObject) {
        if ( variantContainingObject instanceof rodDbSNP )
            return dbsnpToVariantContext(name, (rodDbSNP)variantContainingObject);
        else if ( variantContainingObject instanceof RodVCF )
            return vcfToVariantContext(name, ((RodVCF)variantContainingObject).getRecord());
        else if ( variantContainingObject instanceof VCFRecord )
            return vcfToVariantContext(name, (VCFRecord)variantContainingObject);
        else
            return null;
            //throw new IllegalArgumentException("Cannot convert object " + variantContainingObject + " of class " + variantContainingObject.getClass() + " to a variant context");

    }

    private static VariantContext dbsnpToVariantContext(String name, rodDbSNP dbsnp) {
        if ( dbsnp.isSNP() || dbsnp.isIndel() || dbsnp.varType.contains("mixed") ) {
            // add the reference allele
            if ( ! Allele.acceptableAlleleBases(dbsnp.getReference()) ) {
                //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                return null;
            }

            List<Allele> alleles = new ArrayList<Allele>();
            Allele refAllele = new Allele(dbsnp.getReference(), true);
            alleles.add(refAllele);

            // add all of the alt alleles
            for ( String alt : dbsnp.getAlternateAlleleList() ) {
                if ( ! Allele.acceptableAlleleBases(alt) ) {
                    //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                    return null;
                }
                alleles.add(new Allele(alt, false));
            }

            VariantContext vc = new VariantContext(name, dbsnp.getLocation(), alleles);
            vc.validate();
            return vc;
        } else
            return null; // can't handle anything else
    }

    private static VariantContext vcfToVariantContext(String name, VCFRecord vcf) {
        if ( vcf.isSNP() || vcf.isIndel() ) {
            // add the reference allele
            if ( ! Allele.acceptableAlleleBases(vcf.getReference()) ) {
                System.out.printf("Excluding vcf record %s%n", vcf);
                return null;
            }

            Set<String> filters = vcf.isFiltered() ? new HashSet<String>(Arrays.asList(vcf.getFilteringCodes())) : null;
            Map<String, String> attributes = vcf.getInfoValues();
            attributes.put("ID", vcf.getID());

            // add all of the alt alleles
            List<Allele> alleles = new ArrayList<Allele>();
            Allele refAllele = new Allele(vcf.getReference(), true);
            alleles.add(refAllele);
            for ( String alt : vcf.getAlternateAlleleList() ) {
                if ( ! Allele.acceptableAlleleBases(alt) ) {
                    System.out.printf("Excluding vcf record %s%n", vcf);
                    return null;
                }
                alleles.add(new Allele(alt, false));
            }

            Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
            for ( VCFGenotypeRecord vcfG : vcf.getVCFGenotypeRecords() ) {
                List<Allele> genotypeAlleles = new ArrayList<Allele>();
                for ( VCFGenotypeEncoding s : vcfG.getAlleles() )
                    genotypeAlleles.add(Allele.getMatchingAllele(alleles, s.getBases()));

                double pError = vcfG.getNegLog10PError() == VCFGenotypeRecord.MISSING_GENOTYPE_QUALITY ? InferredGeneticContext.NO_NEG_LOG_10PERROR : vcfG.getNegLog10PError();

                Map<String, String> fields = new HashMap<String, String>();
                for ( Map.Entry<String, String> e : vcfG.getFields().entrySet() ) {
                    // todo -- fixme if we put GQ and GF into key itself
                    if ( ! e.getKey().equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) && ! e.getKey().equals(VCFGenotypeRecord.GENOTYPE_FILTER_KEY) )
                        fields.put(e.getKey(), e.getValue());
                }

                Set<String> genotypeFilters = new HashSet<String>();
                if ( vcfG.isFiltered() ) // setup the FL genotype filter fields
                    genotypeFilters.addAll(Arrays.asList(vcfG.getFields().get(VCFGenotypeRecord.GENOTYPE_FILTER_KEY).split(";")));

                Genotype g = new Genotype(vcfG.getSampleName(), genotypeAlleles, pError, genotypeFilters, fields);
                genotypes.put(g.getSampleName(), g);
            }

            VariantContext vc = new VariantContext(name, vcf.getLocation(), alleles, genotypes, vcf.getNegLog10PError(), filters, attributes);
            vc.validate();
            return vc;
        } else
            return null; // can't handle anything else
    }
}
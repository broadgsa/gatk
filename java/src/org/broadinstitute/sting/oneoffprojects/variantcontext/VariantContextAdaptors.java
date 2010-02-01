package org.broadinstitute.sting.oneoffprojects.variantcontext;

import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeEncoding;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;

import java.util.*;


public class VariantContextAdaptors {
    public static boolean canBeConvertedToVariantContext(Object variantContainingObject) {
        return convertToVariantContext(variantContainingObject) != null;
    }

    public static VariantContext convertToVariantContext(Object variantContainingObject) {
        if ( variantContainingObject instanceof rodDbSNP )
            return dbsnpToVariantContext((rodDbSNP)variantContainingObject);
        else if ( variantContainingObject instanceof RodVCF )
            return vcfToVariantContext(((RodVCF)variantContainingObject).getRecord());
        else if ( variantContainingObject instanceof VCFRecord )
            return vcfToVariantContext((VCFRecord)variantContainingObject);
        else
            return null;
            //throw new IllegalArgumentException("Cannot convert object " + variantContainingObject + " of class " + variantContainingObject.getClass() + " to a variant context");

    }

    private static VariantContext dbsnpToVariantContext(rodDbSNP dbsnp) {
        if ( dbsnp.isSNP() || dbsnp.isIndel() || dbsnp.varType.contains("mixed") ) {
            VariantContext vc = new VariantContext(dbsnp.getLocation());

            // add the reference allele
            if ( ! Allele.acceptableAlleleBases(dbsnp.getReference()) ) {
                //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                return null;
            }

            Allele refAllele = new Allele(dbsnp.getReference(), true);
            vc.addAllele(refAllele);

            // add all of the alt alleles
            for ( String alt : dbsnp.getAlternateAlleleList() ) {
                if ( ! Allele.acceptableAlleleBases(alt) ) {
                    //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                    return null;
                }
                vc.addAllele(new Allele(alt, false));
            }

            vc.validate();
            return vc;
        } else
            return null; // can't handle anything else
    }

    private static VariantContext vcfToVariantContext(VCFRecord vcf) {
        if ( vcf.isSNP() || vcf.isIndel() ) {
            VariantContext vc = new VariantContext(vcf.getLocation());

            // add the reference allele
            if ( ! Allele.acceptableAlleleBases(vcf.getReference()) ) {
                System.out.printf("Excluding vcf record %s%n", vcf);
                return null;
            }

            Allele refAllele = new Allele(vcf.getReference(), true);
            vc.addAllele(refAllele);
            vc.setNegLog10PError(vcf.getNegLog10PError());
            vc.setAttributes(vcf.getInfoValues());
            vc.putAttribute("ID", vcf.getID());
            if ( vcf.isFiltered() ) vc.setFilters(Arrays.asList(vcf.getFilteringCodes()));

            // add all of the alt alleles
            for ( String alt : vcf.getAlternateAlleleList() ) {
                if ( ! Allele.acceptableAlleleBases(alt) ) {
                    System.out.printf("Excluding vcf record %s%n", vcf);
                    return null;
                }
                vc.addAllele(new Allele(alt, false));
            }

            for ( VCFGenotypeRecord vcfG : vcf.getVCFGenotypeRecords() ) {
                List<String> alleleStrings = new ArrayList<String>();
                for ( VCFGenotypeEncoding s : vcfG.getAlleles() )
                    alleleStrings.add(s.getBases());

                double pError = vcfG.getNegLog10PError() == VCFGenotypeRecord.MISSING_GENOTYPE_QUALITY ? AttributedObject.NO_NEG_LOG_10PERROR : vcfG.getNegLog10PError();
                Genotype g = new Genotype(vc, alleleStrings, vcfG.getSampleName(), pError);

                for ( Map.Entry<String, String> e : vcfG.getFields().entrySet() ) {
                    if ( ! e.getKey().equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) && ! e.getKey().equals(VCFGenotypeRecord.GENOTYPE_FILTER_KEY) )
                        g.putAttribute(e.getKey(), e.getValue());
                }

                if ( vcfG.isFiltered() ) // setup the FL genotype filter fields
                    g.setFilters(Arrays.asList(vcfG.getFields().get(VCFGenotypeRecord.GENOTYPE_FILTER_KEY.split(";"))));

                vc.addGenotype(g);
            }

            vc.validate();
            return vc;
        } else
            return null; // can't handle anything else
    }
}
package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeEncoding;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;

import java.util.*;

/**
 * A terrible but temporary approach to converting objects to VariantContexts.  If you want to add a converter,
 * you need to create a adaptor object here and register a converter from your class to this object.  When tribble arrives,
 * we'll use a better approach.
 *
 * @author depristo@broadinstitute.org
 */
public class VariantContextAdaptors {
    // --------------------------------------------------------------------------------------------------------------
    //
    // Generic support routines.  Do not modify
    //
    // --------------------------------------------------------------------------------------------------------------

    private static Map<Class, VCAdaptor> adaptors = new HashMap<Class, VCAdaptor>();

    static {
        adaptors.put(rodDbSNP.class, new RodDBSnpAdaptor());
        adaptors.put(RodVCF.class, new RodVCFAdaptor());
        adaptors.put(VCFRecord.class, new VCFRecordAdaptor());
    }

    public static boolean canBeConvertedToVariantContext(Object variantContainingObject) {
        return adaptors.containsKey(variantContainingObject.getClass());
//        return convertToVariantContext(name, variantContainingObject) != null;
    }

    /** generic superclass */
    private static abstract class VCAdaptor {
        abstract VariantContext convert(String name, Object input);
    }

    public static VariantContext toVariantContext(String name, Object variantContainingObject) {
        if ( ! adaptors.containsKey(variantContainingObject.getClass()) )
            return null;
        else {
            return adaptors.get(variantContainingObject.getClass()).convert(name, variantContainingObject);
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // From here below you can add adaptor classes for new rods (or other types) to convert to VCF
    //
    // --------------------------------------------------------------------------------------------------------------



    // --------------------------------------------------------------------------------------------------------------
    //
    // dbSNP to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class RodDBSnpAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            rodDbSNP dbsnp = (rodDbSNP)input;
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
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // VCF to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class RodVCFAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            return vcfToVariantContext(name, ((RodVCF)input).getRecord());
        }
    }

    private static class VCFRecordAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            return vcfToVariantContext(name, (VCFRecord)input);
        }
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
                    // todo -- cleanup
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

                double pError = vcfG.getNegLog10PError() == VCFGenotypeRecord.MISSING_GENOTYPE_QUALITY ? Genotype.NO_NEG_LOG_10PERROR : vcfG.getNegLog10PError();

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
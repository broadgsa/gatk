package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;

import java.util.*;

/**
 * A terrible but temporary approach to converting objects to VariantContexts.  If you want to add a converter,
 * you need to create a adaptor object here and register a converter from your class to this object.  When tribble arrives,
 * we'll use a better approach.
 *
 * To add a new converter:
 *
 *   create a subclass of VCAdaptor, overloading the convert operator
 *   add it to the static map from input type -> converter where the input type is the object.class you want to convert
 *
 * That's it 
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

                Map<String, String> attributes = new HashMap<String, String>();
                attributes.put("ID", dbsnp.getRS_ID());
                Collection<Genotype> genotypes = null;
                VariantContext vc = new VariantContext(name, dbsnp.getLocation(), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, attributes);
                vc.validate();
                return vc;
            } else
                return null; // can't handle anything else
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // VCF to VariantContext and back again
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
                for ( VCFGenotypeEncoding s : vcfG.getAlleles() ) {
                    Allele a = Allele.getMatchingAllele(alleles, s.getBases());
                    if ( a == null )
                        throw new StingException("Invalid VCF genotype allele " + s + " in VCF " + vcf);

                    genotypeAlleles.add(a);
                }

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

                Genotype g = new Genotype(vcfG.getSampleName(), genotypeAlleles, pError, genotypeFilters, fields, vcfG.getPhaseType() == VCFGenotypeRecord.PHASE.PHASED);
                genotypes.put(g.getSampleName(), g);
            }

            double qual = vcf.getQual() == -1 ? VariantContext.NO_NEG_LOG_10PERROR : vcf.getNegLog10PError();
            VariantContext vc = new VariantContext(name, vcf.getLocation(), alleles, genotypes, qual, filters, attributes);
            vc.validate();
            return vc;
        } else
            return null; // can't handle anything else
    }


    public static VCFHeader createVCFHeader(Set<VCFHeaderLine> hInfo, VariantContext vc) {
        HashSet<String> names = new LinkedHashSet<String>();
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            names.add(g.getSampleName());
        }

        return new VCFHeader(hInfo == null ? new HashSet<VCFHeaderLine>() : hInfo, names);
    }

    public static VCFRecord toVCF(VariantContext vc) {
        // deal with the reference
        String referenceBase = new String(vc.getReference().getBases());

        String contig = vc.getLocation().getContig();
        long position = vc.getLocation().getStart();
        String ID = vc.hasAttribute("ID") ? vc.getAttributeAsString("ID") : ".";
        double qual = vc.hasNegLog10PError() ? vc.getPhredScaledQual() : -1;
        
        String filters = vc.isFiltered() ? Utils.join(";", Utils.sorted(vc.getFilters())) : VCFRecord.PASSES_FILTERS;

        Map<Allele, VCFGenotypeEncoding> alleleMap = new HashMap<Allele, VCFGenotypeEncoding>();
        alleleMap.put(Allele.NO_CALL, new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE)); // convenience for lookup
        List<VCFGenotypeEncoding> vcfAltAlleles = new ArrayList<VCFGenotypeEncoding>();
        for ( Allele a : vc.getAlleles() ) {
            VCFGenotypeEncoding encoding = new VCFGenotypeEncoding(new String(a.getBases()));
            if ( a.isNonReference() ) {
                vcfAltAlleles.add(encoding);
            }
            alleleMap.put(a, encoding);
        }

        List<String> vcfGenotypeAttributeKeys = new ArrayList<String>(Arrays.asList(VCFGenotypeRecord.GENOTYPE_KEY));
        List<String> vcGenotypeKeys = calcVCFGenotypeKeys(vc);
        if ( vc.hasGenotypes() ) vcfGenotypeAttributeKeys.addAll(vcGenotypeKeys);
        String genotypeFormatString = Utils.join(VCFRecord.GENOTYPE_FIELD_SEPERATOR, vcfGenotypeAttributeKeys);

        List<VCFGenotypeRecord> genotypeObjects = new ArrayList<VCFGenotypeRecord>(vc.getGenotypes().size());
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            List<VCFGenotypeEncoding> encodings = new ArrayList<VCFGenotypeEncoding>(g.getPloidy());
            for ( Allele a : g.getAlleles() ) {
                encodings.add(alleleMap.get(a));
            }

            VCFGenotypeRecord.PHASE phasing = g.genotypesArePhased() ? VCFGenotypeRecord.PHASE.PHASED : VCFGenotypeRecord.PHASE.UNPHASED;
            VCFGenotypeRecord vcfG = new VCFGenotypeRecord(g.getSampleName(), encodings, phasing);

            if ( ! g.isNoCall() ) {
                for ( String key : vcGenotypeKeys ) {
                    String val = key.equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) ? String.format("%.2f", g.getPhredScaledQual()) : g.getAttribute(key).toString();
                    vcfG.setField(key, val);
                }
            }

            genotypeObjects.add(vcfG);
        }

        // info fields
        Map<String, String> infoFields = new HashMap<String, String>();
        for ( Map.Entry<String, Object> elt : vc.getAttributes().entrySet() ) {
            String key = elt.getKey();
            String val = elt.getValue().toString();
            if ( ! key.equals("ID") ) {
                infoFields.put(key, val);
            }
        }

        return new VCFRecord(referenceBase, contig, position, ID, vcfAltAlleles, qual, filters, infoFields, genotypeFormatString, genotypeObjects);
    }

    private static List<String> calcVCFGenotypeKeys(VariantContext vc) {
        Set<String> keys = new HashSet<String>();

        for ( Genotype g : vc.getGenotypes().values() ) {
            for ( String key : g.getAttributes().keySet() ) {
                keys.add(key);
            }
        }

        keys.add(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);
        return Utils.sorted(new ArrayList<String>(keys));
    }
}

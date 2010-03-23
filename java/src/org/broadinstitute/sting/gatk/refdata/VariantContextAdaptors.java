package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableGenotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.glf.GLFSingleCall;
import org.broadinstitute.sting.utils.genotype.glf.GLFWriter;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

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
        adaptors.put(PlinkRod.class, new PlinkRodAdaptor());
        adaptors.put(RodGLF.class, new GLFAdaptor());
        // adaptors.put(RodGeliText.class, new GeliAdaptor());
    }

    public static boolean canBeConvertedToVariantContext(Object variantContainingObject) {
        return adaptors.containsKey(variantContainingObject.getClass());
//        return convertToVariantContext(name, variantContainingObject) != null;
    }

    /** generic superclass */
    private static abstract class VCAdaptor {
        abstract VariantContext convert(String name, Object input);
        abstract VariantContext convert(String name, Object input, Allele refAllele);
    }

    public static VariantContext toVariantContext(String name, Object variantContainingObject) {
        if ( ! adaptors.containsKey(variantContainingObject.getClass()) )
            return null;
        else {
            return adaptors.get(variantContainingObject.getClass()).convert(name, variantContainingObject);
        }
    }

    public static VariantContext toVariantContext(String name, Object variantContainingObject, Allele refAllele) {
        if ( ! adaptors.containsKey(variantContainingObject.getClass()) )
            return null;
        else {
            return adaptors.get(variantContainingObject.getClass()).convert(name, variantContainingObject, refAllele);
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // From here below you can add adaptor classes for new rods (or other types) to convert to VC
    //
    // --------------------------------------------------------------------------------------------------------------



    // --------------------------------------------------------------------------------------------------------------
    //
    // dbSNP to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class RodDBSnpAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            if ( ! Allele.acceptableAlleleBases(((rodDbSNP)input).getReference()) )
                return null;
            Allele refAllele = new Allele(((rodDbSNP)input).getReference(), true);
            return convert(name, input, refAllele);
        }

        VariantContext convert(String name, Object input, Allele refAllele) {
            rodDbSNP dbsnp = (rodDbSNP)input;
            if ( dbsnp.isSNP() || dbsnp.isIndel() || dbsnp.varType.contains("mixed") ) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
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
            Allele refAllele = new Allele(((RodVCF)input).getReference(), true);
            return vcfToVariantContext(name, ((RodVCF)input).getRecord(), refAllele);
        }

        VariantContext convert(String name, Object input, Allele refAllele) {
            return vcfToVariantContext(name, ((RodVCF)input).getRecord(), refAllele);
        }
    }

    private static class VCFRecordAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            Allele refAllele = new Allele(((VCFRecord)input).getReference(), true);
            return vcfToVariantContext(name, (VCFRecord)input, refAllele);
        }

        VariantContext convert(String name, Object input, Allele refAllele) {
            return vcfToVariantContext(name, (VCFRecord)input, refAllele);
        }
    }


    private static VariantContext vcfToVariantContext(String name, VCFRecord vcf, Allele refAllele) {
        if ( vcf.isReference() || vcf.isSNP() || vcf.isIndel() ) {
            // add the reference allele
            if ( ! Allele.acceptableAlleleBases(vcf.getReference()) ) {
                System.out.printf("Excluding vcf record %s%n", vcf);
                return null;
            }

            Set<String> filters = vcf.isFiltered() ? new HashSet<String>(Arrays.asList(vcf.getFilteringCodes())) : null;
            Map<String, String> attributes = new HashMap<String, String>(vcf.getInfoValues());
            attributes.put("ID", vcf.getID());

            // add all of the alt alleles
            List<Allele> alleles = new ArrayList<Allele>();
            alleles.add(refAllele);
            for ( String alt : vcf.getAlternateAlleleList() ) {
                if ( ! Allele.acceptableAlleleBases(alt) ) {
                    //System.out.printf("Excluding vcf record %s%n", vcf);
                    return null;
                }

                Allele allele = new Allele(alt, false);
                if ( ! allele.isNoCall() )
                    alleles.add(allele);
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

                Map<String, String> fields = new HashMap<String, String>();
                for ( Map.Entry<String, String> e : vcfG.getFields().entrySet() ) {
                    // todo -- fixme if we put GQ and GF into key itself
                    if ( ! e.getKey().equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) && ! e.getKey().equals(VCFGenotypeRecord.GENOTYPE_FILTER_KEY) )
                        fields.put(e.getKey(), e.getValue());
                }

                Set<String> genotypeFilters = new HashSet<String>();
                if ( vcfG.isFiltered() ) // setup the FL genotype filter fields
                    genotypeFilters.addAll(Arrays.asList(vcfG.getFields().get(VCFGenotypeRecord.GENOTYPE_FILTER_KEY).split(";")));

                double qual = vcfG.isMissingQual() ? VariantContext.NO_NEG_LOG_10PERROR : vcfG.getNegLog10PError();
                Genotype g = new Genotype(vcfG.getSampleName(), genotypeAlleles, qual, genotypeFilters, fields, vcfG.getPhaseType() == VCFGenotypeRecord.PHASE.PHASED);
                genotypes.put(g.getSampleName(), g);
            }

            double qual = vcf.isMissingQual() ? VariantContext.NO_NEG_LOG_10PERROR : vcf.getNegLog10PError();
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

    public static VCFRecord toVCF(VariantContext vc, char vcfRefBase) {
        return toVCF(vc, vcfRefBase, null, true, false);
    }

    public static VCFRecord toVCF(VariantContext vc, char vcfRefBase, List<String> allowedGenotypeAttributeKeys, boolean filtersWereAppliedToContext, boolean filtersWereAppliedToGenotypes) {
        // deal with the reference
        String referenceBases = new String(vc.getReference().getBases());

        String contig = vc.getLocation().getContig();
        long position = vc.getLocation().getStart();
        String ID = vc.hasAttribute("ID") ? vc.getAttributeAsString("ID") : VCFRecord.EMPTY_ID_FIELD;
        double qual = vc.hasNegLog10PError() ? vc.getPhredScaledQual() : -1;
        
        String filters = vc.isFiltered() ? Utils.join(";", Utils.sorted(vc.getFilters())) : (filtersWereAppliedToContext ? VCFRecord.PASSES_FILTERS : VCFRecord.UNFILTERED);

        Map<Allele, VCFGenotypeEncoding> alleleMap = new HashMap<Allele, VCFGenotypeEncoding>();
        alleleMap.put(Allele.NO_CALL, new VCFGenotypeEncoding(VCFGenotypeRecord.EMPTY_ALLELE)); // convenience for lookup
        List<VCFGenotypeEncoding> vcfAltAlleles = new ArrayList<VCFGenotypeEncoding>();
        for ( Allele a : vc.getAlleles() ) {

            VCFGenotypeEncoding encoding;

            // This is tricky because the VCF spec states that the reference must be a single
            // character, whereas the reference alleles for deletions of length > 1 are strings.
            // To be safe, we require the reference base be passed in and we use that whenever
            // we know that the given allele is the reference.

            String alleleString = new String(a.getBases());
            if ( vc.getType() == VariantContext.Type.MIXED ) {
                throw new UnsupportedOperationException("Conversion from a mixed type isn't currently supported");
            } else if ( vc.getType() == VariantContext.Type.INDEL ) {
                if ( a.isNull() ) {
                    if ( a.isReference() ) {
                        // ref, where alt is insertion
                        encoding = new VCFGenotypeEncoding(Character.toString(vcfRefBase));
                    } else {
                        // non-ref deletion
                        encoding = new VCFGenotypeEncoding("D" + Integer.toString(referenceBases.length()));
                    }
                } else if ( a.isReference() ) {
                    // ref, where alt is deletion
                    encoding = new VCFGenotypeEncoding(Character.toString(vcfRefBase));
                } else {
                    // non-ref insertion
                    encoding = new VCFGenotypeEncoding("I" + alleleString);
                }
            } else if ( vc.getType() == VariantContext.Type.NO_VARIATION ) {
                // ref
                encoding = new VCFGenotypeEncoding(Character.toString(vcfRefBase));
            } else {
                // ref or alt for snp
                encoding = new VCFGenotypeEncoding(alleleString);
            }

            if ( a.isNonReference() ) {
                vcfAltAlleles.add(encoding);
            }

            alleleMap.put(a, encoding);
        }

        List<String> vcfGenotypeAttributeKeys = new ArrayList<String>();
        if ( vc.hasGenotypes() ) {
            vcfGenotypeAttributeKeys.add(VCFGenotypeRecord.GENOTYPE_KEY);
            for ( String key : calcVCFGenotypeKeys(vc) ) {
                if ( allowedGenotypeAttributeKeys == null || allowedGenotypeAttributeKeys.contains(key) )
                    vcfGenotypeAttributeKeys.add(key);
            }
            if ( filtersWereAppliedToGenotypes )
                vcfGenotypeAttributeKeys.add(VCFGenotypeRecord.GENOTYPE_FILTER_KEY);
        }
        String genotypeFormatString = Utils.join(VCFRecord.GENOTYPE_FIELD_SEPERATOR, vcfGenotypeAttributeKeys);

        List<VCFGenotypeRecord> genotypeObjects = new ArrayList<VCFGenotypeRecord>(vc.getGenotypes().size());
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            List<VCFGenotypeEncoding> encodings = new ArrayList<VCFGenotypeEncoding>(g.getPloidy());

            for ( Allele a : g.getAlleles() ) {
                encodings.add(alleleMap.get(a));
            }

            VCFGenotypeRecord.PHASE phasing = g.genotypesArePhased() ? VCFGenotypeRecord.PHASE.PHASED : VCFGenotypeRecord.PHASE.UNPHASED;
            VCFGenotypeRecord vcfG = new VCFGenotypeRecord(g.getSampleName(), encodings, phasing);

            for ( String key : vcfGenotypeAttributeKeys ) {
                if ( key.equals(VCFGenotypeRecord.GENOTYPE_KEY) )
                    continue;
                
                Object val = g.getAttribute(key);
                // some exceptions
                if ( key.equals(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY) ) {
                    if ( MathUtils.compareDoubles(g.getNegLog10PError(), Genotype.NO_NEG_LOG_10PERROR) == 0 )
                        val = VCFGenotypeRecord.MISSING_GENOTYPE_QUALITY;
                    else
                        val = Math.min(g.getPhredScaledQual(), VCFGenotypeRecord.MAX_QUAL_VALUE);
                } else if ( key.equals(VCFGenotypeRecord.DEPTH_KEY) && val == null ) {
                    ReadBackedPileup pileup = (ReadBackedPileup)g.getAttribute(CalledGenotype.READBACKEDPILEUP_ATTRIBUTE_KEY);
                    if ( pileup != null )
                        val = pileup.size();
                } else if ( key.equals(VCFGenotypeRecord.GENOTYPE_FILTER_KEY) ) {
                    val = g.isFiltered() ? Utils.join(";", Utils.sorted(g.getFilters())) : VCFRecord.PASSES_FILTERS;
                }

                String outputValue = formatVCFField(key, val);
                if ( outputValue != null )
                    vcfG.setField(key, outputValue);
            }

            genotypeObjects.add(vcfG);
        }
        // info fields
        Map<String, String> infoFields = new HashMap<String, String>();
        for ( Map.Entry<String, Object> elt : vc.getAttributes().entrySet() ) {
            String key = elt.getKey();
            if ( key.equals("ID") )
                continue;

            String outputValue = formatVCFField(key, elt.getValue());
            if ( outputValue != null )
                infoFields.put(key, outputValue);
        }

        return new VCFRecord(Character.toString(vcfRefBase), contig, position, ID, vcfAltAlleles, qual, filters, infoFields, genotypeFormatString, genotypeObjects);
    }

    private static String formatVCFField(String key, Object val) {
        String result;
        if ( val == null )
            result = VCFGenotypeRecord.getMissingFieldValue(key);
        else if ( val instanceof Double )
            result = String.format("%.2f", val);
        else
            result = val.toString();

        return result;
    }

    private static List<String> calcVCFGenotypeKeys(VariantContext vc) {
        Set<String> keys = new HashSet<String>();

        for ( Genotype g : vc.getGenotypes().values() ) {
            keys.addAll(g.getAttributes().keySet());
        }

        keys.add(VCFGenotypeRecord.GENOTYPE_QUALITY_KEY);
        return Utils.sorted(new ArrayList<String>(keys));
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Plink to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class PlinkRodAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            throw new UnsupportedOperationException("Conversion from Plink to VC requires a reference allele");
        }

        VariantContext convert(String name, Object input, Allele refAllele) {
            PlinkRod plink = (PlinkRod)input;

            HashMap<String, Allele> alleles = new HashMap<String, Allele>(); // use String keys to help maintain uniqueness
            alleles.put(refAllele.isNull() ? Allele.NULL_ALLELE_STRING : new String(refAllele.getBases()), refAllele);

            Set<Genotype> genotypes = new HashSet<Genotype>();

            Map<String, List<String>> genotypeSets = plink.getGenotypes();
            // for each sample
            for ( Map.Entry<String, List<String>> genotype : genotypeSets.entrySet() ) {
                ArrayList<Allele> myAlleles = new ArrayList<Allele>(2);

                // for each allele
                for ( String alleleString : genotype.getValue() ) {
                    Allele allele;
                    if ( alleleString.equals(Allele.NO_CALL_STRING) ) {
                        allele = Allele.NO_CALL;
                    } else {                    
                        if ( !plink.isIndel() ) {
                            allele = new Allele(alleleString, refAllele.basesMatch(alleleString));
                        } else if ( alleleString.equals(Allele.NULL_ALLELE_STRING) ) {
                            allele = new Allele(alleleString, plink.isInsertion());
                        } else {
                            allele = new Allele(alleleString, !plink.isInsertion());
                        }

                        if ( !alleles.containsKey(alleleString) )
                            alleles.put(alleleString, allele);
                    }

                    myAlleles.add(allele);
                }

                // create the genotype
                genotypes.add(new Genotype(genotype.getKey(), myAlleles));
            }

            // create the variant context
            try {
                GenomeLoc loc = GenomeLocParser.setStop(plink.getLocation(), plink.getLocation().getStop() + plink.getLength()-1);
                VariantContext vc = new VariantContext(plink.getVariantName(), loc, alleles.values(), genotypes);
                vc.validate();
                return vc;
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException(e.getMessage() + "; please make sure that e.g. a sample isn't present more than one time in your ped file");
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // GLF to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------
    private static class GLFAdaptor extends VCAdaptor {
        /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGLF
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input) {
            if ( ! Allele.acceptableAlleleBases(((RodGLF)input).getReference()) )
                return null;
            Allele refAllele = new Allele(((RodGLF)input).getReference(), true);
            return convert(name, input, refAllele);
        }

        /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGLF
         * @param refAllele the reference base as an Allele object
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input, Allele refAllele) {
            RodGLF glf = (RodGLF)input;

            // make sure we can convert it
            if ( glf.isSNP() || glf.isIndel()) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                alleles.add(refAllele);

                // add all of the alt alleles
                for ( String alt : glf.getAlternateAlleleList() ) {
                    if ( ! Allele.acceptableAlleleBases(alt) ) {
                        return null;
                    }
                    Allele allele = new Allele(alt, false);
                    if (!alleles.contains(allele)) alleles.add(allele);
                }


                Map<String, String> attributes = new HashMap<String, String>();
                Collection<Genotype> genotypes = new ArrayList<Genotype>();
                MutableGenotype call = new MutableGenotype(name, alleles);

                if (glf.mRecord instanceof GLFSingleCall) {
                    // transform the likelihoods from negative log (positive double values) to log values (negitive values)
                    LikelihoodObject obj = new LikelihoodObject(((GLFSingleCall)glf.mRecord).getLikelihoods(), LikelihoodObject.LIKELIHOOD_TYPE.NEGATIVE_LOG);
                    obj.setLikelihoodType(LikelihoodObject.LIKELIHOOD_TYPE.LOG);

                    // set the likelihoods, depth, and RMS mapping quality values
                    call.putAttribute(CalledGenotype.LIKELIHOODS_ATTRIBUTE_KEY,obj.toDoubleArray());
                    call.putAttribute(VCFGenotypeRecord.DEPTH_KEY,(glf.mRecord.getReadDepth()));
                    call.putAttribute(GLFWriter.RMS_MAPPING_QUAL, (double) glf.mRecord.getRmsMapQ());
                } else {
                    throw new UnsupportedOperationException("We don't currenly support indel calls");
                }

                // add the call to the genotype list, and then use this list to create a VariantContext
                genotypes.add(call);
                VariantContext vc = new VariantContext(name, glf.getLocation(), alleles, genotypes, glf.getNegLog10PError(), null, attributes);
                vc.validate();
                return vc;
            } else
                return null; // can't handle anything else
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // GELI to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------
/*
    private static class GeliAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input) {
            if (!Allele.acceptableAlleleBases(((RodGeliText) input).getReference()))
                return null;
            Allele refAllele = new Allele(((RodGeliText) input).getReference(), true);
            return convert(name, input, refAllele);
        }

        VariantContext convert(String name, Object input, Allele refAllele) {
            RodGeliText geliText = (RodGeliText) input;
            if (geliText.isSNP() || geliText.isIndel()) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                alleles.add(refAllele);

                // add all of the alt alleles
                for (String alt : geliText.getAlternateAlleleList()) {
                    if (!Allele.acceptableAlleleBases(alt)) {
                        return null;
                    }
                    alleles.add(new Allele(alt, false));
                }

                Map<String, String> attributes = new HashMap<String, String>();
                attributes.put("ID", geliText.getName());
                Collection<Genotype> genotypes = null;
                VariantContext vc = new VariantContext(name, geliText.getLocation(), alleles, genotypes, geliText.getNegLog10PError(), null, attributes);
                vc.validate();
                return vc;
            } else
                return null; // can't handle anything else
        }
    }*/
}

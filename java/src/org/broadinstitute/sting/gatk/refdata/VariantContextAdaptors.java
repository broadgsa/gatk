package org.broadinstitute.sting.gatk.refdata;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableGenotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.genotype.LikelihoodObject;
import org.broadinstitute.sting.utils.genotype.geli.GeliTextWriter;
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
        adaptors.put(HapMapROD.class, new HapMapAdaptor());
        adaptors.put(RodGLF.class, new GLFAdaptor());
        adaptors.put(RodGeliText.class, new GeliTextAdaptor());
        adaptors.put(rodGELI.class, new GeliAdaptor());
    }

    public static boolean canBeConvertedToVariantContext(Object variantContainingObject) {
        return adaptors.containsKey(variantContainingObject.getClass());
//        return convertToVariantContext(name, variantContainingObject) != null;
    }

    /** generic superclass */
    private static abstract class VCAdaptor {
        abstract VariantContext convert(String name, Object input);
        abstract VariantContext convert(String name, Object input, ReferenceContext ref);
    }

    public static VariantContext toVariantContext(String name, Object variantContainingObject) {
        if ( ! adaptors.containsKey(variantContainingObject.getClass()) )
            return null;
        else {
            return adaptors.get(variantContainingObject.getClass()).convert(name, variantContainingObject);
        }
    }

    public static VariantContext toVariantContext(String name, Object variantContainingObject, ReferenceContext ref) {
        if ( ! adaptors.containsKey(variantContainingObject.getClass()) )
            return null;
        else {
            return adaptors.get(variantContainingObject.getClass()).convert(name, variantContainingObject, ref);
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
            return convert(name, input, null);
        }

        VariantContext convert(String name, Object input, ReferenceContext ref) {
            rodDbSNP dbsnp = (rodDbSNP)input;
            if ( ! Allele.acceptableAlleleBases(dbsnp.getReference()) )
                return null;
            Allele refAllele = new Allele(dbsnp.getReference(), true);

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
        // WARNING: do not use this method if you have anything other than point mutations in your VCF
        VariantContext convert(String name, Object input) {
            return vcfToVariantContext(name, ((RodVCF)input).getRecord(), null);
        }

        VariantContext convert(String name, Object input, ReferenceContext ref) {
            return vcfToVariantContext(name, ((RodVCF)input).getRecord(), ref);
        }
    }

    private static class VCFRecordAdaptor extends VCAdaptor {
        // WARNING: do not use this method if you have anything other than point mutations in your VCF
        VariantContext convert(String name, Object input) {
            return vcfToVariantContext(name, (VCFRecord)input, null);
        }

        VariantContext convert(String name, Object input, ReferenceContext ref) {
            return vcfToVariantContext(name, (VCFRecord)input, ref);
        }
    }

    private static VariantContext vcfToVariantContext(String name, VCFRecord vcf, ReferenceContext ref) {
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
            Allele refAllele = determineRefAllele(vcf, ref);
            alleles.add(refAllele);

            for ( VCFGenotypeEncoding alt : vcf.getAlternateAlleles() ) {
                if ( ! Allele.acceptableAlleleBases(alt.getBases()) ) {
                    //System.out.printf("Excluding vcf record %s%n", vcf);
                    return null;
                }

                Allele allele;
                // special case: semi-deletion
                if ( vcf.isDeletion() && refAllele.length() > alt.getLength() ) {
                    char[] semiDeletion = new char[refAllele.length() - alt.getLength()];
                    System.arraycopy(ref.getBases(), alt.getLength(), semiDeletion, 0, refAllele.length() - alt.getLength());
                    allele = new Allele(String.valueOf(semiDeletion), false);
                } else {
                    allele = new Allele(alt.getBases(), false);
                }
                if ( ! allele.isNoCall() )
                    alleles.add(allele);
            }

            Map<String, Genotype> genotypes = new HashMap<String, Genotype>();
            for ( VCFGenotypeRecord vcfG : vcf.getVCFGenotypeRecords() ) {
                List<Allele> genotypeAlleles = new ArrayList<Allele>();
                for ( VCFGenotypeEncoding s : vcfG.getAlleles() ) {
                    Allele a = Allele.getMatchingAllele(alleles, s.getBases());
                    if ( a == null ) {
                        if ( vcf.isIndel() )
                            genotypeAlleles.add(refAllele);
                        else
                            throw new StingException("Invalid VCF genotype allele " + s + " in VCF " + vcf);
                    } else {
                        genotypeAlleles.add(a);
                    }
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

            GenomeLoc loc = vcf.getLocation();
            if ( vcf.isDeletion() )
                loc = GenomeLocParser.createGenomeLoc(loc.getContig(), loc.getStart(), loc.getStart()+refAllele.length()-1);

            VariantContext vc = new VariantContext(name, loc, alleles, genotypes, qual, filters, attributes);
            vc.validate();
            return vc;
        } else
            return null; // can't handle anything else
    }

    private static Allele determineRefAllele(VCFRecord vcf, ReferenceContext ref) {
        Allele refAllele;
        if ( vcf.isInsertion() ) {
            refAllele = new Allele(Allele.NULL_ALLELE_STRING, true);
        } else if ( ref == null ) {
            refAllele = new Allele(vcf.getReference(), true);
        } else if ( !vcf.isIndel() ) {
            refAllele = new Allele(Character.toString(ref.getBase()), true);
        } else if ( vcf.isDeletion() ) {
            int start = (int)(ref.getLocus().getStart() - ref.getWindow().getStart() + 1);
            int delLength = 0;
            for ( VCFGenotypeEncoding enc : vcf.getAlternateAlleles() ) {
                if ( enc.getLength() > delLength )
                    delLength = enc.getLength();
            }
            if ( delLength > ref.getWindow().getStop() - ref.getLocus().getStop() )
                throw new IllegalArgumentException("Length of deletion is larger than reference context provided at " + ref.getLocus());

            char[] deletion = new char[delLength];
            System.arraycopy(ref.getBases(), start, deletion, 0, delLength);
            refAllele = new Allele(String.valueOf(deletion), true);
        } else {
            throw new UnsupportedOperationException("Conversion of VCF type " + vcf.getType() + " is not supported.");
        }

        return refAllele;
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

        boolean sawGoodQual = false;
        for ( Genotype g : vc.getGenotypes().values() ) {
            keys.addAll(g.getAttributes().keySet());
            if ( g.hasNegLog10PError() )
                sawGoodQual = true;
        }

        if ( sawGoodQual )
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
            return convert(name, input, null);
        }

        VariantContext convert(String name, Object input, ReferenceContext ref) {
            if ( ref == null )
                throw new UnsupportedOperationException("Conversion from Plink to VariantContext requires a reference context");

            PlinkRod plink = (PlinkRod)input;

            HashSet<Allele> VCAlleles = new HashSet<Allele>();
            Allele refAllele = determineRefAllele(plink, ref);
            VCAlleles.add(refAllele);

            // mapping from Plink Allele to VC Allele (since the PlinkRod does not annotate any of the Alleles as reference)
            HashMap<Allele, Allele> plinkAlleleToVCAllele = new HashMap<Allele, Allele>();

            Set<Genotype> genotypes = new HashSet<Genotype>();

            Map<String, Allele[]> genotypeSets = plink.getGenotypes();

            // We need to iterate through this list and recreate the Alleles since the
            //  PlinkRod does not annotate any of the Alleles as reference.
            // for each sample...
            for ( Map.Entry<String, Allele[]> genotype : genotypeSets.entrySet() ) {
                ArrayList<Allele> myVCAlleles = new ArrayList<Allele>(2);

                // for each allele...
                for ( Allele myPlinkAllele : genotype.getValue() ) {
                    Allele VCAllele;
                    if ( myPlinkAllele.isNoCall() ) {
                        VCAllele = Allele.NO_CALL;
                    } else {                    
                        if ( plinkAlleleToVCAllele.containsKey(myPlinkAllele) ) {
                            VCAllele = plinkAlleleToVCAllele.get(myPlinkAllele);
                        } else {
                            if ( !plink.isIndel() ) {
                                VCAllele = new Allele(myPlinkAllele.getBases(), refAllele.equals(myPlinkAllele, true));
                            } else if ( myPlinkAllele.isNull() ) {
                                VCAllele = new Allele(Allele.NULL_ALLELE_STRING, plink.isInsertion());
                            } else {
                                VCAllele = new Allele(myPlinkAllele.getBases(), !plink.isInsertion());
                            }
                            plinkAlleleToVCAllele.put(myPlinkAllele, VCAllele);
                            VCAlleles.add(VCAllele);
                        }
                    }

                    myVCAlleles.add(VCAllele);
                }

                // create the genotype
                genotypes.add(new Genotype(genotype.getKey(), myVCAlleles));
            }

            // create the variant context
            try {
                GenomeLoc loc = GenomeLocParser.setStop(plink.getLocation(), plink.getLocation().getStop() + plink.getLength()-1);
                VariantContext vc = new VariantContext(plink.getVariantName(), loc, VCAlleles, genotypes);
                vc.validate();
                return vc;
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException(e.getMessage() + "; please make sure that e.g. a sample isn't present more than one time in your ped file");
            }
        }

        private Allele determineRefAllele(PlinkRod plink, ReferenceContext ref) {
            Allele refAllele;
            if ( !plink.isIndel() ) {
                refAllele = new Allele(Character.toString(ref.getBase()), true);
            } else if ( plink.isInsertion() ) {
                refAllele = new Allele(Allele.NULL_ALLELE_STRING, true);
            } else {
                long maxLength = ref.getWindow().getStop() - ref.getLocus().getStop();                
                if ( plink.getLength() > maxLength )
                    throw new UnsupportedOperationException("Plink conversion currently can only handle indels up to length " + maxLength);
                char[] deletion = new char[plink.getLength()];
                System.arraycopy(ref.getBases(), 1, deletion, 0, plink.getLength());
                refAllele = new Allele(String.valueOf(deletion), true);
            }
            return refAllele;
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
            return convert(name, input, null);
        }

        /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGLF
         * @param ref the reference context
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            RodGLF glf = (RodGLF)input;

            if ( ! Allele.acceptableAlleleBases(glf.getReference()) )
                return null;
            Allele refAllele = new Allele(glf.getReference(), true);

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

    private static class GeliTextAdaptor extends VCAdaptor {
          /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input) {
            return convert(name, input, null);
        }

        /**
         * convert to a Variant Context, given:
         * @param name  the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @param ref   the reference context
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            RodGeliText geli = (RodGeliText)input;
            if ( ! Allele.acceptableAlleleBases(geli.getReference()) )
                return null;
            Allele refAllele = new Allele(geli.getReference(), true);

            // make sure we can convert it
            if ( geli.isSNP() || geli.isIndel()) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                alleles.add(refAllele);

                // add all of the alt alleles
                for ( String alt : geli.getAlternateAlleleList() ) {
                    if ( ! Allele.acceptableAlleleBases(alt) ) {
                        return null;
                    }
                    Allele allele = new Allele(alt, false);
                    if (!alleles.contains(allele)) alleles.add(allele);
                }


                Map<String, String> attributes = new HashMap<String, String>();
                Collection<Genotype> genotypes = new ArrayList<Genotype>();
                MutableGenotype call = new MutableGenotype(name, alleles);

                // set the likelihoods, depth, and RMS mapping quality values
                call.putAttribute(CalledGenotype.POSTERIORS_ATTRIBUTE_KEY,geli.genotypePosteriors);
                call.putAttribute(GeliTextWriter.MAXIMUM_MAPPING_QUALITY_ATTRIBUTE_KEY,geli.getMaxMappingQuality());
                call.putAttribute(GeliTextWriter.READ_COUNT_ATTRIBUTE_KEY,geli.depth);

                // add the call to the genotype list, and then use this list to create a VariantContext
                genotypes.add(call);
                VariantContext vc = new VariantContext(name, geli.getLocation(), alleles, genotypes, geli.getNegLog10PError(), null, attributes);
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

    private static class GeliAdaptor extends VCAdaptor {
          /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input) {
            return convert(name, input, null);
        }

        /**
         * convert to a Variant Context, given:
         * @param name  the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @param ref   the reference context
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            GenotypeLikelihoods geli = ((rodGELI)input).getGeliRecord();
            if ( ! Allele.acceptableAlleleBases(String.valueOf((char)geli.getReferenceBase())) )
                return null;
            Allele refAllele = new Allele(String.valueOf((char)geli.getReferenceBase()), true);

            // add the reference allele
            List<Allele> alleles = new ArrayList<Allele>();
            alleles.add(refAllele);

            // add the two alt alleles
            if (!Allele.acceptableAlleleBases(String.valueOf((char) geli.getBestGenotype().getAllele1()))) {
                return null;
            }
            Allele allele1 = new Allele(String.valueOf((char) geli.getBestGenotype().getAllele1()), false);
            if (!alleles.contains(allele1) && !allele1.equals(refAllele, true)) alleles.add(allele1);

            // add the two alt alleles
            if (!Allele.acceptableAlleleBases(String.valueOf((char) geli.getBestGenotype().getAllele2()))) {
                return null;
            }
            Allele allele2 = new Allele(String.valueOf((char) geli.getBestGenotype().getAllele2()), false);
            if (!alleles.contains(allele2) && !allele2.equals(refAllele, true)) alleles.add(allele2);


            Map<String, String> attributes = new HashMap<String, String>();
            Collection<Genotype> genotypes = new ArrayList<Genotype>();
            MutableGenotype call = new MutableGenotype(name, alleles);

            // setup the genotype likelihoods
            double[] post = new double[10];
            String[] gTypes = {"AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT"};
            for (int x = 0; x < 10; x++)
                post[x] = Double.valueOf(geli.getLikelihood(DiploidGenotype.fromBases((byte) gTypes[x].charAt(0), (byte) gTypes[x].charAt(1))));

            // set the likelihoods, depth, and RMS mapping quality values
            call.putAttribute(CalledGenotype.POSTERIORS_ATTRIBUTE_KEY, post);
            call.putAttribute(GeliTextWriter.MAXIMUM_MAPPING_QUALITY_ATTRIBUTE_KEY, (int) geli.getMaxMappingQuality());
            call.putAttribute(GeliTextWriter.READ_COUNT_ATTRIBUTE_KEY, geli.getNumReads());

            // add the call to the genotype list, and then use this list to create a VariantContext
            genotypes.add(call);
            VariantContext vc = new VariantContext(name, ((rodGELI) input).getLocation(), alleles, genotypes, geli.getBestToReferenceLod(), null, attributes);
            vc.validate();
            return vc;

        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // HapMap to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class HapMapAdaptor extends VCAdaptor {
          /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input) {
            return convert(name, input, null);
        }

        /**
         * convert to a Variant Context, given:
         * @param name  the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @param ref   the reference context
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            if ( ref == null )
                throw new UnsupportedOperationException("Conversion from HapMap to VariantContext requires a reference context");

            HapMapROD hapmap = (HapMapROD)input;

            // add the reference allele
            HashSet<Allele> alleles = new HashSet<Allele>();
            Allele refAllele = new Allele(Character.toString(ref.getBase()), true);
            alleles.add(refAllele);

            // make a mapping from sample to genotype
            String[] samples = hapmap.getSampleIDs();
            String[] genotypeStrings = hapmap.getGenotypes();

            Map<String, Genotype> genotypes = new HashMap<String, Genotype>(samples.length);
            for ( int i = 0; i < samples.length; i++ ) {
                // ignore bad genotypes
                if ( genotypeStrings[i].contains("N") )
                    continue;

                String a1 = genotypeStrings[i].substring(0,1);
                String a2 = genotypeStrings[i].substring(1);

                Allele allele1 = new Allele(a1, refAllele.basesMatch(a1));
                Allele allele2 = new Allele(a2, refAllele.basesMatch(a2));

                ArrayList<Allele> myAlleles = new ArrayList<Allele>(2);
                myAlleles.add(allele1);
                myAlleles.add(allele2);
                alleles.add(allele1);
                alleles.add(allele2);

                Genotype g = new Genotype(samples[i], myAlleles);
                genotypes.put(samples[i], g);
            }

            VariantContext vc = new VariantContext(name, hapmap.getLocation(), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, new HashMap<String, String>());
            vc.validate();
            return vc;
       }
    }
}

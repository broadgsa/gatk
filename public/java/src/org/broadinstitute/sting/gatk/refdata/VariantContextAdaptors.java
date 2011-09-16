package org.broadinstitute.sting.gatk.refdata;

import net.sf.samtools.util.SequenceUtil;
import org.broad.tribble.Feature;
import org.broad.tribble.annotation.Strand;
import org.broad.tribble.dbsnp.OldDbSNPFeature;
import org.broad.tribble.gelitext.GeliTextFeature;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.codecs.hapmap.RawHapMapFeature;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeader;
import org.broadinstitute.sting.utils.codecs.vcf.VCFHeaderLine;
import org.broadinstitute.sting.utils.variantcontext.*;

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

    private static Map<Class<? extends Feature>,VCAdaptor> adaptors = new HashMap<Class<? extends Feature>,VCAdaptor>();

    static {
        PluginManager<VCAdaptor> vcAdaptorManager = new PluginManager<VCAdaptor>(VCAdaptor.class);
        List<VCAdaptor> adaptorInstances = vcAdaptorManager.createAllTypes();
        for(VCAdaptor adaptor: adaptorInstances)
            adaptors.put(adaptor.getAdaptableFeatureType(),adaptor);
    }

    public static boolean canBeConvertedToVariantContext(Object variantContainingObject) {
        return adaptors.containsKey(variantContainingObject.getClass());
    }

    /** generic superclass */
    public interface VCAdaptor {
        /**
         * Gets the type of feature that this adaptor can 'adapt' into a VariantContext.
         * @return Type of adaptable feature.  Must be a Tribble feature class.
         */
        Class<? extends Feature> getAdaptableFeatureType();
        VariantContext convert(String name, Object input, ReferenceContext ref);
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
    private static class VariantContextAdaptor implements VCAdaptor {
        /**
         * 'Null' adaptor; adapts variant contexts to variant contexts.
         * @return VariantContext.
         */
        @Override
        public Class<? extends Feature> getAdaptableFeatureType() { return VariantContext.class; }

        // already a VC, just cast and return it
        @Override        
        public VariantContext convert(String name, Object input, ReferenceContext ref) {
            return (VariantContext)input;
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // dbSNP to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class DBSnpAdaptor implements VCAdaptor {
        private static boolean isSNP(OldDbSNPFeature feature) {
            return feature.getVariantType().contains("single") && feature.getLocationType().contains("exact");
        }

        private static boolean isMNP(OldDbSNPFeature feature) {
            return feature.getVariantType().contains("mnp") && feature.getLocationType().contains("range");
        }

        private static boolean isInsertion(OldDbSNPFeature feature) {
            return feature.getVariantType().contains("insertion");
        }

        private static boolean isDeletion(OldDbSNPFeature feature) {
            return feature.getVariantType().contains("deletion");
        }

        private static boolean isIndel(OldDbSNPFeature feature) {
            return isInsertion(feature) || isDeletion(feature) || isComplexIndel(feature);
        }

        public static boolean isComplexIndel(OldDbSNPFeature feature) {
            return feature.getVariantType().contains("in-del");
        }

        /**
         * gets the alternate alleles.  This method should return all the alleles present at the location,
         * NOT including the reference base.  This is returned as a string list with no guarantee ordering
         * of alleles (i.e. the first alternate allele is not always going to be the allele with the greatest
         * frequency).
         *
         * @return an alternate allele list
         */
        public static List<String> getAlternateAlleleList(OldDbSNPFeature feature) {
            List<String> ret = new ArrayList<String>();
            for (String allele : getAlleleList(feature))
                if (!allele.equals(String.valueOf(feature.getNCBIRefBase()))) ret.add(allele);
            return ret;
        }

        /**
         * gets the alleles.  This method should return all the alleles present at the location,
         * including the reference base.  The first allele should always be the reference allele, followed
         * by an unordered list of alternate alleles.
         *
         * @return an alternate allele list
         */
        public static List<String> getAlleleList(OldDbSNPFeature feature) {
            List<String> alleleList = new ArrayList<String>();
            // add ref first
            if ( feature.getStrand() == Strand.POSITIVE )
                alleleList = Arrays.asList(feature.getObserved());
            else
                for (String str : feature.getObserved())
                    alleleList.add(SequenceUtil.reverseComplement(str));
            if ( alleleList.size() > 0 && alleleList.contains(feature.getNCBIRefBase())
                    && !alleleList.get(0).equals(feature.getNCBIRefBase()) )
                Collections.swap(alleleList, alleleList.indexOf(feature.getNCBIRefBase()), 0);

            return alleleList;
        }

        /**
         * Converts non-VCF formatted dbSNP records to VariantContext. 
         * @return OldDbSNPFeature.
         */
        @Override
        public Class<? extends Feature> getAdaptableFeatureType() { return OldDbSNPFeature.class; }

        @Override        
        public VariantContext convert(String name, Object input, ReferenceContext ref) {
            OldDbSNPFeature dbsnp = (OldDbSNPFeature)input;
            if ( ! Allele.acceptableAlleleBases(dbsnp.getNCBIRefBase()) )
                return null;
            Allele refAllele = Allele.create(dbsnp.getNCBIRefBase(), true);

            if ( isSNP(dbsnp) || isIndel(dbsnp) || isMNP(dbsnp) || dbsnp.getVariantType().contains("mixed") ) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                alleles.add(refAllele);

                // add all of the alt alleles
                boolean sawNullAllele = refAllele.isNull();
                for ( String alt : getAlternateAlleleList(dbsnp) ) {
                    if ( ! Allele.acceptableAlleleBases(alt) ) {
                        //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                        return null;
                    }
                    Allele altAllele = Allele.create(alt, false);
                    alleles.add(altAllele);
                    if ( altAllele.isNull() )
                        sawNullAllele = true;
                }

                Map<String, Object> attributes = new HashMap<String, Object>();
                attributes.put(VariantContext.ID_KEY, dbsnp.getRsID());

                int index = dbsnp.getStart() - ref.getWindow().getStart() - 1;
                if ( index < 0 )
                    return null; // we weren't given enough reference context to create the VariantContext
                Byte refBaseForIndel = new Byte(ref.getBases()[index]);

                Map<String, Genotype> genotypes = null;
                VariantContext vc = new VariantContext(name, dbsnp.getChr(), dbsnp.getStart() - (sawNullAllele ? 1 : 0), dbsnp.getEnd() - (refAllele.isNull() ? 1 : 0), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, attributes, refBaseForIndel);
                return vc;
            } else
                return null; // can't handle anything else
        }
    }

    public static VCFHeader createVCFHeader(Set<VCFHeaderLine> hInfo, VariantContext vc) {
        HashSet<String> names = new LinkedHashSet<String>();
        for ( Genotype g : vc.getGenotypesSortedByName() ) {
            names.add(g.getSampleName());
        }

        return new VCFHeader(hInfo == null ? new HashSet<VCFHeaderLine>() : hInfo, names);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // GELI to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class GeliTextAdaptor implements VCAdaptor {
        /**
         * Converts Geli text records to VariantContext. 
         * @return GeliTextFeature.
         */
        @Override
        public Class<? extends Feature> getAdaptableFeatureType() { return GeliTextFeature.class; }

        /**
         * convert to a Variant Context, given:
         * @param name  the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @param ref   the reference context
         * @return a VariantContext object
         */
        @Override
        public VariantContext convert(String name, Object input, ReferenceContext ref) {
            GeliTextFeature geli = (GeliTextFeature)input;
            if ( ! Allele.acceptableAlleleBases(String.valueOf(geli.getRefBase())) )
                return null;
            Allele refAllele = Allele.create(String.valueOf(geli.getRefBase()), true);

            // make sure we can convert it
            if ( geli.getGenotype().isHet() || !geli.getGenotype().containsBase(geli.getRefBase())) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                List<Allele> genotypeAlleles = new ArrayList<Allele>();
                // add all of the alt alleles
                for ( char alt : geli.getGenotype().toString().toCharArray() ) {
                    if ( ! Allele.acceptableAlleleBases(String.valueOf(alt)) ) {
                        return null;
                    }
                    Allele allele = Allele.create(String.valueOf(alt), false);
                    if (!alleles.contains(allele) && !refAllele.basesMatch(allele.getBases())) alleles.add(allele);

                    // add the allele, first checking if it's reference or not
                    if (!refAllele.basesMatch(allele.getBases())) genotypeAlleles.add(allele);
                    else genotypeAlleles.add(refAllele);
                }

                Map<String, String> attributes = new HashMap<String, String>();
                Collection<Genotype> genotypes = new ArrayList<Genotype>();
                MutableGenotype call = new MutableGenotype(name, genotypeAlleles);

                // set the likelihoods, depth, and RMS mapping quality values
                //call.putAttribute(CalledGenotype.POSTERIORS_ATTRIBUTE_KEY,geli.getLikelihoods());
                //call.putAttribute(GeliTextWriter.MAXIMUM_MAPPING_QUALITY_ATTRIBUTE_KEY,geli.getMaximumMappingQual());
                //call.putAttribute(GeliTextWriter.READ_COUNT_ATTRIBUTE_KEY,geli.getDepthOfCoverage());

                // add the call to the genotype list, and then use this list to create a VariantContext
                genotypes.add(call);
                alleles.add(refAllele);
                VariantContext vc = VariantContextUtils.toVC(name, ref.getGenomeLocParser().createGenomeLoc(geli.getChr(),geli.getStart()), alleles, genotypes, geli.getLODBestToReference(), null, attributes);
                return vc;
            } else
                return null; // can't handle anything else
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // HapMap to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class HapMapAdaptor implements VCAdaptor {
        /**
         * Converts HapMap records to VariantContext. 
         * @return HapMapFeature.
         */
        @Override
        public Class<? extends Feature> getAdaptableFeatureType() { return RawHapMapFeature.class; }

        /**
         * convert to a Variant Context, given:
         * @param name  the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @param ref   the reference context
         * @return a VariantContext object
         */
        @Override        
        public VariantContext convert(String name, Object input, ReferenceContext ref) {
            if ( ref == null )
                throw new UnsupportedOperationException("Conversion from HapMap to VariantContext requires a reference context");

            RawHapMapFeature hapmap = (RawHapMapFeature)input;

            int index = hapmap.getStart() - ref.getWindow().getStart();
            if ( index < 0 )
                return null; // we weren't given enough reference context to create the VariantContext
            Byte refBaseForIndel = new Byte(ref.getBases()[index]);

            HashSet<Allele> alleles = new HashSet<Allele>();
            Allele refSNPAllele = Allele.create(ref.getBase(), true);
            int deletionLength = -1;

            Map<String, Allele> alleleMap = hapmap.getActualAlleles();
            // use the actual alleles, if available
            if ( alleleMap != null ) {
                alleles.addAll(alleleMap.values());
                Allele deletionAllele = alleleMap.get(RawHapMapFeature.INSERTION);  // yes, use insertion here (since we want the reference bases)
                if ( deletionAllele != null && deletionAllele.isReference() )
                    deletionLength = deletionAllele.length();
            } else {
                // add the reference allele for SNPs
                alleles.add(refSNPAllele);
            }

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
                ArrayList<Allele> myAlleles = new ArrayList<Allele>(2);

                // use the mapping to actual alleles, if available
                if ( alleleMap != null ) {
                    myAlleles.add(alleleMap.get(a1));
                    myAlleles.add(alleleMap.get(a2));
                } else {
                    // ignore indels (which we can't handle without knowing the alleles)
                    if ( genotypeStrings[i].contains("I") || genotypeStrings[i].contains("D") )
                        continue;

                    Allele allele1 = Allele.create(a1, refSNPAllele.basesMatch(a1));
                    Allele allele2 = Allele.create(a2, refSNPAllele.basesMatch(a2));

                    myAlleles.add(allele1);
                    myAlleles.add(allele2);
                    alleles.add(allele1);
                    alleles.add(allele2);
                }

                Genotype g = new Genotype(samples[i], myAlleles);
                genotypes.put(samples[i], g);
            }

            HashMap<String, Object> attrs = new HashMap<String, Object>(1);
            attrs.put(VariantContext.ID_KEY, hapmap.getName());

            long end = hapmap.getEnd();
            if ( deletionLength > 0 )
                end += deletionLength;
            VariantContext vc = new VariantContext(name, hapmap.getChr(), hapmap.getStart(), end, alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, attrs, refBaseForIndel);
            return vc;
       }
    }
}

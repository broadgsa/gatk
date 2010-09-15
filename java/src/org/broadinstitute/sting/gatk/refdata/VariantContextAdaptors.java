package org.broadinstitute.sting.gatk.refdata;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.gelitext.GeliTextFeature;
import org.broad.tribble.hapmap.HapMapFeature;
import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.MutableGenotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.utils.*;

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
        adaptors.put(DbSNPFeature.class, new DBSnpAdaptor());
        adaptors.put(HapMapFeature.class, new HapMapAdaptor());
        adaptors.put(GeliTextFeature.class, new GeliTextAdaptor());
        adaptors.put(VariantContext.class, new VariantContextAdaptor());
    }

    public static boolean canBeConvertedToVariantContext(Object variantContainingObject) {
        return adaptors.containsKey(variantContainingObject.getClass());
    }

    /** generic superclass */
    private static abstract class VCAdaptor {
        abstract VariantContext convert(String name, Object input, ReferenceContext ref);
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

    private static class VariantContextAdaptor extends VCAdaptor {
        // already a VC, just cast and return it
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            return (VariantContext)input;
        }
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // dbSNP to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class DBSnpAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            DbSNPFeature dbsnp = (DbSNPFeature)input;
            if ( ! Allele.acceptableAlleleBases(DbSNPHelper.getReference(dbsnp),true) )
                return null;
            Allele refAllele = Allele.create(DbSNPHelper.getReference(dbsnp), true);

            if ( DbSNPHelper.isSNP(dbsnp) || DbSNPHelper.isIndel(dbsnp) || dbsnp.getVariantType().contains("mixed") ) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                alleles.add(refAllele);

                // add all of the alt alleles
                for ( String alt : DbSNPHelper.getAlternateAlleleList(dbsnp) ) {
                    if ( ! Allele.acceptableAlleleBases(alt,false) ) {
                        //System.out.printf("Excluding dbsnp record %s%n", dbsnp);
                        return null;
                    }
                    alleles.add(Allele.create(alt, false));
                }

                Map<String, String> attributes = new HashMap<String, String>();
                attributes.put(VariantContext.ID_KEY, dbsnp.getRsID());
                Collection<Genotype> genotypes = null;
                VariantContext vc = new VariantContext(name, dbsnp.getChr(),dbsnp.getStart(),dbsnp.getEnd(), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, attributes);
                return vc;
            } else
                return null; // can't handle anything else
        }
    }



    private static Allele deletionAllele(ReferenceContext ref, int start, int len) {
        byte[] deletion = new byte[len];
        System.arraycopy(ref.getBases(), start, deletion, 0, len);
        return Allele.create(deletion, true);
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

    private static class GeliTextAdaptor extends VCAdaptor {
          /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @return a VariantContext object
         */
//        VariantContext convert(String name, Object input) {
//            return convert(name, input, null);
//        }

        /**
         * convert to a Variant Context, given:
         * @param name  the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @param ref   the reference context
         * @return a VariantContext object
         */
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            GeliTextFeature geli = (GeliTextFeature)input;
            if ( ! Allele.acceptableAlleleBases(String.valueOf(geli.getRefBase()),true) )
                return null;
            Allele refAllele = Allele.create(String.valueOf(geli.getRefBase()), true);

            // make sure we can convert it
            if ( geli.getGenotype().isHet() || !geli.getGenotype().containsBase(geli.getRefBase())) {
                // add the reference allele
                List<Allele> alleles = new ArrayList<Allele>();
                List<Allele> genotypeAlleles = new ArrayList<Allele>();
                // add all of the alt alleles
                for ( char alt : geli.getGenotype().toString().toCharArray() ) {
                    if ( ! Allele.acceptableAlleleBases(String.valueOf(alt),false) ) {
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
                VariantContext vc = VariantContextUtils.toVC(name, GenomeLocParser.createGenomeLoc(geli.getChr(),geli.getStart()), alleles, genotypes, geli.getLODBestToReference(), null, attributes);
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

    private static class HapMapAdaptor extends VCAdaptor {
          /**
         * convert to a Variant Context, given:
         * @param name the name of the ROD
         * @param input the Rod object, in this case a RodGeliText
         * @return a VariantContext object
         */
//        VariantContext convert(String name, Object input) {
//            return convert(name, input, null);
//        }

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

            HapMapFeature hapmap = (HapMapFeature)input;

            // add the reference allele
            HashSet<Allele> alleles = new HashSet<Allele>();
            Allele refAllele = Allele.create(ref.getBase(), true);
            alleles.add(refAllele);

            // make a mapping from sample to genotype
            String[] samples = hapmap.getSampleIDs();
            String[] genotypeStrings = hapmap.getGenotypes();

            Map<String, Genotype> genotypes = new HashMap<String, Genotype>(samples.length);
            for ( int i = 0; i < samples.length; i++ ) {
                // ignore bad genotypes
                if ( genotypeStrings[i].contains("N") || genotypeStrings[i].contains("I") || genotypeStrings[i].contains("D") )
                    continue;

                String a1 = genotypeStrings[i].substring(0,1);
                String a2 = genotypeStrings[i].substring(1);

                Allele allele1 = Allele.create(a1, refAllele.basesMatch(a1));
                Allele allele2 = Allele.create(a2, refAllele.basesMatch(a2));

                ArrayList<Allele> myAlleles = new ArrayList<Allele>(2);
                myAlleles.add(allele1);
                myAlleles.add(allele2);
                alleles.add(allele1);
                alleles.add(allele2);

                Genotype g = new Genotype(samples[i], myAlleles);
                genotypes.put(samples[i], g);
            }

            VariantContext vc = new VariantContext(name, hapmap.getChr(), hapmap.getStart(), hapmap.getEnd(), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, new HashMap<String, String>());
            return vc;
       }
    }
}

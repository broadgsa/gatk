package org.broadinstitute.sting.gatk.refdata;

import edu.mit.broad.picard.genotype.DiploidGenotype;
import edu.mit.broad.picard.genotype.geli.GenotypeLikelihoods;
import org.broad.tribble.dbsnp.DbSNPFeature;
import org.broad.tribble.gelitext.GeliTextFeature;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.MutableGenotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.utils.helpers.DbSNPHelper;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.CalledGenotype;
import org.broadinstitute.sting.utils.genotype.geli.GeliTextWriter;

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
        adaptors.put(VCFRecord.class, new VCFRecordAdaptor());
        adaptors.put(PlinkRod.class, new PlinkRodAdaptor());
        adaptors.put(HapMapROD.class, new HapMapAdaptor());
        adaptors.put(GeliTextFeature.class, new GeliTextAdaptor());
        adaptors.put(rodGELI.class, new GeliAdaptor());
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
                VariantContext vc = new VariantContext(name, GenomeLocParser.createGenomeLoc(dbsnp.getChr(),dbsnp.getStart(),dbsnp.getEnd()), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, attributes);
                return vc;
            } else
                return null; // can't handle anything else
        }
    }

    private static class VCFRecordAdaptor extends VCAdaptor {
        VariantContext convert(String name, Object input, ReferenceContext ref) {
            return vcfToVariantContext(name, (VCFRecord)input, ref);
        }
    }

    private static VariantContext vcfToVariantContext(String name, VCFRecord vcf, ReferenceContext ref) {
        if ( vcf.isReference() || vcf.isSNP() || vcf.isIndel() ) {
            // add the reference allele
            if ( ! Allele.acceptableAlleleBases(vcf.getReference(),true) ) {
                System.out.printf("Excluding vcf record %s%n", vcf);
                return null;
            }

            Set<String> filters = vcf.isFiltered() ? new HashSet<String>(Arrays.asList(vcf.getFilteringCodes())) : null;
            Map<String, String> attributes = new HashMap<String, String>(vcf.getInfoValues());
            attributes.put(VariantContext.ID_KEY, vcf.getID());

            // add all of the alt alleles
            List<Allele> alleles = new ArrayList<Allele>();
            Allele refAllele = determineRefAllele(vcf, ref);
            alleles.add(refAllele);

            for ( VCFGenotypeEncoding alt : vcf.getAlternateAlleles() ) {
                if ( ! Allele.acceptableAlleleBases(alt.getBases(),false) ) {
                    //System.out.printf("Excluding vcf record %s%n", vcf);
                    return null;
                }

                Allele allele;
                // special case: semi-deletion
                if ( vcf.isDeletion() && refAllele.length() > alt.getLength() ) {
                    byte[] semiDeletion = new byte[refAllele.length() - alt.getLength()];
                    System.arraycopy(ref.getBases(), alt.getLength(), semiDeletion, 0, refAllele.length() - alt.getLength());
                    allele = Allele.create(new String(semiDeletion), false);
                } else {
                    allele = Allele.create(alt.getBases(), false);
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
                    // todo -- fixme if we put GQ and FT into key itself
                    if ( ! e.getKey().equals(VCFConstants.GENOTYPE_QUALITY_KEY) && ! e.getKey().equals(VCFConstants.GENOTYPE_FILTER_KEY) )
                        fields.put(e.getKey(), e.getValue());
                }

                Set<String> genotypeFilters = new HashSet<String>();
                if ( vcfG.isFiltered() ) // setup the genotype filter fields
                    genotypeFilters.addAll(Arrays.asList(vcfG.getFields().get(VCFConstants.GENOTYPE_FILTER_KEY).split(";")));

                double qual = vcfG.isMissingQual() ? VariantContext.NO_NEG_LOG_10PERROR : vcfG.getNegLog10PError();
                Genotype g = new Genotype(vcfG.getSampleName(), genotypeAlleles, qual, genotypeFilters, fields, vcfG.getPhaseType() == VCFGenotypeRecord.PHASE.PHASED);
                genotypes.put(g.getSampleName(), g);
            }

            double qual = vcf.isMissingQual() ? VariantContext.NO_NEG_LOG_10PERROR : vcf.getNegLog10PError();

            GenomeLoc loc = GenomeLocParser.createGenomeLoc(vcf.getChr(),vcf.getStart());
            if ( vcf.isDeletion() )
                loc = GenomeLocParser.createGenomeLoc(loc.getContig(), loc.getStart(), loc.getStart()+refAllele.length()-1);

            VariantContext vc = new VariantContext(name, loc, alleles, genotypes, qual, filters, attributes);
            return vc;
        } else
            return null; // can't handle anything else
    }

    private static Allele determineRefAllele(VCFRecord vcf, ReferenceContext ref) {
        if ( ref == null )
            throw new StingException("Illegal determineRefAllele call!");

        Allele refAllele;
        if ( vcf.isInsertion() ) {
            refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
//        } else if ( ref == null ) {
//            refAllele = Allele.create(vcf.getReference(), true);
        } else if ( !vcf.isIndel() ) {
            refAllele = Allele.create(ref.getBase(), true);
            if ( (char)ref.getBase() != vcf.getReference().charAt(0) )
                throw new StingException("The VCF reference base (" + vcf.getReference().charAt(0) + ") doesn't match the actual reference base (" + (char)ref.getBase() + "); please check that you are using the appropriate reference file");
        } else if ( vcf.isDeletion() ) {
            int start = vcf.getPosition() - (int)ref.getWindow().getStart() + 1;
            int delLength = 0;
            for ( VCFGenotypeEncoding enc : vcf.getAlternateAlleles() ) {
                if ( enc.getLength() > delLength )
                    delLength = enc.getLength();
            }
            if ( delLength > ref.getWindow().getStop() - vcf.getPosition() )
                throw new IllegalArgumentException("Length of deletion is larger than reference context provided at " + ref.getLocus());

            refAllele = deletionAllele(ref, start, delLength);
        } else {
            throw new UnsupportedOperationException("Conversion of VCF type " + vcf.getType() + " is not supported.");
        }

        return refAllele;
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
    // Plink to VariantContext
    //
    // --------------------------------------------------------------------------------------------------------------

    private static class PlinkRodAdaptor extends VCAdaptor {
//        VariantContext convert(String name, Object input) {
//            return convert(name, input, null);
//        }

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
                                VCAllele = Allele.create(myPlinkAllele.getBases(), refAllele.equals(myPlinkAllele, true));
                            } else if ( myPlinkAllele.isNull() ) {
                                VCAllele = Allele.create(Allele.NULL_ALLELE_STRING, plink.isInsertion());
                            } else {
                                VCAllele = Allele.create(myPlinkAllele.getBases(), !plink.isInsertion());
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
                return vc;
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException(e.getMessage() + "; please make sure that e.g. a sample isn't present more than one time in your ped file");
            }
        }

        private Allele determineRefAllele(PlinkRod plink, ReferenceContext ref) {
            Allele refAllele;
            if ( !plink.isIndel() ) {
                refAllele = Allele.create(ref.getBase(), true);
            } else if ( plink.isInsertion() ) {
                refAllele = Allele.create(Allele.NULL_ALLELE_STRING, true);
            } else {
                long maxLength = ref.getWindow().getStop() - ref.getLocus().getStop();                
                if ( plink.getLength() > maxLength )
                    throw new UnsupportedOperationException("Plink conversion currently can only handle indels up to length " + maxLength);
                refAllele = deletionAllele(ref, 1, plink.getLength());

            }
            return refAllele;
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
                call.putAttribute(CalledGenotype.POSTERIORS_ATTRIBUTE_KEY,geli.getLikelihoods());
                call.putAttribute(GeliTextWriter.MAXIMUM_MAPPING_QUALITY_ATTRIBUTE_KEY,geli.getMaximumMappingQual());
                call.putAttribute(GeliTextWriter.READ_COUNT_ATTRIBUTE_KEY,geli.getDepthOfCoverage());

                // add the call to the genotype list, and then use this list to create a VariantContext
                genotypes.add(call);
                alleles.add(refAllele);
                VariantContext vc = new VariantContext(name, GenomeLocParser.createGenomeLoc(geli.getChr(),geli.getStart()), alleles, genotypes, geli.getLODBestToReference(), null, attributes);
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
            GenotypeLikelihoods geli = ((rodGELI)input).getGeliRecord();
            if ( ! Allele.acceptableAlleleBases(String.valueOf((char)geli.getReferenceBase()),true) )
                return null;
            Allele refAllele = Allele.create(String.valueOf((char)geli.getReferenceBase()), true);

            // add the reference allele
            List<Allele> alleles = new ArrayList<Allele>();
            alleles.add(refAllele);

            // add the two alt alleles
            if (!Allele.acceptableAlleleBases(String.valueOf((char) geli.getBestGenotype().getAllele1()),false)) {
                return null;
            }
            Allele allele1 = Allele.create(String.valueOf((char) geli.getBestGenotype().getAllele1()), false);
            if (!alleles.contains(allele1) && !allele1.equals(refAllele, true)) alleles.add(allele1);

            // add the two alt alleles
            if (!Allele.acceptableAlleleBases(String.valueOf((char) geli.getBestGenotype().getAllele2()),false)) {
                return null;
            }
            Allele allele2 = Allele.create(String.valueOf((char) geli.getBestGenotype().getAllele2()), false);
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

            HapMapROD hapmap = (HapMapROD)input;

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

            VariantContext vc = new VariantContext(name, hapmap.getLocation(), alleles, genotypes, VariantContext.NO_NEG_LOG_10PERROR, null, new HashMap<String, String>());
            return vc;
       }
    }
}

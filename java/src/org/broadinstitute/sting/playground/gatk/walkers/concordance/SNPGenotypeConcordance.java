package org.broadinstitute.sting.playground.gatk.walkers.concordance;

import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.utils.StingException;

import java.util.HashMap;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * Split up two call sets into their various concordance sets
 */
public class SNPGenotypeConcordance implements ConcordanceType {

    private double LOD = 5.0;

    private PrintWriter sameVarWriter = null;
    private PrintWriter oneVarWriter = null;
    private PrintWriter sameAlleleWriter = null;
    private PrintWriter diffAlleleWriter = null;
    private PrintWriter refVsVar1Writer = null;
    private PrintWriter refVsVar2Writer = null;
    private PrintWriter comboVarWriter = null;
    private PrintWriter coverageVar1Writer = null;
    private PrintWriter coverageVar2Writer = null;
    private String prefix;

    public SNPGenotypeConcordance(String prefix) {
        this.prefix = prefix;
    }

    public void initialize(HashMap<String,String> args) {
        try {
            sameVarWriter = new PrintWriter(prefix + ".snp_concordance.sameConfidentVariant.calls");
            oneVarWriter = new PrintWriter(prefix + ".snp_concordance.bothVariantOnlyOneIsConfident.calls");
            sameAlleleWriter = new PrintWriter(prefix + ".snp_concordance.sameVariantAlleleDifferentGenotype.calls");
            diffAlleleWriter = new PrintWriter(prefix + ".snp_concordance.differentVariantAllele.calls");
            refVsVar1Writer = new PrintWriter(prefix + ".snp_concordance.set1VariantSet2Ref.calls");
            refVsVar2Writer = new PrintWriter(prefix + ".snp_concordance.set1RefSet2Variant.calls");
            comboVarWriter = new PrintWriter(prefix + ".snp_concordance.confidentVariantWhenCombined.calls");
            coverageVar1Writer = new PrintWriter(prefix + ".snp_concordance.set1ConfidentVariantSet2NoCoverage.calls");
            coverageVar2Writer = new PrintWriter(prefix + ".snp_concordance.set1NoCoverageSet2ConfidentVariant.calls");
        } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }

        if ( args.get("lod") != null )
            LOD = Double.valueOf(args.get("lod"));
    }

    public void computeConcordance(RefMetaDataTracker tracker, ReferenceContext ref) {
        AllelicVariant call1 = (AllelicVariant)tracker.lookup("callset1", null);
        AllelicVariant call2 = (AllelicVariant)tracker.lookup("callset2", null);

        // the only reason they would be null is a lack of coverage
        if ( call1 == null || call2 == null ) {
            if ( call1 != null && call1.isSNP() && call1.getVariationConfidence() >= LOD )
                printVariant(coverageVar1Writer, call1);
            else if ( call2 != null && call2.isSNP() && call2.getVariationConfidence() >= LOD )
                printVariant(coverageVar2Writer, call2);
            return;
        }

        double bestVsRef1 = call1.getVariationConfidence();
        double bestVsRef2 = call2.getVariationConfidence();
        //double bestVsNext1 = call1.getConsensusConfidence();
        //double bestVsNext2 = call2.getConsensusConfidence();
        String genotype1 = call1.getGenotype().get(0);
        String genotype2 = call2.getGenotype().get(0);

        // are they both variant SNPs?
        if ( call1.isSNP() && call2.isSNP() ) {

            // are they both confident calls?
            if ( bestVsRef1 >= LOD && bestVsRef2 >= LOD ) {
                // same genotype
                if ( genotype1.equals(genotype2) )
                    printVariant(sameVarWriter, call1);

                // same allele, different genotype
                else if ( sameVariantAllele(genotype1, genotype2, ref.getBase()) )
                    printVariant(sameAlleleWriter, call1);

                // different variant allele
                else
                    printVariant(diffAlleleWriter, call1);
            }

            // confident only when combined
            else if ( bestVsRef1 < LOD && bestVsRef2 < LOD && bestVsRef1 + bestVsRef2 >= LOD ) {
                printVariant(comboVarWriter, call1);
            }

            // only one is confident variant
            else if ( (bestVsRef1 < LOD && bestVsRef2 >= LOD) || (bestVsRef1 >= LOD && bestVsRef2 < LOD) ) {
                printVariant(oneVarWriter, call1);                
            }
        }

        // one is variant and the other is ref
        else if ( call1.isSNP() && call2.isReference() && bestVsRef1 >= LOD )
             printVariant(refVsVar1Writer, call1);
        else if ( call2.isSNP() && call1.isReference() && bestVsRef2 >= LOD )
            printVariant(refVsVar2Writer, call1);
    }

    private boolean sameVariantAllele(String genotype1, String genotype2, char ref) {
        if ( genotype1.length() < 2 || genotype2.length() < 2 )
            return genotype1.equals(genotype2);
        char altAllele1 = genotype1.charAt(0) != ref ? genotype1.charAt(0) : genotype1.charAt(1);
        char altAllele2 = genotype2.charAt(0) != ref ? genotype2.charAt(0) : genotype2.charAt(1);
        return altAllele1 == altAllele2;
    }

    private static void printVariant(PrintWriter writer, AllelicVariant variant) {
        writer.println(variant.toString());
    }

    public void cleanup() {
        sameVarWriter.close();
        oneVarWriter.close();
        sameAlleleWriter.close();
        diffAlleleWriter.close();
        refVsVar1Writer.close();
        refVsVar2Writer.close();
        comboVarWriter.close();
        coverageVar1Writer.close();
        coverageVar2Writer.close();
    }
}
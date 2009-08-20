package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * Split up two call sets into their various concordance sets
 */
@Requires(value={DataSource.REFERENCE},referenceMetaData={@RMD(name="callset1",type=AllelicVariant.class),@RMD(name="callset2",type=AllelicVariant.class)})
public class CallSetsConcordanceWalker extends RefWalker<Integer, Integer> {
    @Argument(fullName="outputPrefix", shortName="prefix", doc="File prefix for output; turns on ALL output options below", required=false)
    private String PREFIX = null;
    @Argument(fullName="bothConfidentSameVar", shortName="confSameVar", doc="File to which calls that are confidently called var (identically) in both sets should be written", required=false)
    private String CONF_SAME_VAR = null;
    @Argument(fullName="oneConfidentVar", shortName="oneConfVar", doc="File to which calls that are confidently called var in only one of the sets should be written", required=false)
    private String CONF_ONE_VAR = null;
    @Argument(fullName="bothConfidentSameAllele", shortName="confSameAllele", doc="File to which calls that are confidently called var, with different genotypes, but with the same alt allele should be written", required=false)
    private String CONF_SAME_ALLELE = null;
    @Argument(fullName="bothConfidentDifferentAllele", shortName="confDiffAllele", doc="File to which calls that are confidently called var, with different genotypes, and with different alt alleles should be written", required=false)
    private String CONF_DIFF_ALLELE = null;
    @Argument(fullName="confidentRefVsConfidentVar", shortName="confRefVar", doc="File to which calls that are confidently called var in one set and ref in the other should be written", required=false)
    private String CONF_REF_VAR = null;
    @Argument(fullName="confidentVarWhenCombined", shortName="confComboVar", doc="File to which calls that are not confidently called var in each set separately but are when combined should be written", required=false)
    private String CONF_COMBO_VAR = null;
    @Argument(fullName="noCoverageVar", shortName="coverageVar", doc="File to which calls that are confidently called var in one set and for which there is no coverage in the other set should be written", required=false)
    private String NO_COVERAGE_VAR = null;
    @Argument(fullName="LOD_threshold", shortName="LOD", doc="LOD score for determining confidence", required=false)
    private Double LOD = 5.0;

    private PrintWriter sameVarWriter = null;
    private PrintWriter oneVarWriter = null;
    private PrintWriter sameAlleleWriter = null;
    private PrintWriter diffAlleleWriter = null;
    private PrintWriter refVsVarWriter = null;
    private PrintWriter comboVarWriter = null;
    private PrintWriter coverageVarWriter = null;

    /**
     * Prepare the output file and the list of available features.
     */
    public void initialize() {
        try {
            if ( CONF_SAME_VAR != null )
                sameVarWriter = new PrintWriter(CONF_SAME_VAR);
            else if ( PREFIX != null )
                sameVarWriter = new PrintWriter(PREFIX + ".sameConfidentVariant.calls");
            if ( CONF_ONE_VAR != null )
                oneVarWriter = new PrintWriter(CONF_ONE_VAR);
            else if ( PREFIX != null )
                oneVarWriter = new PrintWriter(PREFIX + ".bothVariantOnlyOneIsConfident.calls");
            if ( CONF_SAME_ALLELE != null )
                sameAlleleWriter = new PrintWriter(CONF_SAME_ALLELE);
            else if ( PREFIX != null )
                sameAlleleWriter = new PrintWriter(PREFIX + ".sameVariantAlleleDifferentGenotype.calls");
            if ( CONF_DIFF_ALLELE != null )
                diffAlleleWriter = new PrintWriter(CONF_DIFF_ALLELE);
            else if ( PREFIX != null )
                diffAlleleWriter = new PrintWriter(PREFIX + ".differentVariantAllele.calls");
            if ( CONF_REF_VAR != null )
                refVsVarWriter = new PrintWriter(CONF_REF_VAR);
            else if ( PREFIX != null )
                refVsVarWriter = new PrintWriter(PREFIX + ".oneRefOneVariant.calls");
            if ( CONF_COMBO_VAR != null )
                comboVarWriter = new PrintWriter(CONF_COMBO_VAR);
            else if ( PREFIX != null )
                comboVarWriter = new PrintWriter(PREFIX + ".confidentVariantWhenCombined.calls");
            if ( NO_COVERAGE_VAR != null )
                coverageVarWriter = new PrintWriter(NO_COVERAGE_VAR);
            else if ( PREFIX != null )
                coverageVarWriter = new PrintWriter(PREFIX + ".oneConfidentVariantOneNoCoverage.calls");
         } catch (FileNotFoundException e) {
            throw new StingException(String.format("Could not open file(s) for writing"));
        }
    }

    /**
     * Initialize the number of loci processed to zero.
     *
     * @return 0
     */
    public Integer reduceInit() { return 0; }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        AllelicVariant call1 = (AllelicVariant)tracker.lookup("callset1", null);
        AllelicVariant call2 = (AllelicVariant)tracker.lookup("callset2", null);

        if ( call1 == null || call2 == null ) {
            if ( call1 != null && call1.isSNP() && call1.getVariationConfidence() >= LOD )
                printVariant(coverageVarWriter, call1);
            else if ( call2 != null && call2.isSNP() && call2.getVariationConfidence() >= LOD )
                printVariant(coverageVarWriter, call2);
            return 1;
        }

        double bestVsRef1 = call1.getVariationConfidence();
        double bestVsRef2 = call2.getVariationConfidence();
        //double bestVsNext1 = call1.getConsensusConfidence();
        //double bestVsNext2 = call2.getConsensusConfidence();
        String genotype1 = call1.getGenotype().get(0);
        String genotype2 = call2.getGenotype().get(0);

        // are they both variants?
        if ( call1.isSNP() && call2.isSNP() ) {
            // are they both confident calls?
            if ( bestVsRef1 >= LOD && bestVsRef2 >= LOD ) {
                if ( genotype1.equals(genotype2) )
                    printVariant(sameVarWriter, call1);
                else if ( sameVariantAllele(genotype1, genotype2, ref.getBase()) )
                    printVariant(sameAlleleWriter, call1);
                else
                    printVariant(diffAlleleWriter, call1);
            } else if ( bestVsRef1 < LOD && bestVsRef2 < LOD && bestVsRef1 + bestVsRef2 >= LOD ) {
                printVariant(comboVarWriter, call1);
            } else if ( (bestVsRef1 < LOD && bestVsRef2 >= LOD) || (bestVsRef1 >= LOD && bestVsRef2 < LOD) ) {
                printVariant(oneVarWriter, call1);                
            }
        } else if ( (call1.isSNP() && call2.isReference() && bestVsRef1 >= LOD) ||
                    (call2.isSNP() && call1.isReference() && bestVsRef2 >= LOD) ) {
            printVariant(refVsVarWriter, call1);
        }

        return 1;
    }

    private boolean sameVariantAllele(String genotype1, String genotype2, char ref) {
        if ( genotype1.length() < 2 || genotype2.length() < 2 )
            return genotype1.equals(genotype2);
        char altAllele1 = genotype1.charAt(0) != ref ? genotype1.charAt(0) : genotype1.charAt(1);
        char altAllele2 = genotype2.charAt(0) != ref ? genotype2.charAt(0) : genotype2.charAt(1);
        return altAllele1 == altAllele2;
    }

    static void printVariant(PrintWriter writer, AllelicVariant variant) {
        if ( writer != null && variant != null )
            writer.println(variant.toString());
    }

    /**
     * Increment the number of loci processed.
     *
     * @param value result of the map.
     * @param sum   accumulator for the reduce.
     * @return the new number of loci processed.
     */
    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        out.printf("Processed %d loci.\n", result);

        if ( sameVarWriter != null )
            sameVarWriter.close();
        if ( oneVarWriter != null )
            oneVarWriter.close();
        if ( sameAlleleWriter != null )
            sameAlleleWriter.close();
        if ( diffAlleleWriter != null )
            diffAlleleWriter.close();
        if ( refVsVarWriter != null )
            refVsVarWriter.close();
        if ( comboVarWriter != null )
            comboVarWriter.close();
        if ( coverageVarWriter != null )
            coverageVarWriter.close();
    }
}
/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFWriter;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;

/**
 * Takes a VCF file, selects variants based on sample(s) in which it was found and/or on various annotation criteria,
 * recompute the value of certain annotations based on the new sample set, and output a new VCF with the results.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
public class SelectVariants extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="sample", shortName="sn", doc="Sample(s) to include.  Can be a single sample, specified multiple times for many samples, a file containing sample names, a regular expression to select many samples, or any combination thereof.", required=false)
    public Set<String> SAMPLE_EXPRESSIONS;

    @Argument(shortName="select", doc="One or more criteria to use when selecting the data.  Evaluated *after* the specified samples are extracted and the INFO-field annotations are updated.", required=false)
    public ArrayList<String> SELECT_EXPRESSIONS = new ArrayList<String>();

    @Argument(fullName="excludeNonVariants", shortName="env", doc="Don't include loci found to be non-variant after the subsetting procedure.", required=false)
    private boolean EXCLUDE_NON_VARIANTS = false;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered loci.", required=false)
    private boolean EXCLUDE_FILTERED = false;

    @Argument(fullName="discordance", shortName =  "disc", doc="Output variants that were not called on a ROD comparison track. Use -disc ROD_NAME", required=false)
    private String discordanceRodName = "";

    @Argument(fullName="concordance", shortName =  "conc", doc="Output variants that were also called on a ROD comparison track. Use -conc ROD_NAME", required=false)
    private String concordanceRodName = "";

    @Argument(fullName="family_structure", shortName="family", doc="USE YAML FILE INSTEAD (-SM) !!! string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
    private String FAMILY_STRUCTURE = "";

    @Argument(fullName="mendelianViolation", shortName="mv", doc="output mendelian violation sites only. Sample metadata information will be taken from YAML file (passed with -SM)", required=false)
    private Boolean MENDELIAN_VIOLATIONS = false;

    @Argument(fullName="mendelianViolationQualThreshold", shortName="mvq", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    private double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 0;

    @Argument(fullName="select_random_number", shortName="number", doc="Selects a number of variants at random from the variant track. Variants are kept in memory to guarantee that n variants will be output, so use it only for a reasonable number of variants. Use select_random_fraction for larger numbers of variants", required=false)
    private int numRandom = 0;

    @Argument(fullName="select_random_fraction", shortName="fraction", doc="Selects a fraction (a number between 0 and 1) of the total variants at random from the variant track. Routine is based on probability, so the final result is not guaranteed to carry the exact fraction. Can be used for large fractions", required=false)
    private double fractionRandom = 0;


    /* Private class used to store the intermediate variants in the integer random selection process */
    private class RandomVariantStructure {
        private VariantContext vc;
        private byte refBase;

        RandomVariantStructure(VariantContext vcP, byte refBaseP) {
            vc = vcP;
            refBase = refBaseP;
        }

        public void set (VariantContext vcP, byte refBaseP) {
            vc = vcP;
            refBase = refBaseP;
        }

    }


    private ArrayList<String> selectNames = new ArrayList<String>();
    private List<VariantContextUtils.JexlVCMatchExp> jexls = null;

    private Set<String> samples = new HashSet<String>();

    private boolean DISCORDANCE_ONLY = false;
    private boolean CONCORDANCE_ONLY = false;

    private MendelianViolation mv;

    /* default name for the variant dataset (VCF) */
    private final String variantRodName = "variant";


    /* variables used by the SELECT RANDOM modules */
    private boolean SELECT_RANDOM_NUMBER = false;
    private boolean SELECT_RANDOM_FRACTION = false;
    private int variantNumber = 0;
    private int nVariantsAdded = 0;
    private int positionToAdd = 0;
    private RandomVariantStructure [] variantArray;





    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        // Get list of samples to include in the output
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantRodName);

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        samples = SampleUtils.getSamplesFromCommandLineInput(vcfSamples, SAMPLE_EXPRESSIONS);
        for (String sample : samples) {
            logger.info("Including sample '" + sample + "'");
        }

        // Initialize VCF header
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        headerLines.add(new VCFHeaderLine("source", "SelectVariants"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));

        for (int i = 0; i < SELECT_EXPRESSIONS.size(); i++) {
            // It's not necessary that the user supply select names for the JEXL expressions, since those
            // expressions will only be needed for omitting records.  Make up the select names here.
            selectNames.add(String.format("select-%d", i));
        }

        jexls = VariantContextUtils.initializeMatchExps(selectNames, SELECT_EXPRESSIONS);

        // Look at the parameters to decide which analysis to perform
        DISCORDANCE_ONLY = discordanceRodName.length() > 0;
        if (DISCORDANCE_ONLY) logger.info("Selecting only variants discordant with the track: " + discordanceRodName);

        CONCORDANCE_ONLY = concordanceRodName.length() > 0;
        if (CONCORDANCE_ONLY) logger.info("Selecting only variants concordant with the track: " + concordanceRodName);

        if (MENDELIAN_VIOLATIONS)
            mv = new MendelianViolation(getToolkit(), MENDELIAN_VIOLATION_QUAL_THRESHOLD);
        else if (!FAMILY_STRUCTURE.isEmpty()) {
            mv = new MendelianViolation(FAMILY_STRUCTURE, MENDELIAN_VIOLATION_QUAL_THRESHOLD);
            MENDELIAN_VIOLATIONS = true;
        }

        SELECT_RANDOM_NUMBER = numRandom > 0;
        if (SELECT_RANDOM_NUMBER) {
            logger.info("Selecting " + numRandom + " variants at random from the variant track");
            variantArray = new RandomVariantStructure[numRandom];
        }

        SELECT_RANDOM_FRACTION = fractionRandom > 0;
        if (SELECT_RANDOM_FRACTION) logger.info("Selecting " + fractionRandom + " variants at random from the variant track");
    }

    /**
     * Subset VC record if necessary and emit the modified record (provided it satisfies criteria for printing)
     *
     * @param  tracker   the ROD tracker
     * @param  ref       reference information
     * @param  context   alignment info
     * @return 1 if the record was printed to the output file, 0 if otherwise
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> vcs = tracker.getVariantContexts(ref, variantRodName, null, context.getLocation(), true, false);

        // If a discordance rod name was provided, we select the variants present in the discordance that
        // are not present in the variant.
        if (DISCORDANCE_ONLY) {
            if (vcs == null || vcs.size() == 0) {
                VariantContext compVC = tracker.getVariantContext(ref, discordanceRodName, context.getLocation());
                if (compVC != null) {
                    // if the genotype of any of the samples we're looking it is HET or HomVar, we write this variant
                    for (String sample : samples) {
                        Genotype g = compVC.getGenotype(sample);
                        if ( g != null && (g.isHet() || g.isHomVar()) ) {
                            vcfWriter.add(compVC, ref.getBase());
                            return 1;
                        }
                    }
                }
            }
        }

        if ( vcs == null || vcs.size() == 0) {
            return 0;
        }

        for (VariantContext vc : vcs) {
            if (MENDELIAN_VIOLATIONS) {
                if (mv.isViolation(vc)) {
                    vcfWriter.add(vc, ref.getBase());
                }
            }

            else if (CONCORDANCE_ONLY) {
                VariantContext compVC = tracker.getVariantContext(ref, concordanceRodName, context.getLocation());
                if (compVC != null) {
                    vcfWriter.add(vc, ref.getBase());
                    break; // we only want one line per event in this case.
                }
            }

            else {
                VariantContext sub = subsetRecord(vc, samples);
                if ( (sub.isPolymorphic() || !EXCLUDE_NON_VARIANTS) && (!sub.isFiltered() || !EXCLUDE_FILTERED) ) {
                        //System.out.printf("%s%n",sub.toString());
                    for ( VariantContextUtils.JexlVCMatchExp jexl : jexls ) {
                        if ( !VariantContextUtils.match(sub, jexl) ) {
                            return 0;
                        }
                    }
                    if (SELECT_RANDOM_NUMBER) {
                        randomlyAddVariant(++variantNumber, sub, ref.getBase());
                    }
                    else if (!SELECT_RANDOM_FRACTION || GenomeAnalysisEngine.getRandomGenerator().nextDouble() < fractionRandom) {
                        vcfWriter.add(sub, ref.getBase());
                    }
                }
            }
        }

        return 1;
    }


    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    public void onTraversalDone(Integer result) {
        logger.info(result + " records processed.");

        if (SELECT_RANDOM_NUMBER) {
            int positionToPrint = positionToAdd;
            for (int i=0; i<numRandom; i++) {
                vcfWriter.add(variantArray[positionToPrint].vc, variantArray[positionToPrint].refBase);
                positionToPrint = nextCircularPosition(positionToPrint, numRandom);
            }
        }
    }



    /**
     * Helper method to subset a VC record, modifying some metadata stored in the INFO field (i.e. AN, AC, AF).
     *
     * @param vc       the VariantContext record to subset
     * @param samples  the samples to extract
     * @return the subsetted VariantContext
     */                                                             
    private VariantContext subsetRecord(VariantContext vc, Set<String> samples) {
        if ( samples == null || samples.isEmpty() )
            return vc;

        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
        for ( Map.Entry<String, Genotype> genotypePair : vc.getGenotypes().entrySet() ) {
            if ( samples.contains(genotypePair.getKey()) )
                genotypes.add(genotypePair.getValue());
        }
        
        VariantContext sub = vc.subContextFromGenotypes(genotypes);

        HashMap<String, Object> attributes = new HashMap<String, Object>(sub.getAttributes());

        int depth = 0;
        for (String sample : sub.getSampleNames()) {
            Genotype g = sub.getGenotype(sample);

            if (g.isNotFiltered() && g.isCalled()) {
                
                String dp = (String) g.getAttribute("DP");
                if (dp != null && ! dp.equals(VCFConstants.MISSING_DEPTH_v3) && ! dp.equals(VCFConstants.MISSING_VALUE_v4) ) {
                    depth += Integer.valueOf(dp);
                }
            }
        }


        VariantContextUtils.calculateChromosomeCounts(sub,attributes,false);
        attributes.put("DP", depth);

        sub = VariantContext.modifyAttributes(sub, attributes);

        return sub;
    }

    private void randomlyAddVariant(int rank, VariantContext vc, byte refBase) {
        if (nVariantsAdded < numRandom)
            variantArray[nVariantsAdded++] = new RandomVariantStructure(vc, refBase);

        else {
            double v = GenomeAnalysisEngine.getRandomGenerator().nextDouble();
            double t = (1.0/(rank-numRandom+1));
            if ( v < t) {
                variantArray[positionToAdd].set(vc, refBase);
                nVariantsAdded++;
                positionToAdd = nextCircularPosition(positionToAdd, numRandom);
            }
        }
    }

    private int nextCircularPosition(int cur, int size) {
        if ((cur + 1) == size)
            return 0;
        return cur + 1;
    }
}

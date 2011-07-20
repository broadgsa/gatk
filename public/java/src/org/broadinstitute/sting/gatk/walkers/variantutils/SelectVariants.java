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

import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.variantcontext.*;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.annotation.AnnotationFormatError;
import java.util.*;

/**
 * Takes a VCF file, selects variants based on sample(s) in which it was found and/or on various annotation criteria,
 * recompute the value of certain annotations based on the new sample set, and output a new VCF with the results.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
public class SelectVariants extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="sample_name", shortName="sn", doc="Sample name to be included in the analysis. Can be specified multiple times.", required=false)
    public Set<String> sampleNames;

    @Argument(fullName="sample_expressions", shortName="se", doc="Regular expression to select many samples from the ROD tracks provided. Can be specified multiple times.", required=false)
    public Set<String> sampleExpressions;

    @Argument(fullName="sample_file", shortName="sf", doc="File containing a list of samples (one per line). Can be specified multiple times", required=false)
    public Set<File> sampleFiles;

    @Argument(shortName="select", doc="One or more criteria to use when selecting the data.  Evaluated *after* the specified samples are extracted and the INFO-field annotations are updated.", required=false)
    public ArrayList<String> SELECT_EXPRESSIONS = new ArrayList<String>();

    @Argument(fullName="excludeNonVariants", shortName="env", doc="Don't include loci found to be non-variant after the subsetting procedure.", required=false)
    private boolean EXCLUDE_NON_VARIANTS = false;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered loci in the analysis.", required=false)
    private boolean EXCLUDE_FILTERED = false;

    @Argument(fullName="keepOriginalAC", shortName="keepOriginalAC", doc="Don't include filtered loci.", required=false)
    private boolean KEEP_ORIGINAL_CHR_COUNTS = false;

    @Argument(fullName="discordance", shortName =  "disc", doc="Output variants that were not called on a ROD comparison track. Use -disc ROD_NAME", required=false)
    private String discordanceRodName = "";

    @Argument(fullName="concordance", shortName =  "conc", doc="Output variants that were also called on a ROD comparison track. Use -conc ROD_NAME", required=false)
    private String concordanceRodName = "";

    @Hidden
    @Argument(fullName="inputAF", shortName =  "inputAF", doc="", required=false)
    private String inputAFRodName = "";

    @Hidden
    @Argument(fullName="keepAFSpectrum", shortName="keepAF", doc="Don't include loci found to be non-variant after the subsetting procedure.", required=false)
    private boolean KEEP_AF_SPECTRUM = false;

    @Hidden
    @Argument(fullName="afFile", shortName="afFile", doc="The output recal file used by ApplyRecalibration", required=false)
    private File AF_FILE = new File("");

    @Hidden
    @Argument(fullName="family_structure_file", shortName="familyFile", doc="USE YAML FILE INSTEAD (-SM) !!! string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
    private File FAMILY_STRUCTURE_FILE = null;

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

    @Argument(fullName="selectSNPs", shortName="snps", doc="Select only SNPs.", required=false)
    private boolean SELECT_SNPS = false;

    @Argument(fullName="selectIndels", shortName="indels", doc="Select only Indels.", required=false)
    private boolean SELECT_INDELS = false;

    @Hidden
     @Argument(fullName="outMVFile", shortName="outMVFile", doc="USE YAML FILE INSTEAD (-SM) !!! string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
      private String outMVFile = null;

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

    private TreeSet<String> samples = new TreeSet<String>();
    private boolean NO_SAMPLES_SPECIFIED = false;

    private boolean DISCORDANCE_ONLY = false;
    private boolean CONCORDANCE_ONLY = false;

    private Set<MendelianViolation> mvSet = new HashSet<MendelianViolation>();

    /* default name for the variant dataset (VCF) */
    private final String variantRodName = "variant";


    /* variables used by the SELECT RANDOM modules */
    private boolean SELECT_RANDOM_NUMBER = false;
    private boolean SELECT_RANDOM_FRACTION = false;
    private int variantNumber = 0;
    private int nVariantsAdded = 0;
    private int positionToAdd = 0;
    private RandomVariantStructure [] variantArray;


    /* Variables used for random selection with AF boosting */
    private ArrayList<Double> afBreakpoints = null;
    private ArrayList<Double> afBoosts = null;
    double bkDelta = 0.0;


        private PrintStream outMVFileStream = null;


    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        // Get list of samples to include in the output
        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add(variantRodName);

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

        Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
        Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

        samples.addAll(samplesFromFile);
        samples.addAll(samplesFromExpressions);
        if (sampleNames != null)
            samples.addAll(sampleNames);

        if(samples.isEmpty()) {
            samples.addAll(vcfSamples);
            NO_SAMPLES_SPECIFIED = true;
        }

        for (String sample : samples) {
            logger.info("Including sample '" + sample + "'");
        }

        // Initialize VCF header
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        headerLines.add(new VCFHeaderLine("source", "SelectVariants"));

        if (KEEP_ORIGINAL_CHR_COUNTS) {
            headerLines.add(new VCFFormatHeaderLine("AC_Orig", 1, VCFHeaderLineType.Integer, "Original AC"));
            headerLines.add(new VCFFormatHeaderLine("AF_Orig", 1, VCFHeaderLineType.Float, "Original AF"));
            headerLines.add(new VCFFormatHeaderLine("AN_Orig", 1, VCFHeaderLineType.Integer, "Original AN"));
        }
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

        if (MENDELIAN_VIOLATIONS) {
            if ( FAMILY_STRUCTURE_FILE != null) {
                try {
                    for ( final String line : new XReadLines( FAMILY_STRUCTURE_FILE ) ) {
                        MendelianViolation mv = new MendelianViolation(line, MENDELIAN_VIOLATION_QUAL_THRESHOLD);
                        if (samples.contains(mv.getSampleChild()) &&  samples.contains(mv.getSampleDad()) && samples.contains(mv.getSampleMom()))
                            mvSet.add(mv);
                    }
                } catch ( FileNotFoundException e ) {
                    throw new UserException.CouldNotReadInputFile(AF_FILE, e);
                }
                if (outMVFile != null)
                    try {
                        outMVFileStream = new PrintStream(outMVFile);
                    }
                    catch (FileNotFoundException e) {
                        throw new UserException.CouldNotCreateOutputFile(outMVFile, "Can't open output file", e);   }
            }
            else
                mvSet.add(new MendelianViolation(getToolkit(), MENDELIAN_VIOLATION_QUAL_THRESHOLD));
        }
        else if (!FAMILY_STRUCTURE.isEmpty()) {
            mvSet.add(new MendelianViolation(FAMILY_STRUCTURE, MENDELIAN_VIOLATION_QUAL_THRESHOLD));
            MENDELIAN_VIOLATIONS = true;
        }

        SELECT_RANDOM_NUMBER = numRandom > 0;
        if (SELECT_RANDOM_NUMBER) {
            logger.info("Selecting " + numRandom + " variants at random from the variant track");
            variantArray = new RandomVariantStructure[numRandom];
        }

        SELECT_RANDOM_FRACTION = fractionRandom > 0;
        if (SELECT_RANDOM_FRACTION) logger.info("Selecting approximately " + fractionRandom + "% of the variants at random from the variant track");


        if (KEEP_AF_SPECTRUM) {
            try {
                afBreakpoints = new ArrayList<Double>();
                afBoosts = new ArrayList<Double>();
                logger.info("Reading in AF boost table...");
                boolean firstLine = false;
                for ( final String line : new XReadLines( AF_FILE ) ) {
                    if (!firstLine) {
                        firstLine = true;
                        continue;
                    }
                    final String[] vals = line.split(" ");

                    double bkp = Double.valueOf(vals[0]);
                    double afb = Double.valueOf(vals[1]);
                    afBreakpoints.add(bkp);
                    afBoosts.add(afb);

                }
                bkDelta = afBreakpoints.get(0);
            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotReadInputFile(AF_FILE, e);
            }

        }
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

        if ( vcs == null || vcs.size() == 0) {
            return 0;
        }

        for (VariantContext vc : vcs) {
            if (MENDELIAN_VIOLATIONS) {
                boolean foundMV = false;
                for (MendelianViolation mv : mvSet) {
                    if (mv.isViolation(vc)) {
                        foundMV = true;
                        //System.out.println(vc.toString());
                        if (outMVFile != null)
                            outMVFileStream.format("MV@%s:%d. REF=%s, ALT=%s, AC=%d, momID=%s, dadID=%s, childID=%s, momG=%s, momGL=%s, dadG=%s, dadGL=%s, " +
                                "childG=%s childGL=%s\n",vc.getChr(), vc.getStart(),
                                vc.getReference().getDisplayString(), vc.getAlternateAllele(0).getDisplayString(),  vc.getChromosomeCount(vc.getAlternateAllele(0)),
                                mv.getSampleMom(), mv.getSampleDad(), mv.getSampleChild(),
                                vc.getGenotype(mv.getSampleMom()).toBriefString(), vc.getGenotype(mv.getSampleMom()).getLikelihoods().getAsString(),
                                vc.getGenotype(mv.getSampleDad()).toBriefString(), vc.getGenotype(mv.getSampleMom()).getLikelihoods().getAsString(),
                                vc.getGenotype(mv.getSampleChild()).toBriefString(),vc.getGenotype(mv.getSampleChild()).getLikelihoods().getAsString()  );
                    }
                }

                if (!foundMV)
                    break;
            }
            if (DISCORDANCE_ONLY) {
                Collection<VariantContext> compVCs = tracker.getVariantContexts(ref, discordanceRodName, null, context.getLocation(), true, false);
                if (!isDiscordant(vc, compVCs))
                    return 0;
            }
            if (CONCORDANCE_ONLY) {
                Collection<VariantContext> compVCs = tracker.getVariantContexts(ref, concordanceRodName, null, context.getLocation(), true, false);
                if (!isConcordant(vc, compVCs))
                    return 0;
            }

            // TODO - add ability to also select MNPs
            // TODO - move variant selection arguments to the engine so other walkers can also do this
            if (SELECT_INDELS && !(vc.isIndel() || vc.isMixed()))
                continue;

            if (SELECT_SNPS && !vc.isSNP())
                continue;

            VariantContext sub = subsetRecord(vc, samples);
            if ( (sub.isPolymorphic() || !EXCLUDE_NON_VARIANTS) && (!sub.isFiltered() || !EXCLUDE_FILTERED) ) {
                for ( VariantContextUtils.JexlVCMatchExp jexl : jexls ) {
                    if ( !VariantContextUtils.match(sub, jexl) ) {
                        return 0;
                    }
                }
                if (SELECT_RANDOM_NUMBER) {
                    randomlyAddVariant(++variantNumber, sub, ref.getBase());
                }
                else if (!SELECT_RANDOM_FRACTION || (!KEEP_AF_SPECTRUM && GenomeAnalysisEngine.getRandomGenerator().nextDouble() < fractionRandom)) {
                    vcfWriter.add(sub, ref.getBase());
                }
                else {
                    if (SELECT_RANDOM_FRACTION && KEEP_AF_SPECTRUM ) {
                        // ok we have a comp VC and we need to match the AF spectrum of inputAFRodName.
                        // We then pick a variant with probablity AF*desiredFraction
                        if ( sub.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) )  {
                            String afo = sub.getAttributeAsString(VCFConstants.ALLELE_FREQUENCY_KEY);

                            double af;
                            double afBoost = 1.0;
                            if (afo.contains(",")) {
                                String[] afs = afo.split(",");
                                afs[0] = afs[0].substring(1,afs[0].length());
                                afs[afs.length-1] = afs[afs.length-1].substring(0,afs[afs.length-1].length()-1);

                                double[] afd = new double[afs.length];

                                for (int k=0; k < afd.length; k++)
                                    afd[k] = Double.valueOf(afs[k]);

                                af = MathUtils.arrayMax(afd);
                                //af = Double.valueOf(afs[0]);

                            }
                            else
                                af = Double.valueOf(afo);

                            // now boost af by table read from file if desired
                            //double bkpt = 0.0;
                            int bkidx = 0;
                            if (!afBreakpoints.isEmpty()) {
                                for ( Double bkpt : afBreakpoints) {
                                    if (af < bkpt + bkDelta)
                                        break;
                                    else bkidx++;
                                }
                                if (bkidx >=afBoosts.size())
                                    bkidx = afBoosts.size()-1;
                                afBoost = afBoosts.get(bkidx);
                                //System.out.formatPrin("af:%f bkidx:%d afboost:%f\n",af,bkidx,afBoost);



                            }

                            //System.out.format("%s .. %4.4f\n",afo.toString(), af);
                            if (GenomeAnalysisEngine.getRandomGenerator().nextDouble() < fractionRandom * afBoost *   afBoost)
                                vcfWriter.add(sub, ref.getBase());
                        }


                    }
                }
            }

        }

        return 1;
    }

    /**
     * Checks if vc has a variant call for (at least one of) the samples.
     * @param vc the variant rod VariantContext. Here, the variant is the dataset you're looking for discordances to (e.g. HapMap)
     * @param compVCs the comparison VariantContext (discordance
     * @return
     */
    private boolean isDiscordant (VariantContext vc, Collection<VariantContext> compVCs) {
        if (vc == null)
            return false;

        // if we're not looking at specific samples then the absense of a compVC means discordance
        if (NO_SAMPLES_SPECIFIED && (compVCs == null || compVCs.isEmpty()))
            return true;

        // check if we find it in the variant rod
        Map<String, Genotype> genotypes = vc.getGenotypes(samples);
        for (Genotype g : genotypes.values()) {
            if (sampleHasVariant(g)) {
                // There is a variant called (or filtered with not exclude filtered option set) that is not HomRef for at least one of the samples.
                if (compVCs == null)
                    return true;
                // Look for this sample in the all vcs of the comp ROD track.
                boolean foundVariant = false;
                for (VariantContext compVC : compVCs) {
                    if (sampleHasVariant(compVC.getGenotype(g.getSampleName()))) {
                        foundVariant = true;
                        break;
                    }
                }
                // if (at least one sample) was not found in all VCs of the comp ROD, we have discordance
                if (!foundVariant)
                    return true;
            }
        }
        return false; // we only get here if all samples have a variant in the comp rod.
    }

    private boolean isConcordant (VariantContext vc, Collection<VariantContext> compVCs) {
        if (vc == null || compVCs == null || compVCs.isEmpty())
            return false;

        // if we're not looking for specific samples then the fact that we have both VCs is enough to call it concordant.
        if (NO_SAMPLES_SPECIFIED)
            return true;

        // make a list of all samples contained in this variant VC that are being tracked by the user command line arguments.
        Set<String> variantSamples = vc.getSampleNames();
        variantSamples.retainAll(samples);

        // check if we can find all samples from the variant rod in the comp rod.
        for (String sample : variantSamples) {
            boolean foundSample = false;
            for (VariantContext compVC : compVCs) {
                Genotype varG = vc.getGenotype(sample);
                Genotype compG = compVC.getGenotype(sample);
                if (haveSameGenotypes(varG, compG)) {
                    foundSample = true;
                    break;
                }
            }
            // if at least one sample doesn't have the same genotype, we don't have concordance
            if (!foundSample) {
                return false;
            }
        }
        return true;
    }

    private boolean sampleHasVariant(Genotype g) {
        return (g !=null && !g.isHomRef() && (g.isCalled() || (g.isFiltered() && !EXCLUDE_FILTERED)));
    }

    private boolean haveSameGenotypes(Genotype g1, Genotype g2) {
        if ((g1.isCalled() && g2.isFiltered()) ||
                (g2.isCalled() && g1.isFiltered()) ||
                (g1.isFiltered() && g2.isFiltered() && EXCLUDE_FILTERED))
            return false;

        List<Allele> a1s = g1.getAlleles();
        List<Allele> a2s = g2.getAlleles();
        return (a1s.containsAll(a2s) && a2s.containsAll(a1s));
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
                positionToPrint = nextCircularPosition(positionToPrint);
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

        VariantContext sub = vc.subContextFromGenotypes(genotypes, vc.getAlleles());

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


        if (KEEP_ORIGINAL_CHR_COUNTS) {
            if ( attributes.containsKey(VCFConstants.ALLELE_COUNT_KEY) )
                attributes.put("AC_Orig",attributes.get(VCFConstants.ALLELE_COUNT_KEY));
            if ( attributes.containsKey(VCFConstants.ALLELE_FREQUENCY_KEY) )
                attributes.put("AF_Orig",attributes.get(VCFConstants.ALLELE_FREQUENCY_KEY));
            if ( attributes.containsKey(VCFConstants.ALLELE_NUMBER_KEY) )
                attributes.put("AN_Orig",attributes.get(VCFConstants.ALLELE_NUMBER_KEY));

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
                positionToAdd = nextCircularPosition(positionToAdd);
            }
        }
    }

    private int nextCircularPosition(int cur) {
        if ((cur + 1) == variantArray.length)
            return 0;
        return cur + 1;
    }
}

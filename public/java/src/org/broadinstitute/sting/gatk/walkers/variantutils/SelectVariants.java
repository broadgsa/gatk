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

import org.apache.poi.hpsf.Variant;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.Genotype;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

/**
 * Selects variants from a VCF source.
 *
 * <p>
 * Often, a VCF containing many samples and/or variants will need to be subset in order to facilitate certain analyses
 * (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain
 * requirements, displaying just a few samples in a browser like IGV, etc.). SelectVariants can be used for this purpose.
 * Given a single VCF file, one or more samples can be extracted from the file (based on a complete sample name or a
 * pattern match).  Variants can be further selected by specifying criteria for inclusion, i.e. "DP > 1000" (depth of
 * coverage greater than 1000x), "AF < 0.25" (sites with allele frequency less than 0.25).  These JEXL expressions are
 * documented in the Using JEXL expressions section (http://www.broadinstitute.org/gsa/wiki/index.php/Using_JEXL_expressions).
 * One can optionally include concordance or discordance tracks for use in selecting overlapping variants.
 *
 * <h2>Input</h2>
 * <p>
 * A variant set to select from.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A selected VCF.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * Select two samples out of a VCF with many samples:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_A_PARC \
 *   -sn SAMPLE_B_ACTG
 *
 * Select two samples and any sample that matches a regular expression:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_1_PARC \
 *   -sn SAMPLE_1_ACTG \
 *   -se 'SAMPLE.+PARC'
 *
 * Select any sample that matches a regular expression and sites where the QD annotation is more than 10:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -se 'SAMPLE.+PARC'
 *   -select "QD > 10.0"
 *
 * Select a sample and exclude non-variant loci and filtered loci:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_1_ACTG \
 *   -env \
 *   -ef
 *
 * Select a sample and restrict the output vcf to a set of intervals:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -L /path/to/my.interval_list \
 *   -sn SAMPLE_1_ACTG
 *
 * Select all calls missed in my vcf, but present in HapMap (useful to take a look at why these variants weren't called by this dataset):
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant hapmap.vcf \
 *   --discordance myCalls.vcf
 *   -o output.vcf \
 *   -sn mySample
 *
 * Select all calls made by both myCalls and hisCalls (useful to take a look at what is consistent between the two callers):
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant myCalls.vcf \
 *   --concordance hisCalls.vcf
 *   -o output.vcf \
 *   -sn mySample
 *
 * Generating a VCF of all the variants that are mendelian violations:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -family NA12891+NA12892=NA12878 \
 *   -mvq 50 \
 *   -o violations.vcf
 *
 * Creating a sample of exactly 1000 variants randomly chosen with equal probability from the variant VCF:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -number 1000
 *
 * Creating a set with 50% of the total number of variants in the variant VCF:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -fraction 0.5
 *
 * Select only indels from a VCF:
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -selectType INDEL
 *
 * Select only multi-allelic SNPs and MNPs from a VCF (i.e. SNPs with more than one allele listed in the ALT column):
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -selectType SNP -selectType MNP \
 *   -restrictAllelesTo MULTIALLELIC
 *
 * </pre>
 *
 */
public class SelectVariants extends RodWalker<Integer, Integer> {
    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * A site is considered discordant if there exists some sample in the variant track that has a non-reference genotype
     * and either the site isn't present in this track, the sample isn't present in this track,
     * or the sample is called reference in this track.
     */
    @Input(fullName="discordance", shortName = "disc", doc="Output variants that were not called in this comparison track", required=false)
    private RodBinding<VariantContext> discordanceTrack;

    /**
     * A site is considered concordant if (1) we are not looking for specific samples and there is a variant called
     * in both the variant and concordance tracks or (2) every sample present in the variant track is present in the
     * concordance track and they have the sample genotype call.
     */
    @Input(fullName="concordance", shortName = "conc", doc="Output variants that were also called in this comparison track", required=false)
    private RodBinding<VariantContext> concordanceTrack;

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    @Argument(fullName="sample_name", shortName="sn", doc="Include genotypes from this sample. Can be specified multiple times", required=false)
    public Set<String> sampleNames = new HashSet<String>(0);

    @Argument(fullName="sample_expressions", shortName="se", doc="Regular expression to select many samples from the ROD tracks provided. Can be specified multiple times", required=false)
    public Set<String> sampleExpressions ;

    @Input(fullName="sample_file", shortName="sf", doc="File containing a list of samples (one per line) to include. Can be specified multiple times", required=false)
    public Set<File> sampleFiles;

    /**
     * Note that sample exclusion takes precedence over inclusion, so that if a sample is in both lists it will be excluded.
     */
    @Argument(fullName="exclude_sample_name", shortName="xl_sn", doc="Exclude genotypes from this sample. Can be specified multiple times", required=false)
    public Set<String> XLsampleNames = new HashSet<String>(0);

    /**
     * Note that sample exclusion takes precedence over inclusion, so that if a sample is in both lists it will be excluded.
     */
    @Input(fullName="exclude_sample_file", shortName="xl_sf", doc="File containing a list of samples (one per line) to exclude. Can be specified multiple times", required=false)
    public Set<File> XLsampleFiles = new HashSet<File>(0);

    /**
     * Note that these expressions are evaluated *after* the specified samples are extracted and the INFO field annotations are updated.
     */
    @Argument(shortName="select", doc="One or more criteria to use when selecting the data", required=false)
    public ArrayList<String> SELECT_EXPRESSIONS = new ArrayList<String>();

    @Argument(fullName="excludeNonVariants", shortName="env", doc="Don't include loci found to be non-variant after the subsetting procedure", required=false)
    private boolean EXCLUDE_NON_VARIANTS = false;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered loci in the analysis", required=false)
    private boolean EXCLUDE_FILTERED = false;


    /**
     * When this argument is used, we can choose to include only multiallelic or biallelic sites, depending on how many alleles are listed in the ALT column of a vcf.
     * For example, a multiallelic record such as:
     * 1    100 .   A   AAA,AAAAA
     * will be excluded if "-restrictAllelesTo BIALLELIC" is included, because there are two alternate alleles, whereas a record such as:
     * 1    100 .   A  T
     * will be included in that case, but would be excluded if "-restrictAllelesTo MULTIALLELIC
     */
    @Argument(fullName="restrictAllelesTo", shortName="restrictAllelesTo", doc="Select only variants of a particular allelicity. Valid options are ALL (default), MULTIALLELIC or BIALLELIC", required=false)
    private  NumberAlleleRestriction alleleRestriction = NumberAlleleRestriction.ALL;

    @Argument(fullName="keepOriginalAC", shortName="keepOriginalAC", doc="Don't update the AC, AF, or AN values in the INFO field after selecting", required=false)
    private boolean KEEP_ORIGINAL_CHR_COUNTS = false;

    @Hidden
    @Argument(fullName="family_structure_file", shortName="familyFile", doc="use -family unless you know what you're doing", required=false)
    private File FAMILY_STRUCTURE_FILE = null;

    /**
     * String formatted as dad+mom=child where these parameters determine which sample names are examined.
     */
    @Argument(fullName="family_structure", shortName="family", doc="string formatted as dad+mom=child where these parameters determine which sample names are examined", required=false)
    private String FAMILY_STRUCTURE = "";

    /**
     * This activates the mendelian violation module that will select all variants that correspond to a mendelian violation following the rules given by the family structure.
     */
    @Argument(fullName="mendelianViolation", shortName="mv", doc="output mendelian violation sites only", required=false)
    private Boolean MENDELIAN_VIOLATIONS = false;

    @Argument(fullName="mendelianViolationQualThreshold", shortName="mvq", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    private double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 0;

    /**
     * Variants are kept in memory to guarantee that exactly n variants will be chosen randomly, so use it only for a reasonable
     * number of variants.  Use --select_random_fraction for larger numbers of variants.
     */
    @Argument(fullName="select_random_number", shortName="number", doc="Selects a number of variants at random from the variant track", required=false)
    private int numRandom = 0;

    /**
     * This routine is based on probability, so the final result is not guaranteed to carry the exact fraction.  Can be used for large fractions.
     */
    @Argument(fullName="select_random_fraction", shortName="fraction", doc="Selects a fraction (a number between 0 and 1) of the total variants at random from the variant track", required=false)
    private double fractionRandom = 0;

    /**
      * This argument select particular kinds of variants out of a list. If left empty, there is no type selection and all variant types are considered for other selection criteria.
     * When specified one or more times, a particular type of variant is selected.
     *
      */
    @Argument(fullName="selectTypeToInclude", shortName="selectType", doc="Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times", required=false)
    private List<VariantContext.Type> TYPES_TO_INCLUDE = new ArrayList<VariantContext.Type>();


    @Hidden
    @Argument(fullName="outMVFile", shortName="outMVFile", doc="", required=false)
    private String outMVFile = null;

    /* Private class used to store the intermediate variants in the integer random selection process */
    private class RandomVariantStructure {
        private VariantContext vc;

        RandomVariantStructure(VariantContext vcP) {
            vc = vcP;
        }

        public void set (VariantContext vcP) {
            vc = vcP;
        }

    }

    public enum NumberAlleleRestriction {
           ALL,
           BIALLELIC,
           MULTIALLELIC
    }

    private ArrayList<VariantContext.Type> selectedTypes = new ArrayList<VariantContext.Type>();
    private ArrayList<String> selectNames = new ArrayList<String>();
    private List<VariantContextUtils.JexlVCMatchExp> jexls = null;

    private TreeSet<String> samples = new TreeSet<String>();
    private boolean NO_SAMPLES_SPECIFIED = false;

    private boolean DISCORDANCE_ONLY = false;
    private boolean CONCORDANCE_ONLY = false;

    private Set<MendelianViolation> mvSet = new HashSet<MendelianViolation>();


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
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

        Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
        Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

        // first, add any requested samples
        samples.addAll(samplesFromFile);
        samples.addAll(samplesFromExpressions);
        samples.addAll(sampleNames);

        // if none were requested, we want all of them
        if ( samples.isEmpty() ) {
            samples.addAll(vcfSamples);
            NO_SAMPLES_SPECIFIED = true;
        }

        // now, exclude any requested samples
        Collection<String> XLsamplesFromFile = SampleUtils.getSamplesFromFiles(XLsampleFiles);
        samples.removeAll(XLsamplesFromFile);
        samples.removeAll(XLsampleNames);

        if ( samples.size() == 0 && !NO_SAMPLES_SPECIFIED )
            throw new UserException("All samples requested to be included were also requested to be excluded.");

        for ( String sample : samples )
            logger.info("Including sample '" + sample + "'");



        // if user specified types to include, add these, otherwise, add all possible variant context types to list of vc types to include
        if (TYPES_TO_INCLUDE.isEmpty()) {

            for (VariantContext.Type t : VariantContext.Type.values())
                selectedTypes.add(t);

        }
        else {
            for (VariantContext.Type t : TYPES_TO_INCLUDE)
                selectedTypes.add(t);

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
        DISCORDANCE_ONLY = discordanceTrack.isBound();
        if (DISCORDANCE_ONLY) logger.info("Selecting only variants discordant with the track: " + discordanceTrack.getName());

        CONCORDANCE_ONLY = concordanceTrack.isBound();
        if (CONCORDANCE_ONLY) logger.info("Selecting only variants concordant with the track: " + concordanceTrack.getName());

        if (MENDELIAN_VIOLATIONS) {
            if ( FAMILY_STRUCTURE_FILE != null) {
                try {
                    for ( final String line : new XReadLines( FAMILY_STRUCTURE_FILE ) ) {
                        MendelianViolation mv = new MendelianViolation(line, MENDELIAN_VIOLATION_QUAL_THRESHOLD);
                        if (samples.contains(mv.getSampleChild()) &&  samples.contains(mv.getSampleDad()) && samples.contains(mv.getSampleMom()))
                            mvSet.add(mv);
                    }
                } catch ( FileNotFoundException e ) {
                    throw new UserException.CouldNotReadInputFile(FAMILY_STRUCTURE_FILE, e);
                }
                if (outMVFile != null)
                    try {
                        outMVFileStream = new PrintStream(outMVFile);
                    }
                    catch (FileNotFoundException e) {
                        throw new UserException.CouldNotCreateOutputFile(outMVFile, "Can't open output file", e);   }
            }
            else
                mvSet.add(new MendelianViolation(FAMILY_STRUCTURE, MENDELIAN_VIOLATION_QUAL_THRESHOLD));
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
        if (SELECT_RANDOM_FRACTION) logger.info("Selecting approximately " + 100.0*fractionRandom + "% of the variants at random from the variant track");



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

        Collection<VariantContext> vcs = tracker.getValues(variantCollection.variants, context.getLocation());

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
                Collection<VariantContext> compVCs = tracker.getValues(discordanceTrack, context.getLocation());
                if (!isDiscordant(vc, compVCs))
                    return 0;
            }
            if (CONCORDANCE_ONLY) {
                Collection<VariantContext> compVCs = tracker.getValues(concordanceTrack, context.getLocation());
                if (!isConcordant(vc, compVCs))
                    return 0;
            }

            if (alleleRestriction.equals(NumberAlleleRestriction.BIALLELIC) && !vc.isBiallelic())
                continue;

            if (alleleRestriction.equals(NumberAlleleRestriction.MULTIALLELIC) && vc.isBiallelic())
                continue;

            if (!selectedTypes.contains(vc.getType()))
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
                else if (!SELECT_RANDOM_FRACTION || ( GenomeAnalysisEngine.getRandomGenerator().nextDouble() < fractionRandom)) {
                    vcfWriter.add(sub);
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

        // if we're not looking at specific samples then the absence of a compVC means discordance
        if (NO_SAMPLES_SPECIFIED)
            return (compVCs == null || compVCs.isEmpty());

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
                vcfWriter.add(variantArray[positionToPrint].vc);
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

        // if we have fewer alternate alleles in the selected VC than in the original VC, we need to strip out the GL/PLs (because they are no longer accurate)
        if ( vc.getAlleles().size() != sub.getAlleles().size() )
            sub = VariantContext.modifyGenotypes(sub, VariantContextUtils.stripPLs(vc.getGenotypes()));

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
            variantArray[nVariantsAdded++] = new RandomVariantStructure(vc);

        else {
            double v = GenomeAnalysisEngine.getRandomGenerator().nextDouble();
            double t = (1.0/(rank-numRandom+1));
            if ( v < t) {
                variantArray[positionToAdd].set(vc);
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

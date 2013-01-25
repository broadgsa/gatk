/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.GenomeAnalysisEngine;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.annotator.ChromosomeCountConstants;
import org.broadinstitute.sting.utils.MendelianViolation;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.FileNotFoundException;
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
 *   -bed family.ped \
 *   -mvq 50 \
 *   -o violations.vcf
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
@DocumentedGATKFeature( groupName = "Variant Evaluation and Manipulation Tools", extraDocs = {CommandLineGATK.class} )
public class SelectVariants extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * A site is considered discordant if there exists some sample in the variant track that has a non-reference genotype
     * and either the site isn't present in this track, the sample isn't present in this track,
     * or the sample is called reference in this track.
     */
    @Input(fullName="discordance", shortName = "disc", doc="Output variants that were not called in this comparison track", required=false)
    protected RodBinding<VariantContext> discordanceTrack;

    /**
     * A site is considered concordant if (1) we are not looking for specific samples and there is a variant called
     * in both the variant and concordance tracks or (2) every sample present in the variant track is present in the
     * concordance track and they have the sample genotype call.
     */
    @Input(fullName="concordance", shortName = "conc", doc="Output variants that were also called in this comparison track", required=false)
    protected RodBinding<VariantContext> concordanceTrack;

    @Output(doc="File to which variants should be written",required=true)
    protected VariantContextWriter vcfWriter = null;

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
    protected boolean EXCLUDE_NON_VARIANTS = false;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered loci in the analysis", required=false)
    protected boolean EXCLUDE_FILTERED = false;

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

    @Argument(fullName="keepOriginalAC", shortName="keepOriginalAC", doc="Store the original AC, AF, and AN values in the INFO field after selecting (using keys AC_Orig, AF_Orig, and AN_Orig)", required=false)
    private boolean KEEP_ORIGINAL_CHR_COUNTS = false;

    /**
     * This activates the mendelian violation module that will select all variants that correspond to a mendelian violation following the rules given by the family structure.
     */
    @Argument(fullName="mendelianViolation", shortName="mv", doc="output mendelian violation sites only", required=false)
    private Boolean MENDELIAN_VIOLATIONS = false;

    @Argument(fullName="mendelianViolationQualThreshold", shortName="mvq", doc="Minimum genotype QUAL score for each trio member required to accept a site as a violation", required=false)
    protected double MENDELIAN_VIOLATION_QUAL_THRESHOLD = 0;

    /**
     * This routine is based on probability, so the final result is not guaranteed to carry the exact fraction.  Can be used for large fractions.
     */
    @Argument(fullName="select_random_fraction", shortName="fraction", doc="Selects a fraction (a number between 0 and 1) of the total variants at random from the variant track", required=false)
    protected double fractionRandom = 0;

    @Argument(fullName="remove_fraction_genotypes", shortName="fractionGenotypes", doc="Selects a fraction (a number between 0 and 1) of the total genotypes at random from the variant track and sets them to nocall", required=false)
    protected double fractionGenotypes = 0;

    /**
     * This argument select particular kinds of variants out of a list. If left empty, there is no type selection and all variant types are considered for other selection criteria.
     * When specified one or more times, a particular type of variant is selected.
     *
     */
    @Argument(fullName="selectTypeToInclude", shortName="selectType", doc="Select only a certain type of variants from the input file. Valid types are INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION. Can be specified multiple times", required=false)
    private List<VariantContext.Type> TYPES_TO_INCLUDE = new ArrayList<VariantContext.Type>();

    /**
     * If provided, we will only include variants whose ID field is present in this list of ids.  The matching
     * is exact string matching.  The file format is just one ID per line
     *
     */
    @Argument(fullName="keepIDs", shortName="IDs", doc="Only emit sites whose ID is found in this file (one ID per line)", required=false)
    private File rsIDFile = null;


    @Hidden
    @Argument(fullName="fullyDecode", doc="If true, the incoming VariantContext will be fully decoded", required=false)
    private boolean fullyDecode = false;

    @Hidden
    @Argument(fullName="forceGenotypesDecode", doc="If true, the incoming VariantContext will have its genotypes forcibly decoded by computing AC across all genotypes.  For efficiency testing only", required=false)
    private boolean forceGenotypesDecode = false;

    @Hidden
    @Argument(fullName="justRead", doc="If true, we won't actually write the output file.  For efficiency testing only", required=false)
    private boolean justRead = false;

    @Argument(doc="indel size select",required=false,fullName="maxIndelSize")
    private int maxIndelSize = Integer.MAX_VALUE;

    @Argument(doc="Allow a samples other than those in the VCF to be specified on the command line. These samples will be ignored.",required=false,fullName="ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES")
    private boolean ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES = false;


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

    private MendelianViolation mv;


    /* variables used by the SELECT RANDOM modules */
    private boolean SELECT_RANDOM_FRACTION = false;

    //Random number generator for the genotypes to remove
    private Random randomGenotypes = new Random();

    private Set<String> IDsToKeep = null;
    private Map<String, VCFHeader> vcfRods;

    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        // Get list of samples to include in the output
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<String>(SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

        Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
        Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

        // first, check overlap between requested and present samples
        Set<String> commandLineUniqueSamples = new HashSet<String>(samplesFromFile.size()+samplesFromExpressions.size()+sampleNames.size());
        commandLineUniqueSamples.addAll(samplesFromFile);
        commandLineUniqueSamples.addAll(samplesFromExpressions);
        commandLineUniqueSamples.addAll(sampleNames);
        commandLineUniqueSamples.removeAll(vcfSamples);

        // second, add the requested samples
        samples.addAll(sampleNames);
        samples.addAll(samplesFromExpressions);
        samples.addAll(samplesFromFile);

        logger.debug(Utils.join(",",commandLineUniqueSamples));

        if ( commandLineUniqueSamples.size() > 0 && ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES ) {
            logger.warn("Samples present on command line input that are not present in the VCF. These samples will be ignored.");
            samples.removeAll(commandLineUniqueSamples);
        } else if (commandLineUniqueSamples.size() > 0 ) {
            throw new UserException.BadInput(String.format("%s%n%n%s%n%n%s%n%n%s",
                    "Samples entered on command line (through -sf or -sn) that are not present in the VCF.",
                    "A list of these samples:",
                    Utils.join(",",commandLineUniqueSamples),
                    "To ignore these samples, run with --ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES"));
        }


        // if none were requested, we want all of them
        if ( samples.isEmpty() ) {
            samples.addAll(vcfSamples);
            NO_SAMPLES_SPECIFIED = true;
        }

        // now, exclude any requested samples
        Collection<String> XLsamplesFromFile = SampleUtils.getSamplesFromFiles(XLsampleFiles);
        samples.removeAll(XLsamplesFromFile);
        samples.removeAll(XLsampleNames);
        NO_SAMPLES_SPECIFIED = NO_SAMPLES_SPECIFIED && XLsampleNames.isEmpty();

        if ( samples.size() == 0 && !NO_SAMPLES_SPECIFIED )
            throw new UserException("All samples requested to be included were also requested to be excluded.");

        if ( ! NO_SAMPLES_SPECIFIED )
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
            headerLines.add(new VCFInfoHeaderLine("AC_Orig", 1, VCFHeaderLineType.Integer, "Original AC"));
            headerLines.add(new VCFInfoHeaderLine("AF_Orig", 1, VCFHeaderLineType.Float, "Original AF"));
            headerLines.add(new VCFInfoHeaderLine("AN_Orig", 1, VCFHeaderLineType.Integer, "Original AN"));
        }
        headerLines.addAll(Arrays.asList(ChromosomeCountConstants.descriptions));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));

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
            mv = new MendelianViolation(MENDELIAN_VIOLATION_QUAL_THRESHOLD,false,true);
        }

        SELECT_RANDOM_FRACTION = fractionRandom > 0;
        if (SELECT_RANDOM_FRACTION) logger.info("Selecting approximately " + 100.0*fractionRandom + "% of the variants at random from the variant track");

        /** load in the IDs file to a hashset for matching */
        if ( rsIDFile != null ) {
            IDsToKeep = new HashSet<String>();
            try {
                for ( final String line : new XReadLines(rsIDFile).readLines() ) {
                    IDsToKeep.add(line.trim());
                }
                logger.info("Selecting only variants with one of " + IDsToKeep.size() + " IDs from " + rsIDFile);
            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotReadInputFile(rsIDFile, e);
            }
        }

        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
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
            // an option for performance testing only
            if ( fullyDecode )
                vc = vc.fullyDecode(vcfRods.get(vc.getSource()), getToolkit().lenientVCFProcessing() );

            // an option for performance testing only
            if ( forceGenotypesDecode ) {
                final int x = vc.getCalledChrCount();
                //logger.info("forceGenotypesDecode with getCalledChrCount() = " + );
            }

            if ( IDsToKeep != null && ! IDsToKeep.contains(vc.getID()) )
                continue;

            if (MENDELIAN_VIOLATIONS && mv.countViolations(this.getSampleDB().getFamilies(samples),vc) < 1)
                break;

            if (DISCORDANCE_ONLY) {
                Collection<VariantContext> compVCs = tracker.getValues(discordanceTrack, context.getLocation());
                if (!isDiscordant(vc, compVCs))
                    continue;
            }
            if (CONCORDANCE_ONLY) {
                Collection<VariantContext> compVCs = tracker.getValues(concordanceTrack, context.getLocation());
                if (!isConcordant(vc, compVCs))
                    continue;
            }

            if (alleleRestriction.equals(NumberAlleleRestriction.BIALLELIC) && !vc.isBiallelic())
                continue;

            if (alleleRestriction.equals(NumberAlleleRestriction.MULTIALLELIC) && vc.isBiallelic())
                continue;

            if (!selectedTypes.contains(vc.getType()))
                continue;

            if ( badIndelSize(vc) )
                continue;

            VariantContext sub = subsetRecord(vc, EXCLUDE_NON_VARIANTS);

            if ( (!EXCLUDE_NON_VARIANTS || sub.isPolymorphicInSamples()) && (!EXCLUDE_FILTERED || !sub.isFiltered()) ) {
                boolean failedJexlMatch = false;
                for ( VariantContextUtils.JexlVCMatchExp jexl : jexls ) {
                    if ( !VariantContextUtils.match(sub, jexl) ) {
                        failedJexlMatch = true;
                        break;
                    }
                }
                if ( !failedJexlMatch &&
                        !justRead &&
                        ( !SELECT_RANDOM_FRACTION || GenomeAnalysisEngine.getRandomGenerator().nextDouble() < fractionRandom ) ) {
                    vcfWriter.add(sub);
                }
            }
        }

        return 1;
    }

    private boolean badIndelSize(final VariantContext vc) {
        List<Integer> lengths = vc.getIndelLengths();
        if ( lengths == null )
            return false; // VC does not harbor indel
        for ( Integer indelLength : vc.getIndelLengths() ) {
            if ( indelLength > maxIndelSize )
                return true;
        }

        return false;
    }

    /**
     * Checks if vc has a variant call for (at least one of) the samples.
     * @param vc the variant rod VariantContext. Here, the variant is the dataset you're looking for discordances to (e.g. HapMap)
     * @param compVCs the comparison VariantContext (discordance
     * @return true if is discordant
     */
    private boolean isDiscordant (VariantContext vc, Collection<VariantContext> compVCs) {
        if (vc == null)
            return false;

        // if we're not looking at specific samples then the absence of a compVC means discordance
        if (NO_SAMPLES_SPECIFIED)
            return (compVCs == null || compVCs.isEmpty());

        // check if we find it in the variant rod
        GenotypesContext genotypes = vc.getGenotypes(samples);
        for (final Genotype g : genotypes) {
            if (sampleHasVariant(g)) {
                // There is a variant called (or filtered with not exclude filtered option set) that is not HomRef for at least one of the samples.
                if (compVCs == null)
                    return true;
                // Look for this sample in the all vcs of the comp ROD track.
                boolean foundVariant = false;
                for (VariantContext compVC : compVCs) {
                    if (haveSameGenotypes(g, compVC.getGenotype(g.getSampleName()))) {
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

    private boolean haveSameGenotypes(final Genotype g1, final Genotype g2) {
        if ( g1 == null || g2 == null )
            return false;

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

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return lhs + rhs;
    }

    public void onTraversalDone(Integer result) {
        logger.info(result + " records processed.");
    }



    /**
     * Helper method to subset a VC record, modifying some metadata stored in the INFO field (i.e. AN, AC, AF).
     *
     * @param vc       the VariantContext record to subset
     * @return the subsetted VariantContext
     */
    private VariantContext subsetRecord(final VariantContext vc, final boolean excludeNonVariants) {
        if ( NO_SAMPLES_SPECIFIED || samples.isEmpty() )
            return vc;

        final VariantContext sub = vc.subContextFromSamples(samples, excludeNonVariants); // strip out the alternate alleles that aren't being used

        VariantContextBuilder builder = new VariantContextBuilder(sub);

        GenotypesContext newGC = sub.getGenotypes();

        // if we have fewer alternate alleles in the selected VC than in the original VC, we need to strip out the GL/PLs and AD (because they are no longer accurate)
        if ( vc.getAlleles().size() != sub.getAlleles().size() )
            newGC = VariantContextUtils.stripPLsAndAD(sub.getGenotypes());

        // if we have fewer samples in the selected VC than in the original VC, we need to strip out the MLE tags
        if ( vc.getNSamples() != sub.getNSamples() ) {
            builder.rmAttribute(VCFConstants.MLE_ALLELE_COUNT_KEY);
            builder.rmAttribute(VCFConstants.MLE_ALLELE_FREQUENCY_KEY);
        }

        // Remove a fraction of the genotypes if needed
        if ( fractionGenotypes > 0 ){
            ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
            for ( Genotype genotype : newGC ) {
                //Set genotype to no call if it falls in the fraction.
                if(fractionGenotypes>0 && randomGenotypes.nextDouble()<fractionGenotypes){
                    List<Allele> alleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);
                    genotypes.add(new GenotypeBuilder(genotype).alleles(alleles).noGQ().make());
                }
                else{
                    genotypes.add(genotype);
                }
            }
            newGC = GenotypesContext.create(genotypes);
        }

        builder.genotypes(newGC);

        addAnnotations(builder, sub);

        return builder.make();
    }

    private void addAnnotations(final VariantContextBuilder builder, final VariantContext originalVC) {
        if ( fullyDecode ) return; // TODO -- annotations are broken with fully decoded data

        if (KEEP_ORIGINAL_CHR_COUNTS) {
            if ( originalVC.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) )
                builder.attribute("AC_Orig", originalVC.getAttribute(VCFConstants.ALLELE_COUNT_KEY));
            if ( originalVC.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) )
                builder.attribute("AF_Orig", originalVC.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY));
            if ( originalVC.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) )
                builder.attribute("AN_Orig", originalVC.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
        }

        VariantContextUtils.calculateChromosomeCounts(builder, false);

        boolean sawDP = false;
        int depth = 0;
        for (String sample : originalVC.getSampleNames()) {
            Genotype g = originalVC.getGenotype(sample);

            if ( ! g.isFiltered() ) {
                if ( g.hasDP() ) {
                    depth += g.getDP();
                    sawDP = true;
                }
            }
        }

        if ( sawDP )
            builder.attribute("DP", depth);
    }
}

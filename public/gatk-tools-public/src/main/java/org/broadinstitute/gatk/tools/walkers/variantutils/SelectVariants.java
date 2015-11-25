/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.variantutils;

import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.variant.ChromosomeCountConstants;
import org.broadinstitute.gatk.engine.samples.MendelianViolation;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVCFHeaderLines;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.text.XReadLines;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Select a subset of variants from a larger callset
 *
 * <p>
 * Often, a VCF containing many samples and/or variants will need to be subset in order to facilitate certain analyses
 * (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain
 * requirements, displaying just a few samples in a browser like IGV, etc.). SelectVariants can be used for this purpose.
 * </p>
 * <p>
 * There are many different options for selecting subsets of variants from a larger callset:
 * <ul>
 *     <li>Extract one or more samples from a callset based on either a complete sample name or a pattern match.</li>  
 *     <li>Specify criteria for inclusion that place thresholds on annotation values, e.g. "DP > 1000" (depth of
 * coverage greater than 1000x), "AF < 0.25" (sites with allele frequency less than 0.25). These criteria are written 
 * as "JEXL expressions", which are documented in the 
 * <a href="http://www.broadinstitute.org/gatk/guide/article?id=1255">article about using JEXL expressions</a>.</li>
 *     <li>Provide concordance or discordance tracks in order to include or exclude variants that are 
 * also present in other given callsets.</li>
 *     <li>Select variants based on criteria like their type 
 * (e.g. INDELs only), evidence of mendelian violation, filtering status, allelicity, and so on.</li>
 * </ul>
 * </p>
 * 
 * <p>There are also several options for recording the original values of certain annotations that are recalculated 
 * when a subsetting the new callset, trimming alleles, and so on.</p>
 * 
 * <h3>Input</h3>
 * <p>
 * A variant call set from which to select a subset.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A new VCF file containing the selected subset of variants.
 * </p>
 *
 * <h3>Usage examples</h3>
 * <h4>Select two samples out of a VCF with many samples</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_A_PARC \
 *   -sn SAMPLE_B_ACTG
 * </pre>
 *
 * <h4>Select two samples and any sample that matches a regular expression</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_1_PARC \
 *   -sn SAMPLE_1_ACTG \
 *   -se 'SAMPLE.+PARC'
 * </pre>
 *
 * <h4>Exclude two samples and any sample that matches a regular expression:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -xl_sn SAMPLE_1_PARC \
 *   -xl_sn SAMPLE_1_ACTG \
 *   -xl_se 'SAMPLE.+PARC'
 * </pre>
 *
 * <h4>Select any sample that matches a regular expression and sites where the QD annotation is more than 10:</h4>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -se 'SAMPLE.+PARC' \
 *   -select "QD > 10.0"
 * </pre>
 *
 * <h4>Select any sample that does not match a regular expression and sites where the QD annotation is more than 10:</h4>
 * <pre>
 * java  -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -se 'SAMPLE.+PARC' \
 *   -select "QD > 10.0"
 *   -invertSelect
 * </pre>
 *
 * <h4>Select a sample and exclude non-variant loci and filtered loci (trim remaining alleles by default):</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_1_ACTG \
 *   -env \
 *   -ef
 * </pre>
 *
 * <h4>Select a sample, subset remaining alleles, but don't trim:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -sn SAMPLE_1_ACTG \
 *   -env \
 *   -noTrim
 *</pre>
 *
 * <h4>Select a sample and restrict the output vcf to a set of intervals:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -L /path/to/my.interval_list \
 *   -sn SAMPLE_1_ACTG
 * </pre>
 *
 * <h4>Select all calls missed in my vcf, but present in HapMap (useful to take a look at why these variants weren't called in my dataset):</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V hapmap.vcf \
 *   --discordance myCalls.vcf \
 *   -o output.vcf \
 *   -sn mySample
 * </pre>
 *
 * <h4>Select all calls made by both myCalls and theirCalls (useful to take a look at what is consistent between two callers):</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V myCalls.vcf \
 *   --concordance theirCalls.vcf \
 *   -o output.vcf \
 *   -sn mySample
 * </pre>
 *
 * <h4>Generating a VCF of all the variants that are mendelian violations. The optional argument '-mvq' restricts the selection to sites that have a QUAL score of 50 or more</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -ped family.ped \
 *   -mv -mvq 50 \
 *   -o violations.vcf
 * </pre>
 *
 * <h4>Generating a VCF of all the variants that are not mendelian violations. The optional argument '-mvq' together with '-invMv' restricts the selection to sites that have a QUAL score of 50 or less</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -ped family.ped \
 *   -mv -mvq 50 -invMv \
 *   -o violations.vcf
 * </pre>
 *
 * <h4>Create a set with 50% of the total number of variants in the variant VCF:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -fraction 0.5
 * </pre>
 *
 * <h4>Select only indels between 2 and 5 bases long from a VCF:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -selectType INDEL
 *   --minIndelSize 2
 *   --maxIndelSize 5
 * </pre>
 *
 * <h4>Exclude indels from a VCF:</h4>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   --selectTypeToExclude INDEL
 * </pre>
 *
 * <h4>Select only multi-allelic SNPs and MNPs from a VCF (i.e. SNPs with more than one allele listed in the ALT column):</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T SelectVariants \
 *   -R reference.fasta \
 *   -V input.vcf \
 *   -o output.vcf \
 *   -selectType SNP -selectType MNP \
 *   -restrictAllelesTo MULTIALLELIC
 * </pre>
 *
 * <h4>Select IDs in fileKeep and exclude IDs in fileExclude:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   -o output.vcf \
 *   -IDs fileKeep \
 *   -excludeIDs fileExclude
 * </pre>
 *
 * <h4>Select sites where there are between 2 and 5 samples and between 10 and 50 percent of the sample genotypes are filtered:</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   --maxFilteredGenotypes 5
 *   --minFilteredGenotypes 2
 *   --maxFractionFilteredGenotypes 0.60
 *   --minFractionFilteredGenotypes 0.10
 * </pre>
 *
 *  <h4>Set filtered genotypes to no-call (./.):</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T SelectVariants \
 *   --variant input.vcf \
 *   --setFilteredGtToNocall
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
public class SelectVariants extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    static final int MAX_FILTERED_GENOTYPES_DEFAULT_VALUE  = Integer.MAX_VALUE;
    static final int MIN_FILTERED_GENOTYPES_DEFAULT_VALUE  = 0;
    static final double MAX_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE = 1.0;
    static final double MIN_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE = 0.0;
    
    @ArgumentCollection protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * A site is considered discordant if there exists some sample in the variant track that has a non-reference genotype
     * and either the site isn't present in this track, the sample isn't present in this track,
     * or the sample is called reference in this track.
     */
    @Input(fullName="discordance", shortName = "disc", doc="Output variants not called in this comparison track", required=false)
    protected RodBinding<VariantContext> discordanceTrack;

    /**
     * A site is considered concordant if (1) we are not looking for specific samples and there is a variant called
     * in both the variant and concordance tracks or (2) every sample present in the variant track is present in the
     * concordance track and they have the sample genotype call.
     */
    @Input(fullName="concordance", shortName = "conc", doc="Output variants also called in this comparison track", required=false)
    protected RodBinding<VariantContext> concordanceTrack;

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    /**
     * This argument can be specified multiple times in order to provide multiple sample names.
     */
    @Argument(fullName="sample_name", shortName="sn", doc="Include genotypes from this sample", required=false)
    public Set<String> sampleNames = new HashSet<>(0);

    /**
     * Using a regular expression allows you to match multiple sample names that have that pattern in common. This
     * argument can be specified multiple times in order to use multiple different matching patterns.
     */
    @Argument(fullName="sample_expressions", shortName="se", doc="Regular expression to select multiple samples", required=false)
    public Set<String> sampleExpressions;

    /**
     * Sample names should be in a plain text file listing one sample name per line. This argument can be specified multiple times in order to provide
     * multiple sample list files.
     */
    @Input(fullName="sample_file", shortName="sf", doc="File containing a list of samples to include", required=false)
    public Set<File> sampleFiles;

    /**
     * Note that sample exclusion takes precedence over inclusion, so that if a sample is in both lists it will be
     * excluded. This argument can be specified multiple times in order to provide multiple sample names.
     */
    @Argument(fullName="exclude_sample_name", shortName="xl_sn", doc="Exclude genotypes from this sample", required=false)
    public Set<String> XLsampleNames = new HashSet<>(0);

    /**
     * Sample names should be in a plain text file listing one sample name per line. Note that sample exclusion takes precedence over inclusion, so that
     * if a sample is in both lists it will be excluded. This argument can be specified multiple times in order to
     * provide multiple sample list files.
     */
    @Input(fullName="exclude_sample_file", shortName="xl_sf", doc="List of samples to exclude", required=false)
    public Set<File> XLsampleFiles = new HashSet<>(0);

    /**
     * Using a regular expression allows you to match multiple sample names that have that pattern in common. Note that sample exclusion takes precedence
     * over inclusion, so that if a sample is in both lists it will be excluded. This  argument can be specified multiple times in order to use multiple
     * different matching patterns.
     */
    @Input(fullName="exclude_sample_expressions", shortName="xl_se", doc="List of sample expressions to exclude", required=false)
    public Set<String> XLsampleExpressions = new HashSet<>(0);

    /**
     * See example commands above for detailed usage examples. Note that these expressions are evaluated *after* the
     * specified samples are extracted and the INFO field annotations are updated.
     */
    @Argument(shortName="select", doc="One or more criteria to use when selecting the data", required=false)
    public ArrayList<String> selectExpressions = new ArrayList<>();

    /**
     * Invert the selection criteria for -select.
     */
    @Argument(shortName="invertSelect", doc="Invert the selection criteria for -select", required=false)
    protected boolean invertSelect = false;

    /*
     * If this flag is enabled, sites that are found to be non-variant after the subsetting procedure (i.e. where none
     * of the selected samples display evidence of variation) will be excluded from the output.
     */
    @Argument(fullName="excludeNonVariants", shortName="env", doc="Don't include non-variant sites", required=false)
    protected boolean XLnonVariants = false;

    /**
     * If this flag is enabled, sites that have been marked as filtered (i.e. have anything other than `.` or `PASS`
     * in the FILTER field) will be excluded from the output.
     */
    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered sites", required=false)
    protected boolean XLfiltered = false;

    /**
     * The default behavior of this tool is to remove bases common to all remaining alleles after subsetting
     * operations have been completed, leaving only their minimal representation. If this flag is enabled, the original
     * alleles will be preserved as recorded in the input VCF.
     */
    @Argument(fullName="preserveAlleles", shortName="noTrim", doc="Preserve original alleles, do not trim", required=false)
    protected boolean preserveAlleles = false;

    /**
     * When this flag is enabled, all alternate alleles that are not present in the (output) samples will be removed.
     * Note that this even extends to biallelic SNPs - if the alternate allele is not present in any sample, it will be
     * removed and the record will contain a '.' in the ALT column. Note also that sites-only VCFs, by definition, do
     * not include the alternate allele in any genotype calls.
     */
    @Argument(fullName="removeUnusedAlternates", shortName="trimAlternates", doc="Remove alternate alleles not present in any genotypes", required=false)
    protected boolean removeUnusedAlternates = false;

    /**
     * When this argument is used, we can choose to include only multiallelic or biallelic sites, depending on how many alleles are listed in the ALT column of a VCF.
     * For example, a multiallelic record such as:
     *     1    100 .   A   AAA,AAAAA
     * will be excluded if `-restrictAllelesTo BIALLELIC` is used, because there are two alternate alleles, whereas a record such as:
     *     1    100 .   A  T
     * will be included in that case, but would be excluded if `-restrictAllelesTo MULTIALLELIC` is used.
     * Valid options are ALL (default), MULTIALLELIC or BIALLELIC.
     */
    @Argument(fullName="restrictAllelesTo", shortName="restrictAllelesTo", doc="Select only variants of a particular allelicity", required=false)
    private  NumberAlleleRestriction alleleRestriction = NumberAlleleRestriction.ALL;

    /**
     * When subsetting a callset, this tool recalculates the AC, AF, and AN values corresponding to the contents of the
     * subset. If this flag is enabled, the original values of those annotations will be stored in new annotations called
     * AC_Orig, AF_Orig, and AN_Orig.
     */
    @Argument(fullName="keepOriginalAC", shortName="keepOriginalAC", doc="Store the original AC, AF, and AN values after subsetting", required=false)
    private boolean keepOriginalChrCounts = false;

    /**
     * When subsetting a callset, this tool recalculates the site-level (INFO field) DP value corresponding to the contents of the
     * subset. If this flag is enabled, the original value of the DP annotation will be stored in a new annotation called
     * DP_Orig.
     */
    @Argument(fullName="keepOriginalDP", shortName="keepOriginalDP", doc="Store the original DP value after subsetting", required=false)
    private boolean keepOriginalDepth = false;

    /**
     * If this flag is enabled, this tool will select only variants that correspond to a mendelian violation as
     * determined on the basis of family structure. Requires passing a pedigree file using the engine-level
     * `-ped` argument.
     */
    @Argument(fullName="mendelianViolation", shortName="mv", doc="Output mendelian violation sites only", required=false)
    private Boolean mendelianViolations = false;

    /**
     * If this flag is enabled, this tool will select only variants that do not correspond to a mendelian violation as
     * determined on the basis of family structure. Requires passing a pedigree file using the engine-level
     * `-ped` argument.
     */
    @Argument(fullName="invertMendelianViolation", shortName="invMv", doc="Output non-mendelian violation sites only", required=false)
    private Boolean invertMendelianViolations = false;

    /**
     * This argument specifies the genotype quality (GQ) threshold that all members of a trio must have in order
     * for a site to be accepted as a mendelian violation. Note that the `-mv` flag must be set for this argument to have an effect.
     */
    @Argument(fullName="mendelianViolationQualThreshold", shortName="mvq", doc="Minimum GQ score for each trio member to accept a site as a violation", required=false)
    protected double medelianViolationQualThreshold = 0;

    /**
     * The value of this argument should be a number between 0 and 1 specifying the fraction of total variants to be
     * randomly selected from the input callset. Note that this is done using a probabilistic function, so the final
     * result is not guaranteed to carry the exact fraction requested. Can be used for large fractions.
     */
    @Argument(fullName="select_random_fraction", shortName="fraction", doc="Select a fraction of variants at random from the input", required=false)
    protected double fractionRandom = 0;

    /**
     * The value of this argument should be a number between 0 and 1 specifying the fraction of total variants to be
     * randomly selected from the input callset and set to no-call (./). Note that this is done using a probabilistic
     * function, so the final result is not guaranteed to carry the exact fraction requested. Can be used for large fractions.
     */
    @Argument(fullName="remove_fraction_genotypes", shortName="fractionGenotypes", doc="Select a fraction of genotypes at random from the input and sets them to no-call", required=false)
    protected double fractionGenotypes = 0;

    /**
     * This argument selects particular kinds of variants out of a list. If left empty, there is no type selection
     * and all variant types are considered for other selection criteria. Valid types are INDEL, SNP, MIXED, MNP,
     * SYMBOLIC, NO_VARIATION. Can be specified multiple times.
     */
    @Argument(fullName="selectTypeToInclude", shortName="selectType", doc="Select only a certain type of variants from the input file", required=false)
    private List<VariantContext.Type> typesToInclude = new ArrayList<>();

    /**
     * This argument excludes particular kinds of variants out of a list. If left empty, there is no type selection
     * and all variant types are considered for other selection criteria. Valid types are INDEL, SNP, MIXED, MNP,
     * SYMBOLIC, NO_VARIATION. Can be specified multiple times.
     */
    @Argument(fullName="selectTypeToExclude", shortName="xlSelectType", doc="Do not select certain type of variants from the input file", required=false)
    private List<VariantContext.Type> typesToExclude = new ArrayList<>();

    /**
     * If a file containing a list of IDs is provided to this argument, the tool will only select variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. The expected file format
     * is simply plain text with one ID per line.
     */
    @Argument(fullName="keepIDs", shortName="IDs", doc="List of variant IDs to select", required=false)
    private File rsIDFile = null;

    /**
     * If a file containing a list of IDs is provided to this argument, the tool will not select variants whose ID
     * field is present in this list of IDs. The matching is done by exact string matching. The expected file format
     * is simply plain text with one ID per line.
     */
    @Argument(fullName="excludeIDs", shortName="xlIDs", doc="List of variant IDs to select", required=false)
    private File XLrsIDFile = null;

    @Hidden
    @Argument(fullName="fullyDecode", doc="If true, the incoming VariantContext will be fully decoded", required=false)
    private boolean fullyDecode = false;

    @Hidden
    @Argument(fullName="justRead", doc="If true, we won't actually write the output file.  For efficiency testing only", required=false)
    private boolean justRead = false;

    /**
     * If this argument is provided, indels that are larger than the specified size will be excluded.
     */
    @Argument(fullName="maxIndelSize", required=false, doc="Maximum size of indels to include")
    private int maxIndelSize = Integer.MAX_VALUE;

    /**
     * If this argument is provided, indels that are smaller than the specified size will be excluded.
     */
    @Argument(fullName="minIndelSize", required=false, doc="Minimum size of indels to include")
    private int minIndelSize = 0;

    /**
     * If this argument is provided, select sites where at most a maximum number of samples are filtered at the genotype level.
     */
    @Argument(fullName="maxFilteredGenotypes", required=false, doc="Maximum number of samples filtered at the genotype level")
    private int maxFilteredGenotypes = MAX_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where at least a minimum number of samples are filtered at the genotype level.
     */
    @Argument(fullName="minFilteredGenotypes", required=false, doc="Minimum number of samples filtered at the genotype level")
    private int minFilteredGenotypes = MIN_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where a fraction or less of the samples are filtered at the genotype level.
     */
    @Argument(fullName="maxFractionFilteredGenotypes", required=false, doc="Maximum fraction of samples filtered at the genotype level")
    private double maxFractionFilteredGenotypes = MAX_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, select sites where a fraction or more of the samples are filtered at the genotype level.
     */
    @Argument(fullName="minFractionFilteredGenotypes", required=false, doc="Maximum fraction of samples filtered at the genotype level")
    private double minFractionFilteredGenotypes = MIN_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE;

    /**
     * If this argument is provided, set filtered genotypes to no-call (./.).
     */
    @Argument(fullName="setFilteredGtToNocall", required=false, doc="Set filtered genotypes to no-call")
    private boolean setFilteredGenotypesToNocall = false;

    @Hidden
    @Argument(fullName="ALLOW_NONOVERLAPPING_COMMAND_LINE_SAMPLES", required=false, doc="Allow samples other than those in the VCF to be specified on the command line. These samples will be ignored.")
    private boolean allowNonOverlappingCommandLineSamples = false;

    /**
     * If this argument is provided, the output will be compliant with the version in the header, however it will also
     * cause the tool to run slower than without the argument. Without the argument the header will be compliant with
     * the up-to-date version, but the output in the body may not be compliant. If an up-to-date input file is used,
     * then the output will also be up-to-date regardless of this argument.
     */
    @Argument(fullName="forceValidOutput", required=false, doc="Forces output VCF to be compliant to up-to-date version")
    private boolean forceValidOutput = false;

    public enum NumberAlleleRestriction {
        ALL,
        BIALLELIC,
        MULTIALLELIC
    }

    private ArrayList<VariantContext.Type> selectedTypes = new ArrayList<>();
    private ArrayList<String> selectNames = new ArrayList<>();
    private List<VariantContextUtils.JexlVCMatchExp> jexls = null;

    private TreeSet<String> samples = new TreeSet<>();
    private boolean noSamplesSpecified = false;

    private boolean discordanceOnly = false;
    private boolean concordanceOnly = false;

    private MendelianViolation mv;


    /* variables used by the SELECT RANDOM modules */
    private boolean selectRandomFraction = false;

    //Random number generator for the genotypes to remove
    private Random randomGenotypes = new Random();

    private Set<String> IDsToKeep = null;
    private Set<String> IDsToRemove = null;
    private Map<String, VCFHeader> vcfRods;

    private final List<Allele> diploidNoCallAlleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        // Get list of samples to include in the output
        List<String> rodNames = Arrays.asList(variantCollection.variants.getName());

        vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        TreeSet<String> vcfSamples = new TreeSet<>(SampleUtils.getSampleList(vcfRods, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE));

        Collection<String> samplesFromFile = SampleUtils.getSamplesFromFiles(sampleFiles);
        Collection<String> samplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, sampleExpressions);

        // first, check overlap between requested and present samples
        Set<String> commandLineUniqueSamples = new HashSet<>(samplesFromFile.size()+samplesFromExpressions.size()+sampleNames.size());
        commandLineUniqueSamples.addAll(samplesFromFile);
        commandLineUniqueSamples.addAll(samplesFromExpressions);
        commandLineUniqueSamples.addAll(sampleNames);
        commandLineUniqueSamples.removeAll(vcfSamples);

        // second, add the requested samples
        samples.addAll(sampleNames);
        samples.addAll(samplesFromExpressions);
        samples.addAll(samplesFromFile);

        logger.debug(Utils.join(",", commandLineUniqueSamples));

            if (!commandLineUniqueSamples.isEmpty()) {
                if (allowNonOverlappingCommandLineSamples) {
                    logger.warn("Samples present on command line input that are not present in the VCF. These samples will be ignored.");
                    samples.removeAll(commandLineUniqueSamples);
                } else {
                    throw new UserException.BadInput(String.format("%s%n%n%s%n%n%s%n%n%s",
                            "Samples entered on command line (through -sf or -sn) that are not present in the VCF.",
                            "A list of these samples:",
                            Utils.join(",", commandLineUniqueSamples),
                            "To ignore these samples, run with --allowNonOverlappingCommandLineSamples"));
                }
            }

        // if none were requested, we want all of them
        if ( samples.isEmpty() ) {
            samples.addAll(vcfSamples);
            noSamplesSpecified = true;
        }

        // now, exclude any requested samples
        final Collection<String> XLsamplesFromFile = SampleUtils.getSamplesFromFiles(XLsampleFiles);
        final Collection<String> XLsamplesFromExpressions = SampleUtils.matchSamplesExpressions(vcfSamples, XLsampleExpressions);
        samples.removeAll(XLsamplesFromFile);
        samples.removeAll(XLsampleNames);
        samples.removeAll(XLsamplesFromExpressions);
        noSamplesSpecified = noSamplesSpecified && XLsampleNames.isEmpty() && XLsamplesFromFile.isEmpty() &&
                XLsamplesFromExpressions.isEmpty();

        if ( samples.isEmpty() && !noSamplesSpecified )
            throw new UserException("All samples requested to be included were also requested to be excluded.");

        if ( ! noSamplesSpecified )
            for ( String sample : samples )
            logger.info("Including sample '" + sample + "'");

        // if user specified types to include, add these, otherwise, add all possible variant context types to list of vc types to include
        if (typesToInclude.isEmpty()) {
            for (VariantContext.Type t : VariantContext.Type.values())
                selectedTypes.add(t);
        }
        else {
            for (VariantContext.Type t : typesToInclude)
                selectedTypes.add(t);
        }

        // remove specified types
        for (VariantContext.Type t : typesToExclude)
            selectedTypes.remove(t);

        // Initialize VCF header
        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        headerLines.add(new VCFHeaderLine("source", "SelectVariants"));

        if (keepOriginalChrCounts) {
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AC_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AF_KEY));
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_AN_KEY));
        }
        if (keepOriginalDepth)
            headerLines.add(GATKVCFHeaderLines.getInfoLine(GATKVCFConstants.ORIGINAL_DP_KEY));
        headerLines.addAll(Arrays.asList(ChromosomeCountConstants.descriptions));
        headerLines.add(VCFStandardHeaderLines.getInfoLine(VCFConstants.DEPTH_KEY));

        for (int i = 0; i < selectExpressions.size(); i++) {
            // It's not necessary that the user supply select names for the JEXL expressions, since those
            // expressions will only be needed for omitting records.  Make up the select names here.
            selectNames.add(String.format("select-%d", i));
        }

        jexls = VariantContextUtils.initializeMatchExps(selectNames, selectExpressions);

        // Look at the parameters to decide which analysis to perform
        discordanceOnly = discordanceTrack.isBound();
        if (discordanceOnly) logger.info("Selecting only variants discordant with the track: " + discordanceTrack.getName());

        concordanceOnly = concordanceTrack.isBound();
        if (concordanceOnly) logger.info("Selecting only variants concordant with the track: " + concordanceTrack.getName());

        if (mendelianViolations) {
            mv = new MendelianViolation(medelianViolationQualThreshold,false,true);
        }

        selectRandomFraction = fractionRandom > 0;
        if (selectRandomFraction) logger.info("Selecting approximately " + 100.0*fractionRandom + "% of the variants at random from the variant track");

        // Get variant IDs to keep and removed
        IDsToKeep = getIDsFromFile(rsIDFile);

        IDsToRemove = getIDsFromFile(XLrsIDFile);

        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));
    }

    /**
     * Get IDs from a file
     *
     * @param file file containing the IDs
     * @return set of IDs or null if the file is null
     * @throws UserException.CouldNotReadInputFile if could not read the file
     */
    private Set<String> getIDsFromFile(final File file){
        /** load in the IDs file to a hashset for matching */
        if ( file != null ) {
            Set<String> ids = new HashSet<>();
            try {
                for ( final java.lang.String line : new XReadLines(file).readLines() ) {
                    ids.add(line.trim());
                }
                logger.info("Selecting only variants with one of " + ids.size() + " IDs from " + file);
            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotReadInputFile(file, e);
            }
            return ids;
        }

        return null;
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

        if ( vcs == null || vcs.isEmpty()) {
            return 0;
        }

        for (VariantContext vc : vcs) {
            // an option for performance testing only
            if (fullyDecode)
                vc = vc.fullyDecode(vcfRods.get(vc.getSource()), getToolkit().lenientVCFProcessing());

            if (IDsToKeep != null && !IDsToKeep.contains(vc.getID()))
                continue;

            if (IDsToRemove != null && IDsToRemove.contains(vc.getID()))
                continue;

            if (mendelianViolations && Utils.invertLogic(mv.countViolations(this.getSampleDB().getFamilies(samples), vc) == 0, invertMendelianViolations))
                break;

            if (discordanceOnly) {
                Collection<VariantContext> compVCs = tracker.getValues(discordanceTrack, context.getLocation());
                if (!isDiscordant(vc, compVCs))
                    continue;
            }
            if (concordanceOnly) {
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

            if (containsIndelLargerOrSmallerThan(vc, maxIndelSize, minIndelSize))
                continue;

            if ( needNumFilteredGenotypes()) {
                int numFilteredSamples = numFilteredGenotypes(vc);
                double fractionFilteredGenotypes = samples.isEmpty() ? 0.0 : numFilteredSamples / samples.size();
                if (numFilteredSamples > maxFilteredGenotypes || numFilteredSamples < minFilteredGenotypes ||
                        fractionFilteredGenotypes > maxFractionFilteredGenotypes || fractionFilteredGenotypes < minFractionFilteredGenotypes)
                    continue;
            }

            VariantContext sub = subsetRecord(vc, preserveAlleles, removeUnusedAlternates);

            VariantContext filteredGenotypeToNocall = setFilteredGenotypeToNocall(sub, setFilteredGenotypesToNocall);

            // Not excluding non-variants or subsetted polymorphic variants AND including filtered loci or subsetted variant is not filtered
            if ( (!XLnonVariants || filteredGenotypeToNocall.isPolymorphicInSamples()) && (!XLfiltered || !filteredGenotypeToNocall.isFiltered()) ) {

                // Write the subsetted variant if it matches all of the expressions
                boolean failedJexlMatch = false;

                try {
                    for (VariantContextUtils.JexlVCMatchExp jexl : jexls) {
                        if ( Utils.invertLogic(!VariantContextUtils.match(filteredGenotypeToNocall, jexl), invertSelect) ){
                            failedJexlMatch = true;
                            break;
                        }
                    }
                } catch (IllegalArgumentException e) {
                    /*The IAE thrown by htsjdk already includes an informative error message ("Invalid JEXL
                      expression detected...")*/
                    throw new UserException(e.getMessage());
                }
                if ( !failedJexlMatch &&
                        !justRead &&
                        ( !selectRandomFraction || Utils.getRandomGenerator().nextDouble() < fractionRandom ) ) {
                    vcfWriter.add(filteredGenotypeToNocall);
                }
            }
        }

        return 1;
    }

    /*
     * Determines if any of the alternate alleles are greater than the max indel size or less than the min indel size
     *
     * @param vc            the variant context to check
     * @param maxIndelSize  the maximum size of allowed indels
     * @param minIndelSize  the minimum size of allowed indels
     * @return true if the VC contains an indel larger than maxIndelSize or less than the minIndelSize, false otherwise
     */
    protected static boolean containsIndelLargerOrSmallerThan(final VariantContext vc, final int maxIndelSize, final int minIndelSize) {
        final List<Integer> lengths = vc.getIndelLengths();
        if ( lengths == null )
            return false;

        for ( Integer indelLength : lengths ) {
            if ( Math.abs(indelLength) > maxIndelSize || Math.abs(indelLength) < minIndelSize )
                return true;
        }

        return false;
    }

    /**
     * Find the number of filtered samples
     *
     * @param vc the variant rod VariantContext
     * @return number of filtered samples
     */
    private int numFilteredGenotypes(final VariantContext vc){
        if (vc == null)
            return 0;

        int numFiltered = 0;
        // check if we find it in the variant rod
        GenotypesContext genotypes = vc.getGenotypes(samples);
        for (final Genotype g : genotypes)
            if ( g.isFiltered() && !g.getFilters().isEmpty())
                numFiltered++;

        return numFiltered;
    }

    /**
     * Checks if vc has a variant call for (at least one of) the samples.
     *
     * @param vc the variant rod VariantContext. Here, the variant is the dataset you're looking for discordances to (e.g. HapMap)
     * @param compVCs the comparison VariantContext (discordance)
     * @return true VariantContexts are discordant, false otherwise
     */
    private boolean isDiscordant (final VariantContext vc, final Collection<VariantContext> compVCs) {
        if (vc == null)
            return false;

        // if we're not looking at specific samples then the absence of a compVC means discordance
        if (noSamplesSpecified)
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

    /**
     * Checks if the two variants have the same genotypes for the selected samples
     *
     * @param vc the variant rod VariantContext.
     * @param compVCs the comparison VariantContext
     * @return true if VariantContexts are concordant, false otherwise
     */
    private boolean isConcordant (final VariantContext vc, final Collection<VariantContext> compVCs) {
        if (vc == null || compVCs == null || compVCs.isEmpty())
            return false;

        // if we're not looking for specific samples then the fact that we have both VCs is enough to call it concordant.
        if (noSamplesSpecified)
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

    private boolean sampleHasVariant(final Genotype g) {
        return (g !=null && !g.isHomRef() && (g.isCalled() || (g.isFiltered() && !XLfiltered)));
    }

    private boolean haveSameGenotypes(final Genotype g1, final Genotype g2) {
        if ( g1 == null || g2 == null )
            return false;

        if ((g1.isCalled() && g2.isFiltered()) ||
                (g2.isCalled() && g1.isFiltered()) ||
                (g1.isFiltered() && g2.isFiltered() && XLfiltered))
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
     * @param preserveAlleles should we trim constant sequence from the beginning and/or end of all alleles, or preserve it?
     * @param removeUnusedAlternates removes alternate alleles with AC=0
     * @return the subsetted VariantContext
     */
    private VariantContext subsetRecord(final VariantContext vc, final boolean preserveAlleles, final boolean removeUnusedAlternates) {
        //subContextFromSamples() always decodes the vc, which is a fairly expensive operation.  Avoid if possible
        if ( noSamplesSpecified && !removeUnusedAlternates && !forceValidOutput )
            return vc;

        // strip out the alternate alleles that aren't being used
        final VariantContext sub = vc.subContextFromSamples(samples, removeUnusedAlternates);

        //If no subsetting happened, exit now
        if ( sub.getNSamples() == vc.getNSamples() && sub.getNAlleles() == vc.getNAlleles() )
            return vc;

        final VariantContextBuilder builder = new VariantContextBuilder(sub);

        // if there are fewer alternate alleles now in the selected VC, we need to fix the PL and AD values
        GenotypesContext newGC = GATKVariantContextUtils.updatePLsSACsAD(sub, vc);

        // since the VC has been subset (either by sample or allele), we need to strip out the MLE tags
        builder.rmAttribute(GATKVCFConstants.MLE_ALLELE_COUNT_KEY);
        builder.rmAttribute(GATKVCFConstants.MLE_ALLELE_FREQUENCY_KEY);

        // Remove a fraction of the genotypes if needed
        if ( fractionGenotypes > 0 ){
            final ArrayList<Genotype> genotypes = new ArrayList<>();
            for ( Genotype genotype : newGC ) {
                //Set genotype to no call if it falls in the fraction.
                if(fractionGenotypes>0 && randomGenotypes.nextDouble()<fractionGenotypes){
                    genotypes.add(new GenotypeBuilder(genotype).alleles(diploidNoCallAlleles).noGQ().make());
                }
                else{
                    genotypes.add(genotype);
                }
            }
            newGC = GenotypesContext.create(genotypes);
        }

        builder.genotypes(newGC);

        addAnnotations(builder, vc, sub.getSampleNames());

        final VariantContext subset = builder.make();

        final VariantContext trimmed = preserveAlleles? subset : GATKVariantContextUtils.trimAlleles(subset,true,true);

        return trimmed;
    }

    /**
     * If --setFilteredGtToNocall, set filtered genotypes to no-call
     *
     * @param vc the VariantContext record to set filtered genotypes to no-call
     * @param filteredGenotypesToNocall  set filtered genotypes to non-call?
     * @return the VariantContext with no-call genotypes if the sample was filtered
     */
    private VariantContext setFilteredGenotypeToNocall(final VariantContext vc, final boolean filteredGenotypesToNocall) {

        if ( !filteredGenotypesToNocall )
            return vc;

        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        final GenotypesContext genotypes = GenotypesContext.create(vc.getGenotypes().size());

        for ( final Genotype g : vc.getGenotypes() ) {
            if ( g.isCalled() && g.isFiltered() )
                genotypes.add(new GenotypeBuilder(g).alleles(diploidNoCallAlleles).make());
            else
                genotypes.add(g);
        }

        return builder.genotypes(genotypes).make();
    }
    /*
     * Add annotations to the new VC
     *
     * @param builder     the new VC to annotate
     * @param originalVC  the original VC
     * @param selectedSampleNames the post-selection list of sample names
     */
    private void addAnnotations(final VariantContextBuilder builder, final VariantContext originalVC, final Set<String> selectedSampleNames) {
        if ( fullyDecode ) return; // TODO -- annotations are broken with fully decoded data

        if ( keepOriginalChrCounts ) {
            final int[] indexOfOriginalAlleleForNewAllele;
            final List<Allele> newAlleles = builder.getAlleles();
            final int numOriginalAlleles = originalVC.getNAlleles();

            // if the alleles already match up, we can just copy the previous list of counts
            if ( numOriginalAlleles == newAlleles.size() ) {
                indexOfOriginalAlleleForNewAllele = null;
            }
            // otherwise we need to parse them and select out the correct ones
            else {
                indexOfOriginalAlleleForNewAllele = new int[newAlleles.size() - 1];
                Arrays.fill(indexOfOriginalAlleleForNewAllele, -1);

                // note that we don't care about the reference allele at position 0
                for ( int newI = 1; newI < newAlleles.size(); newI++ ) {
                    final Allele newAlt = newAlleles.get(newI);
                    for ( int oldI = 0; oldI < numOriginalAlleles - 1; oldI++ ) {
                        if ( newAlt.equals(originalVC.getAlternateAllele(oldI), false) ) {
                            indexOfOriginalAlleleForNewAllele[newI - 1] = oldI;
                            break;
                        }
                    }
                }
            }

            if ( originalVC.hasAttribute(VCFConstants.ALLELE_COUNT_KEY) )
                builder.attribute(GATKVCFConstants.ORIGINAL_AC_KEY, getReorderedAttributes(originalVC.getAttribute(VCFConstants.ALLELE_COUNT_KEY), indexOfOriginalAlleleForNewAllele));
            if ( originalVC.hasAttribute(VCFConstants.ALLELE_FREQUENCY_KEY) )
                builder.attribute(GATKVCFConstants.ORIGINAL_AF_KEY, getReorderedAttributes(originalVC.getAttribute(VCFConstants.ALLELE_FREQUENCY_KEY), indexOfOriginalAlleleForNewAllele));
            if ( originalVC.hasAttribute(VCFConstants.ALLELE_NUMBER_KEY) )
                builder.attribute(GATKVCFConstants.ORIGINAL_AN_KEY, originalVC.getAttribute(VCFConstants.ALLELE_NUMBER_KEY));
        }

        VariantContextUtils.calculateChromosomeCounts(builder, false);

        if (keepOriginalDepth && originalVC.hasAttribute(VCFConstants.DEPTH_KEY))
            builder.attribute(GATKVCFConstants.ORIGINAL_DP_KEY, originalVC.getAttribute(VCFConstants.DEPTH_KEY));

        boolean sawDP = false;
        int depth = 0;
        for ( final String sample : selectedSampleNames ) {
            Genotype g = originalVC.getGenotype(sample);

            if ( ! g.isFiltered() ) {
                if ( g.hasDP() ) {
                    depth += g.getDP();
                    sawDP = true;
                }
            }
        }

        if ( sawDP )
            builder.attribute(VCFConstants.DEPTH_KEY, depth);
    }

    /**
     * Pulls out the appropriate tokens from the old ordering of an attribute to the new ordering
     *
     * @param attribute               the non-null attribute (from the INFO field)
     * @param oldToNewIndexOrdering   the mapping from new to old ordering
     * @return non-null Object attribute
     */
    private Object getReorderedAttributes(final Object attribute, final int[] oldToNewIndexOrdering) {
        // if the ordering is the same, then just use the original attribute
        if ( oldToNewIndexOrdering == null )
            return attribute;

        // break the original attributes into separate tokens; unfortunately, this means being smart about class types
        final Object[] tokens;
        if ( attribute.getClass().isArray() )
            tokens = (Object[])attribute;
        else if ( List.class.isAssignableFrom(attribute.getClass()) )
            tokens = ((List)attribute).toArray();
        else
            tokens = attribute.toString().split(VCFConstants.INFO_FIELD_ARRAY_SEPARATOR);

        final List<Object> result = new ArrayList<>();
        for ( final int index : oldToNewIndexOrdering ) {
            if ( index >= tokens.length )
                throw new IllegalArgumentException("the old attribute has an incorrect number of elements: " + attribute);
            result.add(tokens[index]);
        }
        return result;
    }

    /**
     * Need the number of filtered genotypes samples?
     *
     * @return true if any of the filtered genotype samples arguments is used (not the default value), false otherwise
     */
    private boolean needNumFilteredGenotypes(){
        return maxFilteredGenotypes != MAX_FILTERED_GENOTYPES_DEFAULT_VALUE ||
                minFilteredGenotypes != MIN_FILTERED_GENOTYPES_DEFAULT_VALUE ||
                maxFractionFilteredGenotypes != MAX_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE ||
                minFractionFilteredGenotypes != MIN_FRACTION_FILTERED_GENOTYPES_DEFAULT_VALUE;
    }
}
/*
* Copyright (c) 2012 The Broad Institute
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
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.engine.io.stubs.VariantContextWriterStub;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.variant.ChromosomeCountConstants;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;
import org.broadinstitute.gatk.utils.variant.GATKVariantContextUtils;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextUtils;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.util.*;


/**
 * Combine variant records from different sources
 *
 * <p>CombineVariants reads in variants records from separate ROD (Reference-Ordered Data) sources and combines them into
 * a single VCF. Any (unique) name can be used to bind your ROD and any number of sources can be input. This tool aims
 * to fulfill two main possible use cases, reflected by the two combination options (MERGE and UNION), for merging
 * records at the variant level (the first 8 fields of the VCF) or at the genotype level (the rest).</p>
 *
 * <ul>
 * <li><b>MERGE:</b> combines multiple variant records present at the same site in the different input sources into a
 * single variant record in the output. If sample names overlap, then they are "uniquified" by default, which means a
 * suffix is appended to make them unique. <em>Note that in version 3.3, the automatic uniquifying was disabled
 * (unintentionally), and required setting `-genotypeMergeOptions UNIQUIFY` manually.</em></li>
 *
 * <li><b>UNION:</b> assumes that each ROD source represents the same set of samples (although this is not enforced).
 * It uses the priority list (if provided) to emit a single record instance at every position represented in the input RODs.</li>
 * </ul>
 *
 * <p>CombineVariants will emit a record for every site that was present in any of your input VCF files, and will annotate
 * (in the set attribute in the INFO field) whether the record had a PASS or FILTER status in each input ROD . In effect,
 * CombineVariants always produces a union of the input VCFs.  However, any part of the Venn of the merged VCFs
 * can be extracted using JEXL expressions on the set attribute using SelectVariants.  If you want to extract just
 * the records in common between two VCFs, you would first run CombineVariants on the two files to generate a single
 * VCF and then run SelectVariants to extract the common records with `-select 'set == "Intersection"'`, as worked out
 * in the detailed example in the documentation guide.</p>
 *
 * <h3>Input</h3>
 * <p>
 * Two or more variant sets to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF.
 * </p>
 *
 * <h3>Usage examples</h3>
 * &nbsp;
 * <h4>Merge two separate callsets</h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CombineVariants \
 *   -R reference.fasta \
 *   --variant input1.vcf \
 *   --variant input2.vcf \
 *   -o output.vcf \
 *   -genotypeMergeOptions UNIQUIFY
 * </pre>
 *
 * <h4>Get the union of calls made on the same samples </h4>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CombineVariants \
 *   -R reference.fasta \
 *   --variant:foo input1.vcf \
 *   --variant:bar input2.vcf \
 *   -o output.vcf \
 *   -genotypeMergeOptions PRIORITIZE \
 *   -priority foo,bar
 * </pre>
 *
 * <h3>Caveats</h3>
 * <ul>
 * <li>This tool is not intended to manipulate GVCFS! To combine GVCF files output by HaplotypeCaller, use CombineGVCFs.</li>
 * <li>To join intermediate VCFs produced by running jobs in parallel by interval (e.g. by chromosome), use CatVariants.</li>
 * </ul>
 *
 * <h3>Additional notes</h3>
 * <ul>
 * <li> Using this tool's multi-threaded parallelism capability is particularly useful
 * when converting from VCF to BCF2, which can be time-consuming. In this case each thread spends CPU time
 * doing the conversion, and the GATK engine is smart enough to merge the partial BCF2 blocks together
 * efficiently.  However, since this merge runs in only one thread, you can quickly reach diminishing
 * returns with the number of parallel threads.  In our hands, `-nt 4` works well but `-nt 8` tends to be be too much.</li>
 * <li>Since GATK 2.1, when merging multiple VCF records at a site, the combined VCF record has the QUAL of the first
 * VCF record with a non-MISSING QUAL value.  The previous behavior was to take the max QUAL, which could result
 * in strange downstream confusion</li>
 * </ul>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-50,stop=50))
public class CombineVariants extends RodWalker<Integer, Integer> implements TreeReducible<Integer> {
    /**
     * The VCF files to merge together
     *
     * variants can take any number of arguments on the command line.  Each -V argument
     * will be included in the final merged output VCF.  If no explicit name is provided,
     * the -V arguments will be named using the default algorithm: variants, variants2, variants3, etc.
     * The user can override this by providing an explicit name -V:name,vcf for each -V argument,
     * and each named argument will be labeled as such in the output (i.e., set=name rather than
     * set=variants2).  The order of arguments does not matter unless except for the naming, so
     * if you provide an rod priority list and no explicit names than variants, variants2, etc
     * are technically order dependent.  It is strongly recommended to provide explicit names when
     * a rod priority list is provided.
     */
    @Input(fullName="variant", shortName = "V", doc="Input VCF file", required=true)
    public List<RodBindingCollection<VariantContext>> variantCollections;
    final private List<RodBinding<VariantContext>> variants = new ArrayList<>();

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter vcfWriter = null;

    @Argument(shortName="genotypeMergeOptions", doc="Determines how we should merge genotype records for samples shared across the ROD files", required=false)
    public GATKVariantContextUtils.GenotypeMergeType genotypeMergeOption = null;

    @Argument(shortName="filteredRecordsMergeType", doc="Determines how we should handle records seen at the same site in the VCF, but with different FILTER fields", required=false)
    public GATKVariantContextUtils.FilteredRecordMergeType filteredRecordsMergeType = GATKVariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED;

    @Hidden
    @Argument(shortName="multipleAllelesMergeType", doc="Determines how we should handle records seen at the same site in the VCF, but with different allele types (for example, SNP vs. indel)", required=false)
    public GATKVariantContextUtils.MultipleAllelesMergeType multipleAllelesMergeType = GATKVariantContextUtils.MultipleAllelesMergeType.BY_TYPE;

    /**
     * Used when taking the union of variants that contain genotypes.  A complete priority list MUST be provided.
     */
    @Argument(fullName="rod_priority_list", shortName="priority", doc="A comma-separated string describing the priority ordering for the genotypes as far as which record gets emitted", required=false)
    public String PRIORITY_STRING = null;

    @Argument(fullName="printComplexMerges", shortName="printComplexMerges", doc="Print out interesting sites requiring complex compatibility merging", required=false)
    public boolean printComplexMerges = false;

    @Argument(fullName="filteredAreUncalled", shortName="filteredAreUncalled", doc="If true, then filtered VCFs are treated as uncalled, so that filtered set annotations don't appear in the combined VCF", required=false)
    public boolean filteredAreUncalled = false;

    /**
     * Used to generate a sites-only file.
     */
    @Argument(fullName="minimalVCF", shortName="minimalVCF", doc="If true, then the output VCF will contain no INFO or genotype FORMAT fields", required=false)
    public boolean minimalVCF = false;

    @Argument(fullName="excludeNonVariants", shortName="env", doc="Don't include loci found to be non-variant after the combining procedure", required=false)
    public boolean EXCLUDE_NON_VARIANTS = false;

    /**
     * Set to 'null' if you don't want the set field emitted.
     */
    @Argument(fullName="setKey", shortName="setKey", doc="Key used in the INFO key=value tag emitted describing which set the combined VCF record came from", required=false)
    public String SET_KEY = "set";

    /**
     * This option allows the user to perform a simple merge (concatenation) to combine the VCFs, drastically reducing the runtime.
     */
    @Argument(fullName="assumeIdenticalSamples", shortName="assumeIdenticalSamples", doc="If true, assume input VCFs have identical sample sets and disjoint calls", required=false)
    public boolean ASSUME_IDENTICAL_SAMPLES = false;

    @Argument(fullName="minimumN", shortName="minN", doc="Combine variants and output site only if the variant is present in at least N input files.", required=false)
    public int minimumN = 1;

    /**
     * This option allows the suppression of the command line in the VCF header. This is most often usefully when combining variants for dozens or hundreds of smaller VCFs.
     */
    @Argument(fullName="suppressCommandLineHeader", shortName="suppressCommandLineHeader", doc="If true, do not output the header containing the command line used", required=false)
    public boolean SUPPRESS_COMMAND_LINE_HEADER = false;

    @Argument(fullName="mergeInfoWithMaxAC", shortName="mergeInfoWithMaxAC", doc="If true, when VCF records overlap the info field is taken from the one with the max AC instead of only taking the fields which are identical across the overlapping records.", required=false)
    public boolean MERGE_INFO_WITH_MAX_AC = false;

    private List<String> priority = null;

    /** Optimization to strip out genotypes before merging if we are doing a sites_only output */
    private boolean sitesOnlyVCF = false;
    private Set<String> samples;

    public void initialize() {
        Map<String, VCFHeader> vcfRods = GATKVCFUtils.getVCFHeadersFromRods(getToolkit());

        if ( vcfWriter instanceof VariantContextWriterStub) {
            sitesOnlyVCF = ((VariantContextWriterStub)vcfWriter).getWriterOptions().contains(Options.DO_NOT_WRITE_GENOTYPES);
            if ( sitesOnlyVCF ) logger.info("Pre-stripping genotypes for performance");
        } else
            logger.warn("VCF output file not an instance of VCFWriterStub; cannot enable sites only output option");

        validateAnnotateUnionArguments();

        final boolean sampleNamesAreUnique = SampleUtils.verifyUniqueSamplesNames(vcfRods);

        if (genotypeMergeOption == null && !ASSUME_IDENTICAL_SAMPLES) {
            if (!sampleNamesAreUnique)
                throw new UserException("Duplicate sample names were discovered but no genotypemergeoption was supplied. " +
                    "To combine samples without merging specify --genotypemergeoption UNIQUIFY. Merging duplicate samples " +
                    "without specified priority is unsupported, but can be achieved by specifying --genotypemergeoption UNSORTED.");
            else
                genotypeMergeOption = GATKVariantContextUtils.GenotypeMergeType.UNSORTED;
        }

        if ( PRIORITY_STRING == null && genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE) {
            //PRIORITY_STRING = Utils.join(",", vcfRods.keySet());  Deleted by Ami (7/10/12)
            logger.info("Priority string is not provided, using arbitrary genotyping order: "+priority);
        }

        if (genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE &&
                !sampleNamesAreUnique)
            throw new IllegalStateException("REQUIRE_UNIQUE sample names is true but duplicate names were discovered.");

        samples = sitesOnlyVCF ? Collections.<String>emptySet() : SampleUtils.getSampleList(vcfRods, genotypeMergeOption);

        if ( SET_KEY.toLowerCase().equals("null") )
            SET_KEY = null;

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), true);
        if ( SET_KEY != null )
            headerLines.add(new VCFInfoHeaderLine(SET_KEY, 1, VCFHeaderLineType.String, "Source VCF for the merged record in CombineVariants"));
        if ( !ASSUME_IDENTICAL_SAMPLES )
             headerLines.addAll(Arrays.asList(ChromosomeCountConstants.descriptions));
        VCFHeader vcfHeader = new VCFHeader(headerLines, samples);
        vcfHeader.setWriteCommandLine(!SUPPRESS_COMMAND_LINE_HEADER);
        vcfWriter.writeHeader(vcfHeader);

        // collect the actual rod bindings into a list for use later
        for ( final RodBindingCollection<VariantContext> variantCollection : variantCollections )
            variants.addAll(variantCollection.getRodBindings());
    }

    private void validateAnnotateUnionArguments() {
        Set<String> rodNames = SampleUtils.getRodNamesWithVCFHeader(getToolkit(), null);

        if ( genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.PRIORITIZE && PRIORITY_STRING == null )
            throw new UserException.MissingArgument("rod_priority_list", "Priority string must be provided if you want to prioritize genotypes");

        if ( PRIORITY_STRING != null){
            priority = new ArrayList<>(Arrays.asList(PRIORITY_STRING.split(",")));
            if ( rodNames.size() != priority.size() )
                throw new UserException.BadArgumentValue("rod_priority_list", "The priority list must contain exactly one rod binding per ROD provided to the GATK: rodNames=" + rodNames + " priority=" + priority);

            if ( ! rodNames.containsAll(priority) )
                throw new UserException.BadArgumentValue("rod_priority_list", "Not all priority elements provided as input RODs: " + PRIORITY_STRING);
        }

    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        final Set<String> rodNames = SampleUtils.getRodNamesWithVCFHeader(getToolkit(), null);
        // get all of the vcf rods at this locus
        // Need to provide reference bases to simpleMerge starting at current locus
        Collection<VariantContext> vcs = tracker.getValues(variants, context.getLocation());

        if ( sitesOnlyVCF ) {
            vcs = VariantContextUtils.sitesOnlyVariantContexts(vcs);
        }

        if ( ASSUME_IDENTICAL_SAMPLES ) {
            for ( final VariantContext vc : vcs ) {
                vcfWriter.add(vc);
            }

            return vcs.isEmpty() ? 0 : 1;
        }

        int numFilteredRecords = 0;
        for (final VariantContext vc : vcs) {
            if (vc.filtersWereApplied() && vc.isFiltered())
                numFilteredRecords++;
        }

        if (minimumN > 1 && (vcs.size() - numFilteredRecords < minimumN))
            return 0;

        final List<VariantContext> mergedVCs = new ArrayList<>();

        if (multipleAllelesMergeType == GATKVariantContextUtils.MultipleAllelesMergeType.BY_TYPE) {
            final Map<VariantContext.Type, List<VariantContext>> VCsByType = GATKVariantContextUtils.separateVariantContextsByType(vcs);

            // TODO -- clean this up in a refactoring
            // merge NO_VARIATION into another type of variant (based on the ordering in VariantContext.Type)
            if ( VCsByType.containsKey(VariantContext.Type.NO_VARIATION) && VCsByType.size() > 1 ) {
                final List<VariantContext> refs = VCsByType.remove(VariantContext.Type.NO_VARIATION);
                for ( final VariantContext.Type type : VariantContext.Type.values() ) {
                    if ( VCsByType.containsKey(type) ) {
                        VCsByType.get(type).addAll(refs);
                        break;
                    }
                }
            }

            // iterate over the types so that it's deterministic
            for (final VariantContext.Type type : VariantContext.Type.values()) {
                // make sure that it is a variant or in case it is not, that we want to include the sites with no variants
                if (!EXCLUDE_NON_VARIANTS || !type.equals(VariantContext.Type.NO_VARIATION)) {
                    if (VCsByType.containsKey(type)) {
                        mergedVCs.add(GATKVariantContextUtils.simpleMerge(VCsByType.get(type), priority, rodNames.size(),
                                filteredRecordsMergeType, genotypeMergeOption, true, printComplexMerges,
                                SET_KEY, filteredAreUncalled, MERGE_INFO_WITH_MAX_AC));
                    }
                }
            }
        }
        else if (multipleAllelesMergeType == GATKVariantContextUtils.MultipleAllelesMergeType.MIX_TYPES) {
            mergedVCs.add(GATKVariantContextUtils.simpleMerge(vcs, priority, rodNames.size(), filteredRecordsMergeType,
                    genotypeMergeOption, true, printComplexMerges, SET_KEY, filteredAreUncalled, MERGE_INFO_WITH_MAX_AC));
        }
        else {
            logger.warn("Ignoring all records at site " + ref.getLocus());
        }

        for ( final VariantContext mergedVC : mergedVCs ) {
            // only operate at the start of events
            if ( mergedVC == null )
                continue;

            if ( mergedVC.hasAllele(GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) )
                throw new UserException("CombineVariants should not be used to merge gVCFs produced by the HaplotypeCaller; use CombineGVCFs instead");

            final VariantContextBuilder builder = new VariantContextBuilder(mergedVC);
            // re-compute chromosome counts
            VariantContextUtils.calculateChromosomeCounts(builder, false);

            if ( minimalVCF )
                GATKVariantContextUtils.pruneVariantContext(builder, Arrays.asList(SET_KEY));
            final VariantContext vc = builder.make();
            if( !EXCLUDE_NON_VARIANTS || vc.isPolymorphicInSamples() )
                vcfWriter.add(builder.make());
        }

        return vcs.isEmpty() ? 0 : 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    @Override
    public Integer treeReduce(Integer lhs, Integer rhs) {
        return reduce(lhs, rhs);
    }

    public void onTraversalDone(Integer sum) {}
}

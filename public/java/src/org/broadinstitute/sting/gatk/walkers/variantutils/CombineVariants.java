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

package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.stubs.VariantContextWriterStub;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.TreeReducible;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.gatk.walkers.annotator.ChromosomeCountConstants;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.help.HelpConstants;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.sting.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.VariantContext;
import org.broadinstitute.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.variant.variantcontext.VariantContextUtils;
import org.broadinstitute.variant.variantcontext.writer.Options;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.util.*;

/**
 * Combines VCF records from different sources.
 *
 * <p>
 * CombineVariants combines VCF records from different sources. Any (unique) name can be used to bind your rod data
 * and any number of sources can be input. This tool currently supports two different combination types for each of
 * variants (the first 8 fields of the VCF) and genotypes (the rest).
 * Merge: combines multiple records into a single one; if sample names overlap then they are uniquified.
 * Union: assumes each rod represents the same set of samples (although this is not enforced); using the
 * priority list (if provided), it emits a single record instance at every position represented in the rods.
 *
 * CombineVariants will include a record at every site in all of your input VCF files, and annotate which input ROD
 * bindings the record is present, pass, or filtered in in the set attribute in the INFO field. In effect,
 * CombineVariants always produces a union of the input VCFs.  However, any part of the Venn of the N merged VCFs
 * can be exacted using JEXL expressions on the set attribute using SelectVariants.  If you want to extract just
 * the records in common between two VCFs, you would first run CombineVariants on the two files to generate a single
 * VCF and then run SelectVariants to extract the common records with -select 'set == "Intersection"', as worked out
 * in the detailed example in the documentation guide.
 *
 * Note that CombineVariants supports multi-threaded parallelism (8/15/12).  This is particularly useful
 * when converting from VCF to BCF2, which can be expensive.  In this case each thread spends CPU time
 * doing the conversion, and the GATK engine is smart enough to merge the partial BCF2 blocks together
 * efficiency.  However, since this merge runs in only one thread, you can quickly reach diminishing
 * returns with the number of parallel threads.  -nt 4 works well but -nt 8 may be too much.
 *
 * Some fine details about the merging algorithm:
 *   <ul>
 *   <li> As of GATK 2.1, when merging multiple VCF records at a site, the combined VCF record has the QUAL of
 *      the first VCF record with a non-MISSING QUAL value.  The previous behavior was to take the
 *      max QUAL, which resulted in sometime strange downstream confusion</li>
 *   </ul>
 *
 * <h3>Input</h3>
 * <p>
 * One or more variant sets to combine.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A combined VCF.
 * </p>
 *
 * <h3>Examples</h3>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CombineVariants \
 *   --variant input1.vcf \
 *   --variant input2.vcf \
 *   -o output.vcf \
 *   -genotypeMergeOptions UNIQUIFY
 *
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T CombineVariants \
 *   --variant:foo input1.vcf \
 *   --variant:bar input2.vcf \
 *   -o output.vcf \
 *   -genotypeMergeOptions PRIORITIZE
 *   -priority foo,bar
 * </pre>
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
    public List<RodBinding<VariantContext>> variants;

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

    @Argument(fullName="combineAnnotations", shortName="combineAnnotations", doc="If true, combine the annotation values in some straightforward manner assuming the input callsets are i.i.d.", required=false)
    public boolean COMBINE_ANNOTATIONS = false;

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
        if ( PRIORITY_STRING == null && genotypeMergeOption == null) {
            genotypeMergeOption = GATKVariantContextUtils.GenotypeMergeType.UNSORTED;
            //PRIORITY_STRING = Utils.join(",", vcfRods.keySet());  Deleted by Ami (7/10/12)
            logger.info("Priority string is not provided, using arbitrary genotyping order: "+priority);
        }

        if (genotypeMergeOption == GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE &&
                !SampleUtils.verifyUniqueSamplesNames(vcfRods))
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
        Collection<VariantContext> potentialRefVCs = tracker.getValues(variants);
        potentialRefVCs.removeAll(vcs);

        if ( sitesOnlyVCF ) {
            vcs = VariantContextUtils.sitesOnlyVariantContexts(vcs);
            potentialRefVCs = VariantContextUtils.sitesOnlyVariantContexts(potentialRefVCs);
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
                        mergedVCs.add(GATKVariantContextUtils.simpleMerge(VCsByType.get(type), potentialRefVCs,
                                priority, rodNames.size(), filteredRecordsMergeType, genotypeMergeOption, true, printComplexMerges,
                                SET_KEY, filteredAreUncalled, MERGE_INFO_WITH_MAX_AC, COMBINE_ANNOTATIONS));
                    }
                }
            }
        }
        else if (multipleAllelesMergeType == GATKVariantContextUtils.MultipleAllelesMergeType.MIX_TYPES) {
            mergedVCs.add(GATKVariantContextUtils.simpleMerge(vcs, potentialRefVCs,
                    priority, rodNames.size(), filteredRecordsMergeType, genotypeMergeOption, true, printComplexMerges,
                    SET_KEY, filteredAreUncalled, MERGE_INFO_WITH_MAX_AC, COMBINE_ANNOTATIONS));
        }
        else {
            logger.warn("Ignoring all records at site " + ref.getLocus());
        }

        for ( final VariantContext mergedVC : mergedVCs ) {
            // only operate at the start of events
            if ( mergedVC == null )
                continue;

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

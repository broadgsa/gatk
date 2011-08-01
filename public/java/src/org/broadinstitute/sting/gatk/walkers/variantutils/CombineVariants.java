/*
 * Copyright (c) 2010 The Broad Institute
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

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.io.stubs.VCFWriterStub;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Combines VCF records from different sources; supports both full merges and set unions.
 * Merge: combines multiple records into a single one; if sample names overlap then they are uniquified.
 * Union: assumes each rod represents the same set of samples (although this is not enforced); using the
 *   priority list (if provided), emits a single record instance at every position represented in the rods.
 */
@Reference(window=@Window(start=-50,stop=50))
@Requires(value={})
public class CombineVariants extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfWriter = null;

    // the types of combinations we currently allow
    @Argument(shortName="genotypeMergeOptions", doc="How should we merge genotype records for samples shared across the ROD files?", required=false)
    public VariantContextUtils.GenotypeMergeType genotypeMergeOption = VariantContextUtils.GenotypeMergeType.PRIORITIZE;

    @Argument(shortName="filteredRecordsMergeType", doc="How should we deal with records seen at the same site in the VCF, but with different FILTER fields?  KEEP_IF_ANY_UNFILTERED PASSes the record if any record is unfiltered, KEEP_IF_ALL_UNFILTERED requires all records to be unfiltered", required=false)
    public VariantContextUtils.FilteredRecordMergeType filteredRecordsMergeType = VariantContextUtils.FilteredRecordMergeType.KEEP_IF_ANY_UNFILTERED;

    @Argument(fullName="rod_priority_list", shortName="priority", doc="When taking the union of variants containing genotypes: a comma-separated string describing the priority ordering for the genotypes as far as which record gets emitted; a complete priority list MUST be provided", required=false)
    public String PRIORITY_STRING = null;

    @Argument(fullName="printComplexMerges", shortName="printComplexMerges", doc="Print out interesting sites requiring complex compatibility merging", required=false)
    public boolean printComplexMerges = false;

    @Argument(fullName="filteredAreUncalled", shortName="filteredAreUncalled", doc="If true, then filtered VCFs are treated as uncalled, so that filtered set annotation don't appear in the combined VCF", required=false)
    public boolean filteredAreUncalled = false;

    @Argument(fullName="minimalVCF", shortName="minimalVCF", doc="If true, then the output VCF will contain no INFO or genotype INFO field", required=false)
    public boolean minimalVCF = false;

    @Argument(fullName="setKey", shortName="setKey", doc="Key, by default set, in the INFO key=value tag emitted describing which set the combined VCF record came from.  Set to null if you don't want the set field emitted.", required=false)
    public String SET_KEY = "set";

    @Argument(fullName="assumeIdenticalSamples", shortName="assumeIdenticalSamples", doc="If true, assume input VCFs have identical sample sets and disjoint calls so that one can simply perform a merge sort to combine the VCFs into one, drastically reducing the runtime.", required=false)
    public boolean ASSUME_IDENTICAL_SAMPLES = false;

    @Argument(fullName="minimumN", shortName="minN", doc="Combine variants and output site only if variant is present in at least N input files.", required=false)
    public int minimumN = 1;

    @Hidden
    @Argument(fullName="masterMerge", shortName="master", doc="Master merge mode -- experts only.  You need to look at the code to understand it", required=false)
    public boolean master = false;

    @Hidden
    @Argument(fullName="mergeInfoWithMaxAC", shortName="mergeInfoWithMaxAC", doc="If true, when VCF records overlap the info field is taken from the one with the max AC instead of only taking the fields which are identical across the overlapping records.", required=false)
    public boolean MERGE_INFO_WITH_MAX_AC = false;

    private List<String> priority = null;

    /** Optimization to strip out genotypes before merging if we are doing a sites_only output */
    private boolean sitesOnlyVCF = false;

    public void initialize() {
        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), null);

        if ( PRIORITY_STRING == null ) {
            PRIORITY_STRING = Utils.join(",", vcfRods.keySet());
            logger.info("Priority string not provided, using arbitrary genotyping order: " + PRIORITY_STRING);
        }

        validateAnnotateUnionArguments();
        Set<String> samples = SampleUtils.getSampleList(vcfRods, genotypeMergeOption);

        if ( SET_KEY.toLowerCase().equals("null") )
            SET_KEY = null;

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        if ( SET_KEY != null )
            headerLines.add(new VCFInfoHeaderLine(SET_KEY, 1, VCFHeaderLineType.String, "Source VCF for the merged record in CombineVariants"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, sitesOnlyVCF ? Collections.<String>emptySet() : samples));

        if ( vcfWriter instanceof VCFWriterStub) {
            sitesOnlyVCF = ((VCFWriterStub)vcfWriter).doNotWriteGenotypes();
            if ( sitesOnlyVCF ) logger.info("Pre-stripping genotypes for performance");
        } else
            logger.warn("VCF output file not an instance of VCFWriterStub; cannot enable sites only output option");
    }

    private void validateAnnotateUnionArguments() {
        Set<String> rodNames = SampleUtils.getRodNamesWithVCFHeader(getToolkit(), null);

        if ( genotypeMergeOption == VariantContextUtils.GenotypeMergeType.PRIORITIZE && PRIORITY_STRING == null )
            throw new UserException.MissingArgument("rod_priority_list", "Priority string must be provided if you want to prioritize genotypes");

        if ( genotypeMergeOption == VariantContextUtils.GenotypeMergeType.PRIORITIZE )
            priority = new ArrayList<String>(Arrays.asList(PRIORITY_STRING.split(",")));
        else
            priority = new ArrayList<String>(rodNames);

        if ( rodNames.size() != priority.size() )
            throw new UserException.BadArgumentValue("rod_priority_list", "The priority list must contain exactly one rod binding per ROD provided to the GATK: rodNames=" + rodNames + " priority=" + priority);

        if ( ! rodNames.containsAll(priority) )
            throw new UserException.BadArgumentValue("rod_priority_list", "Not all priority elements provided as input RODs: " + PRIORITY_STRING);
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null ) // RodWalkers can make funky map calls
            return 0;

        // get all of the vcf rods at this locus
        // Need to provide reference bases to simpleMerge starting at current locus
        Collection<VariantContext> vcs = tracker.getValues(VariantContext.class, context.getLocation());

        if ( sitesOnlyVCF ) {
            vcs = VariantContextUtils.sitesOnlyVariantContexts(vcs);
        }

        if ( ASSUME_IDENTICAL_SAMPLES ) {
            for ( final VariantContext vc : vcs ) {
                vcfWriter.add( vc, ref.getBase() );
            }
            
            return vcs.isEmpty() ? 0 : 1;
        }

        int numFilteredRecords = 0;
        for (VariantContext vc : vcs) {
            if (vc.filtersWereApplied() && vc.isFiltered())
                numFilteredRecords++;
        }

        if (minimumN > 1 && (vcs.size() - numFilteredRecords < minimumN))
            return 0;
        
        List<VariantContext> mergedVCs = new ArrayList<VariantContext>();
        if ( master ) {
            mergedVCs.add(VariantContextUtils.masterMerge(vcs, "master"));
        } else {
            Map<VariantContext.Type, List<VariantContext>> VCsByType = VariantContextUtils.separateVariantContextsByType(vcs);
            // iterate over the types so that it's deterministic
            for ( VariantContext.Type type : VariantContext.Type.values() ) {
                if ( VCsByType.containsKey(type) )
                    mergedVCs.add(VariantContextUtils.simpleMerge(getToolkit().getGenomeLocParser(), VCsByType.get(type),
                            priority, filteredRecordsMergeType, genotypeMergeOption, true, printComplexMerges,
                            ref.getBase(), SET_KEY, filteredAreUncalled, MERGE_INFO_WITH_MAX_AC));
            }
        }

        for ( VariantContext mergedVC : mergedVCs ) {
            // only operate at the start of events
            if ( mergedVC == null )
                continue;

            HashMap<String, Object> attributes = new HashMap<String, Object>(mergedVC.getAttributes());
            // re-compute chromosome counts
            VariantContextUtils.calculateChromosomeCounts(mergedVC, attributes, false);
            VariantContext annotatedMergedVC = VariantContext.modifyAttributes(mergedVC, attributes);
            if ( minimalVCF )
                annotatedMergedVC = VariantContextUtils.pruneVariantContext(annotatedMergedVC, Arrays.asList(SET_KEY));
            vcfWriter.add(annotatedMergedVC, ref.getBase());
        }

        return vcs.isEmpty() ? 0 : 1;
    }

    public Integer reduceInit() {
        return 0;
    }

    public Integer reduce(Integer counter, Integer sum) {
        return counter + sum;
    }

    public void onTraversalDone(Integer sum) {}
}

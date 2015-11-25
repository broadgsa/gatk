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

package org.broadinstitute.gatk.tools.walkers.filters;

import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.Reference;
import org.broadinstitute.gatk.engine.walkers.RodWalker;
import org.broadinstitute.gatk.engine.walkers.Window;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.engine.GATKVCFUtils;
import htsjdk.variant.vcf.*;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.util.*;


/**
 * Filter variant calls based on INFO and FORMAT annotations
 *
 * <p>
 * This tool is designed for hard-filtering variant calls based on certain criteria.
 * Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records
 * will be preserved in the output unless their removal is requested in the command line. </p>
 *
 * <h3>Input</h3>
 * <p>
 * A variant set to filter.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A filtered VCF.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T VariantFiltration \
 *   -R reference.fasta \
 *   -o output.vcf \
 *   --variant input.vcf \
 *   --filterExpression "AB < 0.2 || MQ0 > 50" \
 *   --filterName "Nov09filters" \
 *   --mask mask.vcf \
 *   --maskName InDel
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_VARMANIP, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-50,stop=50))
public class VariantFiltration extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * Any variant which overlaps entries from the provided mask rod will be filtered. If the user wants logic to be reversed,
     * i.e. filter variants that do not overlap with provided mask, then argument -filterNotInMask can be used.
     * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
     */
    @Input(fullName="mask", shortName="mask", doc="Input ROD mask", required=false)
    public RodBinding<Feature> mask;

    @Output(doc="File to which variants should be written")
    protected VariantContextWriter writer = null;

    /**
     * VariantFiltration accepts any number of JEXL expressions (so you can have two named filters by using
     * --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2").
     */
    @Argument(fullName="filterExpression", shortName="filter", doc="One or more expression used with INFO fields to filter", required=false)
    protected ArrayList<String> filterExpressions = new ArrayList<String>();

    /**
     * This name is put in the FILTER field for variants that get filtered.  Note that there must be a 1-to-1 mapping between filter expressions and filter names.
     */
    @Argument(fullName="filterName", shortName="filterName", doc="Names to use for the list of filters", required=false)
    protected ArrayList<String> filterNames = new ArrayList<String>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     * VariantFiltration will add the sample-level FT tag to the FORMAT field of filtered samples (this does not affect the record's FILTER tag).
     * One can filter normally based on most fields (e.g. "GQ < 5.0"), but the GT (genotype) field is an exception. We have put in convenience
     * methods so that one can now filter out hets ("isHet == 1"), refs ("isHomRef == 1"), or homs ("isHomVar == 1"). Also available are
     * expressions isCalled, isNoCall, isMixed, and isAvailable, in accordance with the methods of the Genotype object.
     */
    @Argument(fullName="genotypeFilterExpression", shortName="G_filter", doc="One or more expression used with FORMAT (sample/genotype-level) fields to filter (see documentation guide for more info)", required=false)
    protected ArrayList<String> genotypeFilterExpressions = new ArrayList<String>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     */
    @Argument(fullName="genotypeFilterName", shortName="G_filterName", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", required=false)
    protected ArrayList<String> genotypeFilterNames = new ArrayList<String>();

    /**
     * Works together with the --clusterWindowSize argument.
     */
    @Argument(fullName="clusterSize", shortName="cluster", doc="The number of SNPs which make up a cluster", required=false)
    protected Integer clusterSize = 3;

    /**
     * Works together with the --clusterSize argument.  To disable the clustered SNP filter, set this value to less than 1.
     */
    @Argument(fullName="clusterWindowSize", shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs", required=false)
    protected Integer clusterWindow = 0;

    @Argument(fullName="maskExtension", shortName="maskExtend", doc="How many bases beyond records from a provided 'mask' rod should variants be filtered", required=false)
    protected Integer maskExtension = 0;

    /**
     * When using the -mask argument, the maskName will be annotated in the variant record.
     * Note that when using the -filterNotInMask argument to reverse the masking logic,
     * it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
     */
    @Argument(fullName="maskName", shortName="maskName", doc="The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call", required=false)
    protected String maskName = "Mask";

    /**
     * By default, if the -mask argument is used, any variant falling in a mask will be filtered.
     * If this argument is used, logic is reversed, and variants falling outside a given mask will be filtered.
     * Use case is, for example, if we have an interval list or BED file with "good" sites.
     * Note that it is up to the user to adapt the name of the mask to make it clear that the reverse logic was used
     * (e.g. if masking against Hapmap, use -maskName=hapmap for the normal masking and -maskName=not_hapmap for the reverse masking).
     */
    @Argument(fullName="filterNotInMask", shortName="filterNotInMask", doc="Filter records NOT in given input mask.", required=false)
    protected boolean filterRecordsNotInMask = false;

    /**
     * By default, if JEXL cannot evaluate your expression for a particular record because one of the annotations is not present, the whole expression evaluates as PASSing.
     * Use this argument to have it evaluate as failing filters instead for these cases.
     */
    @Argument(fullName="missingValuesInExpressionsShouldEvaluateAsFailing", doc="When evaluating the JEXL expressions, missing values should be considered failing the expression", required=false)
    protected Boolean failMissingValues = false;

    /**
     * Invalidate previous filters applied to the VariantContext, applying only the filters here
     */
    @Argument(fullName="invalidatePreviousFilters",doc="Remove previous filters applied to the VCF",required=false)
    boolean invalidatePrevious = false;

    /**
     * Invert the selection criteria for --filterExpression
     */
    @Argument(fullName="invertFilterExpression", shortName="invfilter", doc="Invert the selection criteria for --filterExpression", required=false)
    protected boolean invertFilterExpression = false;

    /**
     * Invert the selection criteria for --genotypeFilterExpression
     */
    @Argument(fullName="invertGenotypeFilterExpression", shortName="invG_filter", doc="Invert the selection criteria for --genotypeFilterExpression", required=false)
    protected boolean invertGenotypeFilterExpression = false;

    /**
     * If this argument is provided, set filtered genotypes to no-call (./.).
     */
    @Argument(fullName="setFilteredGtToNocall", required=false, doc="Set filtered genotypes to no-call")
    private boolean setFilteredGenotypesToNocall = false;

    // JEXL expressions for the filters
    List<VariantContextUtils.JexlVCMatchExp> filterExps;
    List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps;

    public static final String CLUSTERED_SNP_FILTER_NAME = "SnpCluster";
    private ClusteredSnps clusteredSNPs = null;
    private GenomeLoc previousMaskPosition = null;

    // the structures necessary to initialize and maintain a windowed context
    private FiltrationContextWindow variantContextWindow;
    private static final int WINDOW_SIZE = 10;  // 10 variants on either end of the current one
    private ArrayList<FiltrationContext> windowInitializer = new ArrayList<FiltrationContext>();

    private final List<Allele> diploidNoCallAlleles = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    /**
     * Prepend inverse phrase to description if --invertFilterExpression
     *
     * @param description the description
     * @return the description with inverse prepended if --invert_filter_expression
     */
    private String possiblyInvertFilterExpression( String description ){
        if ( invertFilterExpression )
            description = "Inverse of: " + description;
        return description;
    }

    private void initializeVcfWriter() {

        final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));

        if ( clusterWindow > 0 )
            hInfo.add(new VCFFilterHeaderLine(CLUSTERED_SNP_FILTER_NAME, "SNPs found in clusters"));

        if ( !genotypeFilterExps.isEmpty() )
            hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));

        try {
            for ( VariantContextUtils.JexlVCMatchExp exp : filterExps )
                hInfo.add(new VCFFilterHeaderLine(exp.name, possiblyInvertFilterExpression(exp.exp.toString())));
            for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps )
                hInfo.add(new VCFFilterHeaderLine(exp.name, possiblyInvertFilterExpression(exp.exp.toString())));

            if ( mask.isBound() ) {
                if (filterRecordsNotInMask)
                    hInfo.add(new VCFFilterHeaderLine(maskName, "Doesn't overlap a user-input mask"));
                else hInfo.add(new VCFFilterHeaderLine(maskName, "Overlaps a user-input mask"));
            }
        } catch (IllegalArgumentException e) {
            throw new UserException.BadInput(e.getMessage());
        }

        writer.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));
    }

    public void initialize() {
        if ( clusterWindow > 0 )
            clusteredSNPs = new ClusteredSnps(getToolkit().getGenomeLocParser(),clusterSize, clusterWindow);

        if ( maskExtension < 0 )
             throw new UserException.BadArgumentValue("maskExtension", "negative values are not allowed");

        if (filterRecordsNotInMask && !mask.isBound())
            throw new UserException.BadArgumentValue("filterNotInMask","argument not allowed if mask argument is not provided");
        filterExps = VariantContextUtils.initializeMatchExps(filterNames, filterExpressions);
        genotypeFilterExps = VariantContextUtils.initializeMatchExps(genotypeFilterNames, genotypeFilterExpressions);

        VariantContextUtils.engine.get().setSilent(true);

        initializeVcfWriter();
    }

    public Integer reduceInit() { return 0; }

    /**
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        Collection<VariantContext> VCs = tracker.getValues(variantCollection.variants, context.getLocation());

        // is there a SNP mask present?
        boolean hasMask = (tracker.hasValues(mask) && !filterRecordsNotInMask) || (filterRecordsNotInMask && !tracker.hasValues(mask));
        if ( hasMask )
            previousMaskPosition = ref.getLocus();  // multi-base masks will get triggered over all bases of the mask

        for ( VariantContext vc : VCs ) {

            if ( invalidatePrevious ) {
                vc = (new VariantContextBuilder(vc)).filters(new HashSet<String>()).make();
            }

            // filter based on previous mask position
            vc = addMaskIfCoversVariant(vc, previousMaskPosition, maskName, maskExtension, false);

            FiltrationContext varContext = new FiltrationContext(ref, vc);

            // if we're still initializing the context, do so
            if ( windowInitializer != null ) {

                // if this is a mask position, filter previous records
                if ( hasMask ) {
                    for ( FiltrationContext prevVC : windowInitializer )
                        prevVC.setVariantContext(addMaskIfCoversVariant(prevVC.getVariantContext(), ref.getLocus(), maskName, maskExtension, true));
                }

                windowInitializer.add(varContext);
                if ( windowInitializer.size() == WINDOW_SIZE ) {
                    variantContextWindow = new FiltrationContextWindow(windowInitializer);
                    windowInitializer = null;
                }
            } else {

                // if this is a mask position, filter previous records
                if ( hasMask ) {
                    for ( FiltrationContext prevVC : variantContextWindow.getWindow(10, 10) ) {
                        if ( prevVC != null )
                            prevVC.setVariantContext(addMaskIfCoversVariant(prevVC.getVariantContext(), ref.getLocus(), maskName, maskExtension, true));
                    }
                }

                variantContextWindow.moveWindow(varContext);
                filter();
            }
        }

        return 1;
    }

    /**
     * Helper function to check if a mask covers the variant location.
     *
     * @param vc variant context
     * @param genomeLoc genome location
     * @param maskName name of the mask
     * @param maskExtension bases beyond the mask
     * @param vcBeforeLoc if true, variant context is before the genome location; if false, the converse is true.
     * @return true if the genome location is within the extended mask area, false otherwise
     */
    protected static boolean doesMaskCoverVariant(VariantContext vc, GenomeLoc genomeLoc, String maskName, int maskExtension, boolean vcBeforeLoc) {
        boolean logic = genomeLoc != null &&                                        // have a location
                genomeLoc.getContig().equals(vc.getChr()) &&                        // it's on the same contig
                (vc.getFilters() == null || !vc.getFilters().contains(maskName));   // the filter hasn't already been applied
        if ( logic ) {
            if (vcBeforeLoc)
                return genomeLoc.getStart() - vc.getEnd() <= maskExtension;  // it's within the mask area (multi-base VCs that overlap this site will always give a negative distance)
            else
                return vc.getStart() - genomeLoc.getStop() <= maskExtension;
        } else {
            return false;
        }
    }

    /**
     * Add mask to variant context filters if it covers the it's location
     *
     * @param vc VariantContext
     * @param genomeLoc genome location
     * @param maskName name of the mask
     * @param maskExtension bases beyond the mask
     * @param locStart if true, start at genome location and end at VariantContext. If false, do the opposite.
     * @return VariantContext with the mask added if the VariantContext is within the extended mask area
     */
    private VariantContext addMaskIfCoversVariant(VariantContext vc, GenomeLoc genomeLoc, String maskName, int maskExtension, boolean locStart) {
        if (doesMaskCoverVariant(vc, genomeLoc, maskName, maskExtension, locStart) ) {
            Set<String> filters = new LinkedHashSet<String>(vc.getFilters());
            filters.add(maskName);
            vc = new VariantContextBuilder(vc).filters(filters).make();
        }

        return vc;
    }

    private void filter() {
        // get the current context
        FiltrationContext context = variantContextWindow.getContext();
        if ( context == null )
            return;

        final VariantContext vc = context.getVariantContext();
        final VariantContextBuilder builder = new VariantContextBuilder(vc);

        // make new Genotypes based on filters
        if ( !genotypeFilterExps.isEmpty() || setFilteredGenotypesToNocall ) {
            GenotypesContext genotypes = GenotypesContext.create(vc.getGenotypes().size());

            // for each genotype, check filters then create a new object
            for ( final Genotype g : vc.getGenotypes() ) {
                if ( g.isCalled() ) {
                    final List<String> filters = new ArrayList<String>();
                    if ( g.isFiltered() ) filters.add(g.getFilters());

                    // Add if expression filters the variant context
                    for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps ) {
                        if ( Utils.invertLogic(VariantContextUtils.match(vc, g, exp), invertGenotypeFilterExpression) )
                            filters.add(exp.name);
                    }

                    // if sample is filtered and --setFilteredGtToNocall, set genotype to non-call
                    if ( !filters.isEmpty() && setFilteredGenotypesToNocall )
                        genotypes.add(new GenotypeBuilder(g).filters(filters).alleles(diploidNoCallAlleles).make());
                    else
                        genotypes.add(new GenotypeBuilder(g).filters(filters).make());
                } else {
                    genotypes.add(g);
                }
            }

            builder.genotypes(genotypes);
        }

        // make a new variant context based on filters
        Set<String> filters = new LinkedHashSet<String>(vc.getFilters());

        // test for clustered SNPs if requested
        if ( clusteredSNPs != null && clusteredSNPs.filter(variantContextWindow) )
            filters.add(CLUSTERED_SNP_FILTER_NAME);

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps ) {
            try {
                if ( Utils.invertLogic(VariantContextUtils.match(vc, exp), invertFilterExpression) )
                    filters.add(exp.name);
            } catch (Exception e) {
                // do nothing unless specifically asked to; it just means that the expression isn't defined for this context
                if ( failMissingValues  )
                    filters.add(exp.name);                         
            }
        }

        if ( filters.isEmpty() )
            builder.passFilters();
        else
            builder.filters(filters);

        writer.add(builder.make());
    }

    public Integer reduce(Integer value, Integer sum) {
        return sum + value;
    }

    /**
     * Tell the user the number of loci processed and close out the new variants file.
     *
     * @param result  the number of loci seen.
     */
    public void onTraversalDone(Integer result) {
        // move the window over so that we can filter the last few variants
        if ( windowInitializer != null ) {
            while ( windowInitializer.size() < WINDOW_SIZE )
                windowInitializer.add(null);
            variantContextWindow = new FiltrationContextWindow(windowInitializer);
        }
        for (int i=0; i < WINDOW_SIZE; i++) {
            variantContextWindow.moveWindow(null);
            filter();
        }
    }
}

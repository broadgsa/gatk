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

package org.broadinstitute.sting.gatk.walkers.filters;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.util.*;


/**
 * Filters variant calls using a number of user-selectable, parameterizable criteria.
 *
 * <p>
 * VariantFiltration is a GATK tool for hard-filtering variant calls based on certain criteria.
 * Records are hard-filtered by changing the value in the FILTER field to something other than PASS.
 *
 * <h2>Input</h2>
 * <p>
 * A variant set to filter.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * A filtered VCF.
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantFiltration \
 *   -o output.vcf \
 *   --variant input.vcf \
 *   --filterExpression "AB < 0.2 || MQ0 > 50" \
 *   --filterName "Nov09filters" \
 *   --mask mask.vcf \
 *   --maskName InDel
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = "Variant Evaluation and Manipulation Tools", extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-50,stop=50))
public class VariantFiltration extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * Any variant which overlaps entries from the provided mask rod will be filtered.
     */
    @Input(fullName="mask", doc="Input ROD mask", required=false)
    public RodBinding<Feature> mask;

    @Output(doc="File to which variants should be written", required=true)
    protected VariantContextWriter writer = null;

    /**
     * VariantFiltration accepts any number of JEXL expressions (so you can have two named filters by using
     * --filterName One --filterExpression "X < 1" --filterName Two --filterExpression "X > 2").
     */
    @Argument(fullName="filterExpression", shortName="filter", doc="One or more expression used with INFO fields to filter", required=false)
    protected ArrayList<String> FILTER_EXPS = new ArrayList<String>();

    /**
     * This name is put in the FILTER field for variants that get filtered.  Note that there must be a 1-to-1 mapping between filter expressions and filter names.
     */
    @Argument(fullName="filterName", shortName="filterName", doc="Names to use for the list of filters", required=false)
    protected ArrayList<String> FILTER_NAMES = new ArrayList<String>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     * VariantFiltration will add the sample-level FT tag to the FORMAT field of filtered samples (this does not affect the record's FILTER tag).
     * One can filter normally based on most fields (e.g. "GQ < 5.0"), but the GT (genotype) field is an exception. We have put in convenience
     * methods so that one can now filter out hets ("isHet == 1"), refs ("isHomRef == 1"), or homs ("isHomVar == 1").
     */
    @Argument(fullName="genotypeFilterExpression", shortName="G_filter", doc="One or more expression used with FORMAT (sample/genotype-level) fields to filter (see wiki docs for more info)", required=false)
    protected ArrayList<String> GENOTYPE_FILTER_EXPS = new ArrayList<String>();

    /**
     * Similar to the INFO field based expressions, but used on the FORMAT (genotype) fields instead.
     */
    @Argument(fullName="genotypeFilterName", shortName="G_filterName", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", required=false)
    protected ArrayList<String> GENOTYPE_FILTER_NAMES = new ArrayList<String>();

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
    protected Integer MASK_EXTEND = 0;
    @Argument(fullName="maskName", shortName="maskName", doc="The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call", required=false)
    protected String MASK_NAME = "Mask";

    /**
     * By default, if JEXL cannot evaluate your expression for a particular record because one of the annotations is not present, the whole expression evaluates as PASSing.
     * Use this argument to have it evaluate as failing filters instead for these cases.
     */
    @Argument(fullName="missingValuesInExpressionsShouldEvaluateAsFailing", doc="When evaluating the JEXL expressions, missing values should be considered failing the expression", required=false)
    protected Boolean FAIL_MISSING_VALUES = false;

    /**
     * Invalidate previous filters applied to the VariantContext, applying only the filters here
     */
    @Argument(fullName="invalidatePreviousFilters",doc="Remove previous filters applied to the VCF",required=false)
    boolean invalidatePrevious = false;

    // JEXL expressions for the filters
    List<VariantContextUtils.JexlVCMatchExp> filterExps;
    List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps;

    public static final String CLUSTERED_SNP_FILTER_NAME = "SnpCluster";
    private ClusteredSnps clusteredSNPs = null;
    private GenomeLoc previousMaskPosition = null;

    // the structures necessary to initialize and maintain a windowed context
    private FiltrationContextWindow variantContextWindow;
    private static final int windowSize = 10;  // 10 variants on either end of the current one
    private ArrayList<FiltrationContext> windowInitializer = new ArrayList<FiltrationContext>();

    private void initializeVcfWriter() {

        final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit(), inputNames));

        if ( clusterWindow > 0 )
            hInfo.add(new VCFFilterHeaderLine(CLUSTERED_SNP_FILTER_NAME, "SNPs found in clusters"));

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
        for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));

        if ( genotypeFilterExps.size() > 0 )
            hInfo.add(VCFStandardHeaderLines.getFormatLine(VCFConstants.GENOTYPE_FILTER_KEY));

        if ( mask.isBound() ) {
            hInfo.add(new VCFFilterHeaderLine(MASK_NAME, "Overlaps a user-input mask"));
        }

        writer.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));
    }

    public void initialize() {
        if ( clusterWindow > 0 )
            clusteredSNPs = new ClusteredSnps(getToolkit().getGenomeLocParser(),clusterSize, clusterWindow);

        if ( MASK_EXTEND < 0 )
             throw new UserException.BadArgumentValue("maskExtension", "negative values are not allowed");

        filterExps = VariantContextUtils.initializeMatchExps(FILTER_NAMES, FILTER_EXPS);
        genotypeFilterExps = VariantContextUtils.initializeMatchExps(GENOTYPE_FILTER_NAMES, GENOTYPE_FILTER_EXPS);

        VariantContextUtils.engine.setSilent(true);

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
        boolean hasMask = tracker.hasValues(mask);
        if ( hasMask )
            previousMaskPosition = ref.getLocus();  // multi-base masks will get triggered over all bases of the mask

        for ( VariantContext vc : VCs ) {

            if ( invalidatePrevious ) {
                vc = (new VariantContextBuilder(vc)).filters(new HashSet<String>()).make();
            }
            // filter based on previous mask position
            if ( previousMaskPosition != null &&                                       // we saw a previous mask site
                 previousMaskPosition.getContig().equals(vc.getChr()) &&               // it's on the same contig
                 vc.getStart() - previousMaskPosition.getStop() <= MASK_EXTEND &&      // it's within the mask area (multi-base masks that overlap this site will always give a negative distance)
                 (vc.getFilters() == null || !vc.getFilters().contains(MASK_NAME)) ) { // the filter hasn't already been applied
                Set<String> filters = new LinkedHashSet<String>(vc.getFilters());
                filters.add(MASK_NAME);
                vc = new VariantContextBuilder(vc).filters(filters).make();
            }

            FiltrationContext varContext = new FiltrationContext(ref, vc);

            // if we're still initializing the context, do so
            if ( windowInitializer != null ) {

                // if this is a mask position, filter previous records
                if ( hasMask ) {
                    for ( FiltrationContext prevVC : windowInitializer )
                        prevVC.setVariantContext(checkMaskForPreviousLocation(prevVC.getVariantContext(), ref.getLocus()));
                }

                windowInitializer.add(varContext);
                if ( windowInitializer.size() == windowSize ) {
                    variantContextWindow = new FiltrationContextWindow(windowInitializer);
                    windowInitializer = null;
                }
            } else {

                // if this is a mask position, filter previous records
                if ( hasMask ) {
                    for ( FiltrationContext prevVC : variantContextWindow.getWindow(10, 10) ) {
                        if ( prevVC != null )
                            prevVC.setVariantContext(checkMaskForPreviousLocation(prevVC.getVariantContext(), ref.getLocus()));
                    }
                }

                variantContextWindow.moveWindow(varContext);
                filter();
            }
        }

        return 1;
    }

    private VariantContext checkMaskForPreviousLocation(VariantContext vc, GenomeLoc maskLoc) {
        if ( maskLoc.getContig().equals(vc.getChr()) &&               // it's on the same contig
             maskLoc.getStart() - vc.getEnd() <= MASK_EXTEND &&       // it's within the mask area (multi-base VCs that overlap this site will always give a negative distance)
             (vc.getFilters() == null || !vc.getFilters().contains(MASK_NAME)) ) { // the filter hasn't already been applied
            Set<String> filters = new LinkedHashSet<String>(vc.getFilters());
            filters.add(MASK_NAME);
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
        if ( genotypeFilterExps.size() > 0 ) {
            GenotypesContext genotypes = GenotypesContext.create(vc.getGenotypes().size());

            // for each genotype, check filters then create a new object
            for ( final Genotype g : vc.getGenotypes() ) {
                if ( g.isCalled() ) {
                    final List<String> filters = new ArrayList<String>();
                    if ( g.isFiltered() ) filters.add(g.getFilters());

                    for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps ) {
                        if ( VariantContextUtils.match(vc, g, exp) )
                            filters.add(exp.name);
                    }

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
                if ( VariantContextUtils.match(vc, exp) )
                    filters.add(exp.name);
            } catch (Exception e) {
                // do nothing unless specifically asked to; it just means that the expression isn't defined for this context
                if ( FAIL_MISSING_VALUES )
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
            while ( windowInitializer.size() < windowSize )
                windowInitializer.add(null);
            variantContextWindow = new FiltrationContextWindow(windowInitializer);
        }
        for (int i=0; i < windowSize; i++) {
            variantContextWindow.moveWindow(null);
            filter();
        }
    }
}

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

package org.broadinstitute.sting.gatk.walkers.filters;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.CommandLineUtils;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.util.*;


/**
 * Filters variant calls using a number of user-selectable, parameterizable criteria.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type=VariantContext.class))
@Reference(window=@Window(start=-50,stop=50))
public class VariantFiltrationWalker extends RodWalker<Integer, Integer> {

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter writer = null;

    @Argument(fullName="filterExpression", shortName="filter", doc="One or more expression used with INFO fields to filter (see wiki docs for more info)", required=false)
    protected ArrayList<String> FILTER_EXPS = new ArrayList<String>();
    @Argument(fullName="filterName", shortName="filterName", doc="Names to use for the list of filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", required=false)
    protected ArrayList<String> FILTER_NAMES = new ArrayList<String>();

    @Argument(fullName="genotypeFilterExpression", shortName="G_filter", doc="One or more expression used with FORMAT (sample/genotype-level) fields to filter (see wiki docs for more info)", required=false)
    protected ArrayList<String> GENOTYPE_FILTER_EXPS = new ArrayList<String>();
    @Argument(fullName="genotypeFilterName", shortName="G_filterName", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", required=false)
    protected ArrayList<String> GENOTYPE_FILTER_NAMES = new ArrayList<String>();

    @Argument(fullName="clusterSize", shortName="cluster", doc="The number of SNPs which make up a cluster (see also --clusterWindowSize); [default:3]", required=false)
    protected Integer clusterSize = 3;
    @Argument(fullName="clusterWindowSize", shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs (to disable the clustered SNP filter, set this value to less than 1); [default:0]", required=false)
    protected Integer clusterWindow = 0;

    @Argument(fullName="maskName", shortName="mask", doc="The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call; [default:'Mask']", required=false)
    protected String MASK_NAME = "Mask";

    @Argument(fullName = "NO_HEADER", shortName = "NO_HEADER", doc = "Don't output the usual VCF header tag with the command line. FOR DEBUGGING PURPOSES ONLY. This option is required in order to pass integration tests.", required = false)
    protected Boolean NO_VCF_HEADER_LINE = false;

    // JEXL expressions for the filters
    List<VariantContextUtils.JexlVCMatchExp> filterExps;
    List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps;

    public static final String CLUSTERED_SNP_FILTER_NAME = "SnpCluster";
    private ClusteredSnps clusteredSNPs = null;

    // the structures necessary to initialize and maintain a windowed context
    private FiltrationContextWindow variantContextWindow;
    private static final int windowSize = 10;  // 10 variants on either end of the current one
    private ArrayList<FiltrationContext> windowInitializer = new ArrayList<FiltrationContext>();


    private void initializeVcfWriter(VariantContext vc) {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));

        if ( clusterWindow > 0 )
            hInfo.add(new VCFFilterHeaderLine(CLUSTERED_SNP_FILTER_NAME, "SNPs found in clusters"));

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
        for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));

        if ( genotypeFilterExps.size() > 0 )
            hInfo.add(new VCFFormatHeaderLine(VCFConstants.GENOTYPE_FILTER_KEY, 1, VCFHeaderLineType.String, "Genotype-level filter"));

        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            if ( source.getReferenceOrderedData().getName().equals("mask") ) {
                hInfo.add(new VCFFilterHeaderLine(MASK_NAME, "Overlaps a user-input mask"));
                break;
            }
        }

        if ( !NO_VCF_HEADER_LINE ) {
            Set<Object> args = new HashSet<Object>();
            args.add(this);
            hInfo.add(new VCFHeaderLine("VariantFiltration", "\"" + CommandLineUtils.createApproximateCommandLineArgumentString(getToolkit(), args, getClass()) + "\""));
        }

        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(vc.getSampleNames())));
    }

    public void initialize() {
        if ( clusterWindow > 0 )
            clusteredSNPs = new ClusteredSnps(clusterSize, clusterWindow);

        filterExps = VariantContextUtils.initializeMatchExps(FILTER_NAMES, FILTER_EXPS);
        genotypeFilterExps = VariantContextUtils.initializeMatchExps(GENOTYPE_FILTER_NAMES, GENOTYPE_FILTER_EXPS);

        VariantContextUtils.engine.setSilent(true);
    }

    public Integer reduceInit() { return 0; }

    /**
     * For each site of interest, rescore the genotype likelihoods by applying the specified feature set.
     *
     * @param tracker  the meta-data tracker
     * @param ref      the reference base
     * @param context  the context for the given locus
     * @return 1 if the locus was successfully processed, 0 if otherwise
     */
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        List<Object> rods = tracker.getReferenceMetaData("variant");
        // ignore places where we don't have a variant
        if ( rods.size() == 0 )
            return 0;

        VariantContext vc = VariantContextAdaptors.toVariantContext("variant", rods.get(0), ref);
        FiltrationContext varContext = new FiltrationContext(tracker, ref, vc);

        // if we're still initializing the context, do so
        if ( windowInitializer != null ) {
            windowInitializer.add(varContext);
            if ( windowInitializer.size() == windowSize ) {
                variantContextWindow = new FiltrationContextWindow(windowInitializer);
                windowInitializer = null;
            }
        } else {
            variantContextWindow.moveWindow(varContext);
            filter();
        }

        return 1;
    }

    private void filter() {
        // get the current context
        FiltrationContext context = variantContextWindow.getContext();
        if ( context == null )
            return;

        VariantContext vc = context.getVariantContext();

        // make new Genotypes based on filters
        Map<String, Genotype> genotypes;
        if ( genotypeFilterExps.size() == 0 ) {
            genotypes = vc.getGenotypes();
        } else {
            genotypes = new HashMap<String, Genotype>(vc.getGenotypes().size());

            // for each genotype, check filters then create a new object
            for ( Map.Entry<String, Genotype> genotype : vc.getGenotypes().entrySet() ) {

                Genotype g = genotype.getValue();
                Set<String> filters = new LinkedHashSet<String>(g.getFilters());

                for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps ) {
                    if ( VariantContextUtils.match(vc, g, exp) )
                        filters.add(exp.name);
                }

                genotypes.put(genotype.getKey(), new Genotype(genotype.getKey(), g.getAlleles(), g.getNegLog10PError(), filters, g.getAttributes(), g.genotypesArePhased()));
            }
        }

        // make a new variant context based on filters
        Set<String> filters = new LinkedHashSet<String>(vc.getFilters());

        // test for SNP mask, if present
        List<Object> mask = context.getTracker().getReferenceMetaData("mask");
        if ( mask.size() > 0 )
            filters.add(MASK_NAME);

        // test for clustered SNPs if requested
        if ( clusteredSNPs != null && clusteredSNPs.filter(variantContextWindow) )
            filters.add(CLUSTERED_SNP_FILTER_NAME);

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps ) {
            if ( VariantContextUtils.match(vc, exp) )
                filters.add(exp.name);
        }

        VariantContext filteredVC = new VariantContext(vc.getName(), vc.getChr(), vc.getStart(), vc.getEnd(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), filters, vc.getAttributes());

        writeVCF(filteredVC, context.getReferenceContext().getBase());
    }

    private void writeVCF(VariantContext vc, byte ref) {
        if ( writer == null )
            initializeVcfWriter(vc);
        writer.add(vc, ref);
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

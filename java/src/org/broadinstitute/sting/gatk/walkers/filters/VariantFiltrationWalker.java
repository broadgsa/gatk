package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import org.apache.commons.jexl.*;


/**
 * Filters variant calls using a number of user-selectable, parameterizable criteria.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant",type= ReferenceOrderedDatum.class))
public class VariantFiltrationWalker extends RodWalker<Integer, Integer> {

    @Argument(fullName="filterExpression", shortName="filter", doc="One or more expression used with INFO fields to filter (see wiki docs for more info)", required=false)
    protected String[] FILTER_EXPS = new String[]{};
    @Argument(fullName="filterName", shortName="filterName", doc="Names to use for the list of filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", required=false)
    protected String[] FILTER_NAMES = new String[]{};

    @Argument(fullName="clusterSize", shortName="cluster", doc="The number of SNPs which make up a cluster (see also --clusterWindowSize)", required=false)
    protected Integer clusterSize = 3;
    @Argument(fullName="clusterWindowSize", shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs (to disable the clustered SNP filter, set this value to less than 1)", required=false)
    protected Integer clusterWindow = 0;

    @Argument(fullName="maskName", shortName="mask", doc="The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call", required=false)
    protected String MASK_NAME = "Mask";

    // JEXL expressions for the filters
    List<VariantContextUtils.JexlVCMatchExp> filterExps;

    public static final String CLUSTERED_SNP_FILTER_NAME = "SnpCluster";

    private VCFWriter writer = null;

    private ClusteredSnps clusteredSNPs = null;

    class FilterExp {
        String name;
        String expStr;
        Expression exp;

        public FilterExp(String name, String str, Expression exp) {
            this.name = name;
            this.expStr = str;
            this.exp = exp;
        }
    }

    // the structures necessary to initialize and maintain a windowed context
    private VariantContextWindow variantContextWindow;
    private static final int windowSize = 10;  // 10 variants on either end of the current one
    private ArrayList<FiltrationContext> windowInitializer = new ArrayList<FiltrationContext>();


    private void initializeVcfWriter(VCFRecord rec) {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantFiltration"));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        if ( clusterWindow > 0 )
            hInfo.add(new VCFFilterHeaderLine(CLUSTERED_SNP_FILTER_NAME, "SNPs found in clusters"));

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps ) {
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
        }

        List<ReferenceOrderedDataSource> dataSources = getToolkit().getRodDataSources();
        for ( ReferenceOrderedDataSource source : dataSources ) {
            if ( source.getReferenceOrderedData().getName().equals("mask") ) {
                hInfo.add(new VCFFilterHeaderLine(MASK_NAME, "Overlaps a user-input mask"));
                break;
            }
        }

        writer = new VCFWriter(out);
        writer.writeHeader(new VCFHeader(hInfo, new TreeSet<String>(Arrays.asList(rec.getSampleNames()))));
    }

    public void initialize() {
        if ( clusterWindow > 0 )
            clusteredSNPs = new ClusteredSnps(clusterSize, clusterWindow);

        filterExps = VariantContextUtils.initializeMatchExps(FILTER_NAMES, FILTER_EXPS);
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

        RODRecordList rods = tracker.getTrackData("variant", null);
        // ignore places where we don't have a variant
        if ( rods == null || rods.size() == 0 )
            return 0;

        VariantContext vc = VariantContextAdaptors.toVariantContext("variant", rods.get(0));
        FiltrationContext varContext = new FiltrationContext(tracker, ref, vc);

        // if we're still initializing the context, do so
        if ( windowInitializer != null ) {
            windowInitializer.add(varContext);
            if ( windowInitializer.size() == windowSize ) {
                variantContextWindow = new VariantContextWindow(windowInitializer);
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

        StringBuilder filterString = new StringBuilder();

        // test for SNP mask, if present
        RODRecordList mask = context.getTracker().getTrackData("mask", null);
        if ( mask != null && mask.size() > 0 )
            addFilter(filterString, MASK_NAME);

        // test for clustered SNPs if requested
        if ( clusteredSNPs != null && clusteredSNPs.filter(variantContextWindow) )
            addFilter(filterString, CLUSTERED_SNP_FILTER_NAME);

        VariantContext vc = context.getVariantContext();
        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps ) {
            if ( VariantContextUtils.match(vc, exp) )
                addFilter(filterString, exp.name);
        }

        writeVCF(VariantContextAdaptors.toVCF(vc, context.getReferenceContext().getBase(), null, true), filterString);
    }

    private static void addFilter(StringBuilder sb, String filter) {
        if ( sb.length() > 0 )
            sb.append(";");
        sb.append(filter);
    }

    private void writeVCF(VCFRecord rec, StringBuilder filterString) {

        if ( filterString.length() != 0 ) {
            // if the record is already filtered, don't destroy those filters
            if ( rec.isFiltered() )
                filterString.append(";" + rec.getFilterString());
            rec.setFilterString(filterString.toString());
        }
        // otherwise, if it's not already filtered, set it to "passing filters"
        else if ( !rec.isFiltered() ) {
            rec.setFilterString(VCFRecord.PASSES_FILTERS);    
        }

        if ( writer == null )
            initializeVcfWriter(rec);
        writer.addRecord(rec);
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
            variantContextWindow = new VariantContextWindow(windowInitializer);
        }
        for (int i=0; i < windowSize; i++) {
            variantContextWindow.moveWindow(null);
            filter();
        }

        if ( writer != null )
            writer.close();
    }
}

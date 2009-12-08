package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import org.apache.commons.jexl.*;


/**
 * VariantFiltrationWalker filters variant calls in VCF format.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant",type= RodVCF.class))
public class VariantFiltrationWalker extends RodWalker<Integer, Integer> {
    @Argument(fullName="filterExpression", shortName="filter", doc="Expression used with INFO fields to filter (see wiki docs for more info)", required=false)
    protected String FILTER_STRING = null;
    @Argument(fullName="filterName", shortName="filterName", doc="The text to put in the FILTER field if a filter expression is provided and a variant call matches", required=false)
    protected String FILTER_NAME = "GATK_filter";

    @Argument(fullName="clusterSize", shortName="cluster", doc="The number of SNPs which make up a cluster (see also --clusterWindowSize)", required=false)
    protected Integer clusterSize = 3;
    @Argument(fullName="clusterWindowSize", shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs (to disable the clustered SNP filter, set this value to less than 1)", required=false)
    protected Integer clusterWindow = 0;

    @Argument(fullName="maskName", shortName="mask", doc="The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call", required=false)
    protected String MASK_NAME = "Mask";

    private VCFWriter writer = null;

    private ClusteredSnps clusteredSNPs = null;

    private Expression filterExpression = null;

    // the structures necessary to initialize and maintain a windowed context
    private VariantContextWindow variantContextWindow;
    private static final int windowSize = 10;  // 10 variants on either end of the current one
    private ArrayList<Pair<RefMetaDataTracker, RodVCF>> windowInitializer = new ArrayList<Pair<RefMetaDataTracker, RodVCF>>();


    private void initializeVcfWriter(RodVCF rod) {
        // setup the header fields
        Set<String> hInfo = new HashSet<String>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add("source=" + "VariantFiltration");
        hInfo.add("reference=" + getToolkit().getArguments().referenceFile.getName());

        VCFHeader header = new VCFHeader(hInfo, rod.getHeader().getGenotypeSamples());
        writer = new VCFWriter(header, out);
    }

    public void initialize() {
        if ( clusterWindow > 0 )
            clusteredSNPs = new ClusteredSnps(clusterSize, clusterWindow);

        try {
            if ( FILTER_STRING != null )
                filterExpression = ExpressionFactory.createExpression(FILTER_STRING);
        } catch (Exception e) {
            throw new StingException("Invalid expression used (" + FILTER_STRING + "). Please see the JEXL docs for correct syntax.");
        }
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

        RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData("variant", null);
        // ignore places where we don't have a variant
        if ( rods == null || rods.getRecords().size() == 0 )
            return 0;

        RodVCF variant = (RodVCF)rods.getRecords().get(0);
        Pair<RefMetaDataTracker, RodVCF> varContext = new Pair<RefMetaDataTracker, RodVCF>(tracker, variant);

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
        Pair<RefMetaDataTracker, RodVCF> context = variantContextWindow.getContext();
        if ( context == null )
            return;

        StringBuilder filterString = new StringBuilder();

        // test for SNP mask, if present
        RODRecordList<ReferenceOrderedDatum> mask = context.first.getTrackData("mask", null);
        if ( mask != null && mask.getRecords().size() > 0 )
            addFilter(filterString, MASK_NAME);

        // test for clustered SNPs if requested
        if ( clusteredSNPs != null && clusteredSNPs.filter(variantContextWindow) )
            addFilter(filterString, "SnpCluster");

        if ( filterExpression != null ) {
            Map<String, String> infoMap = new HashMap<String, String>(context.second.mCurrentRecord.getInfoValues());
            infoMap.put("QUAL", String.valueOf(context.second.mCurrentRecord.getQual()));

            JexlContext jContext = JexlHelper.createContext();
            jContext.setVars(infoMap);

            try {
                if ( (Boolean)filterExpression.evaluate(jContext) )
                    addFilter(filterString, FILTER_NAME);
            } catch (Exception e) {
                throw new StingException(e.getMessage());
            }
        }

        writeVCF(context.second, filterString);
    }

    private static void addFilter(StringBuilder sb, String filter) {
        if ( sb.length() > 0 )
            sb.append(";");
        sb.append(filter);
    }

    private void writeVCF(RodVCF variant, StringBuilder filterString) {
        VCFRecord rec = variant.mCurrentRecord;

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
            initializeVcfWriter(variant);
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

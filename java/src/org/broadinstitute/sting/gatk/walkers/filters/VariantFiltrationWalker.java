package org.broadinstitute.sting.gatk.walkers.filters;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
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

    @Argument(fullName="genotypeFilterExpression", shortName="G_filter", doc="One or more expression used with FORMAT (sample/genotype-level) fields to filter (see wiki docs for more info)", required=false)
    protected String[] GENOTYPE_FILTER_EXPS = new String[]{};
    @Argument(fullName="genotypeFilterName", shortName="G_filterName", doc="Names to use for the list of sample/genotype filters (must be a 1-to-1 mapping); this name is put in the FILTER field for variants that get filtered", required=false)
    protected String[] GENOTYPE_FILTER_NAMES = new String[]{};

    @Argument(fullName="clusterSize", shortName="cluster", doc="The number of SNPs which make up a cluster (see also --clusterWindowSize); [default:3]", required=false)
    protected Integer clusterSize = 3;
    @Argument(fullName="clusterWindowSize", shortName="window", doc="The window size (in bases) in which to evaluate clustered SNPs (to disable the clustered SNP filter, set this value to less than 1); [default:0]", required=false)
    protected Integer clusterWindow = 0;

    @Argument(fullName="maskName", shortName="mask", doc="The text to put in the FILTER field if a 'mask' rod is provided and overlaps with a variant call; [default:'Mask']", required=false)
    protected String MASK_NAME = "Mask";

    // JEXL expressions for the filters
    List<VariantContextUtils.JexlVCMatchExp> filterExps;
    List<VariantContextUtils.JexlVCMatchExp> genotypeFilterExps;

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
    private FiltrationContextWindow variantContextWindow;
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

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));
        for ( VariantContextUtils.JexlVCMatchExp exp : genotypeFilterExps )
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.exp.toString()));

        if ( genotypeFilterExps.size() > 0 )
            hInfo.add(new VCFFormatHeaderLine(VCFGenotypeRecord.GENOTYPE_FILTER_KEY, 1, VCFFormatHeaderLine.FORMAT_TYPE.String, "Genotype-level filter"));

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
        genotypeFilterExps = VariantContextUtils.initializeMatchExps(GENOTYPE_FILTER_NAMES, GENOTYPE_FILTER_EXPS);
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
                    if ( VariantContextUtils.match(g, exp) )
                        filters.add(exp.name);
                }

                genotypes.put(genotype.getKey(), new Genotype(genotype.getKey(), g.getAlleles(), g.getNegLog10PError(), filters, g.getAttributes(), g.genotypesArePhased()));
            }
        }

        // make a new variant context based on filters
        Set<String> filters = new LinkedHashSet<String>(vc.getFilters());

        // test for SNP mask, if present
        RODRecordList mask = context.getTracker().getTrackData("mask", null);
        if ( mask != null && mask.size() > 0 )
            filters.add(MASK_NAME);

        // test for clustered SNPs if requested
        if ( clusteredSNPs != null && clusteredSNPs.filter(variantContextWindow) )
            filters.add(CLUSTERED_SNP_FILTER_NAME);

        for ( VariantContextUtils.JexlVCMatchExp exp : filterExps ) {
            if ( VariantContextUtils.match(vc, exp) )
                filters.add(exp.name);
        }

        VariantContext filteredVC = new VariantContext(vc.getName(), vc.getLocation(), vc.getAlleles(), genotypes, vc.getNegLog10PError(), filters, vc.getAttributes());

        writeVCF(VariantContextAdaptors.toVCF(filteredVC, context.getReferenceContext().getBase(), null, true, genotypeFilterExps.size() > 0));
    }

    private void writeVCF(VCFRecord rec) {
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
            variantContextWindow = new FiltrationContextWindow(windowInitializer);
        }
        for (int i=0; i < windowSize; i++) {
            variantContextWindow.moveWindow(null);
            filter();
        }

        if ( writer != null )
            writer.close();
    }
}

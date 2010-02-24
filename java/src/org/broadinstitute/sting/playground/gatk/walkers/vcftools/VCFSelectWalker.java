package org.broadinstitute.sting.playground.gatk.walkers.vcftools;

import org.broadinstitute.sting.gatk.contexts.*;
import org.broadinstitute.sting.gatk.refdata.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.filters.ClusteredSnps;
import org.broadinstitute.sting.gatk.walkers.filters.VariantContextWindow;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.vcf.*;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.util.*;
import org.apache.commons.jexl.*;

/**
 * Selects variant calls for output from a user-supplied VCF file using a number of user-selectable, parameterizable criteria.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant",type= RodVCF.class))
public class VCFSelectWalker extends RodWalker<Integer, Integer> {
    @Argument(fullName="match", shortName="match", doc="Expression used with INFO fields to select VCF records for inclusion in the output VCF(see wiki docs for more info)", required=false)
    protected String[] MATCH_STRINGS = new String[]{null};

    private VCFWriter writer = null;

    class MatchExp {
        String name;
        String expStr;
        Expression exp;

        public MatchExp(String name, String str, Expression exp) {
            this.name = name;
            this.expStr = str;
            this.exp = exp;
        }
    }

    private List<MatchExp> matchExpressions = new ArrayList<MatchExp>();

    private void initializeVcfWriter(RodVCF rod) {
        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "VariantSelect"));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));

        for ( MatchExp exp : matchExpressions ) {
            hInfo.add(new VCFFilterHeaderLine(exp.name, exp.expStr));
        }

        writer = new VCFWriter(out);
        writer.writeHeader(new VCFHeader(hInfo, rod.getHeader().getGenotypeSamples()));
    }

    public void initialize() {
        for ( int i = 0; i < MATCH_STRINGS.length; i++ ) {
            if ( MATCH_STRINGS[i] != null )  {
                try {
                    Expression filterExpression = ExpressionFactory.createExpression(MATCH_STRINGS[i]);
                    matchExpressions.add(new MatchExp(String.format("match-%d", i), MATCH_STRINGS[i], filterExpression));
                } catch (Exception e) {
                    throw new StingException("Invalid expression used (" + MATCH_STRINGS[i] + "). Please see the JEXL docs for correct syntax.");
                }
            }
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

        RODRecordList rods = tracker.getTrackData("variant", null);
        // ignore places where we don't have a variant
        if ( rods == null || rods.getRecords().size() == 0 )
            return 0;

        RodVCF variant = (RodVCF)rods.getRecords().get(0);
        boolean someoneMatched = false;
        for ( MatchExp exp : matchExpressions ) {
            Map<String, String> infoMap = new HashMap<String, String>(variant.mCurrentRecord.getInfoValues());
            infoMap.put("QUAL", String.valueOf(variant.mCurrentRecord.getQual()));

            JexlContext jContext = JexlHelper.createContext();
            jContext.setVars(infoMap);

            try {
                //System.out.printf("Matching %s vs. %s%n", infoMap, exp.expStr);
                if ( (Boolean)exp.exp.evaluate(jContext) ) {
                    //System.out.printf("  => Matched%n");
                    someoneMatched = true;
                    break;
                }
            } catch (Exception e) {
                throw new StingException(e.getMessage());
            }
        }

        if ( someoneMatched )
            writeVCF(variant);

        return 1;
    }

    private void writeVCF(RodVCF variant) {
        VCFRecord rec = variant.mCurrentRecord;

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
        if ( writer != null )
            writer.close();
    }
}
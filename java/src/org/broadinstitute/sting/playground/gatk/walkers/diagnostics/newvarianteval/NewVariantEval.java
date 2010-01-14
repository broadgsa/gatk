package org.broadinstitute.sting.playground.gatk.walkers.diagnostics.newvarianteval;

import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RODRecordList;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.PackageUtils;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.VCFRecord;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeEncoding;
import org.broadinstitute.sting.utils.genotype.vcf.VCFGenotypeRecord;
import org.broadinstitute.sting.playground.gatk.walkers.varianteval.VariantAnalysis;
import org.apache.commons.jexl.ExpressionFactory;
import org.apache.commons.jexl.Expression;
import org.apache.commons.jexl.JexlHelper;
import org.apache.commons.jexl.JexlContext;

import java.util.*;

@Requires(value={},referenceMetaData=@RMD(name="eval",type=RodVCF.class))
public class NewVariantEval extends RodWalker<Integer, Integer> {
    @Argument(fullName="filterExpression", shortName="filter", doc="Expression used with INFO fields to filter (see wiki docs for more info)", required=false)
    private String FILTER_STRING = "FILTER == false";

    @Argument(fullName="sampleName", shortName="sample", doc="If the VCF file has multiple samples, evaluate only the specified sample", required=false)
    private String SAMPLE_NAME = null;

    private Expression filterExpression;
    private HashSet<VariantEvaluation> evals;

    public void initialize() {
        try {
            filterExpression = ExpressionFactory.createExpression(FILTER_STRING);
        } catch (Exception e) {
            throw new StingException("Invalid expression used (" + FILTER_STRING + "). Please see the JEXL docs for correct syntax.");
        }

        try {
            evals = new HashSet<VariantEvaluation>();
            
            List<Class<? extends VariantEvaluation>> cevals = PackageUtils.getClassesImplementingInterface(VariantEvaluation.class);
            for (Class ceval : cevals) {
                out.printf("Analysis: %s%n", ceval.getName());

                VariantEvaluation eval = (VariantEvaluation) ceval.newInstance();
                evals.add(eval);
            }
        } catch (InstantiationException e) {
            throw new StingException(e.getMessage());
        } catch (IllegalAccessException e) {
            throw new StingException(e.getMessage());
        }
    }

    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData("eval", null);
        
        if ( rods != null ) {
            RodVCF variant = (RodVCF) rods.getRecords().get(0);

            Map<String, String> infoMap = new HashMap<String, String>(variant.mCurrentRecord.getInfoValues());
            infoMap.put("QUAL", String.valueOf(variant.mCurrentRecord.getQual()));
            infoMap.put("FILTER", String.valueOf(variant.mCurrentRecord.isFiltered()));
            infoMap.put("ID", String.valueOf(variant.mCurrentRecord.getID()));

            for (String filterCode : variant.mCurrentRecord.getFilteringCodes()) {
                infoMap.put(filterCode, "1");
            }

            JexlContext jContext = JexlHelper.createContext();
            jContext.setVars(infoMap);

            try {
                return (Boolean) filterExpression.evaluate(jContext);
            } catch (Exception e) {
                throw new StingException(e.getMessage());
            }
        }

        return true;
    }

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        RODRecordList<ReferenceOrderedDatum> rods = tracker.getTrackData("eval", null);
        RodVCF variant = (rods != null) ? (RodVCF) rods.getRecords().get(0) : null;

        if (variant != null && SAMPLE_NAME != null && variant.hasGenotypeData()) {
            variant = selectSample(SAMPLE_NAME, variant);

            if (variant == null) {
                throw new StingException(String.format("The sample '%s' is not present in the specified VCF", SAMPLE_NAME));
            }
        }

        for (VariantEvaluation eval : evals) {
            eval.update(variant, tracker, ref, context);
        }

        return 1;
    }

    private RodVCF selectSample(String sample, RodVCF variant) {
        String[] samples = variant.getSampleNames();

        for (int i = 0; i < samples.length; i++) {
            if (samples[i].equalsIgnoreCase(sample)) {
                List<VCFGenotypeRecord> genotypeRecs = new ArrayList<VCFGenotypeRecord>();
                genotypeRecs.add(variant.mCurrentRecord.getVCFGenotypeRecords().get(i));

                variant.mCurrentRecord = new VCFRecord(variant.getReferenceForSNP(),
                                                       variant.getLocation().getContig(),
                                                       (int) variant.getLocation().getStart(),
                                                       variant.getID(),
                                                       variant.mCurrentRecord.getAlternateAlleles(),
                                                       variant.getQual(),
                                                       variant.getFilterString(),
                                                       variant.getInfoValues(),
                                                       variant.mCurrentRecord.getGenotypeFormatString(),
                                                       genotypeRecs);

                return variant;
            }
        }

        return null;
    }

    public Integer reduceInit() {
        return null;
    }

    public Integer reduce(Integer value, Integer sum) {
        return null;
    }

    public void onTraversalDone(Integer sum) {
        for (VariantEvaluation eval : evals) {
            List<String> results = eval.getResult();

            for (String result : results) {
                out.printf("%s:%s%n", eval.getName(), result);
            }
        }
    }
}

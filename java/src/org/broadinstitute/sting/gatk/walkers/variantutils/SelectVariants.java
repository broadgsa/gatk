package org.broadinstitute.sting.gatk.walkers.variantutils;

import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFHeader;
import org.broad.tribble.vcf.VCFHeaderLine;
import org.broad.tribble.vcf.VCFHeaderLineType;
import org.broad.tribble.vcf.VCFInfoHeaderLine;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Takes a VCF file, selects variants based on sample(s) in which it was found and/or on various annotation criteria,
 * recompute the value of certain annotations based on the new sample set, and output a new VCF with the results.
 */
@Requires(value={},referenceMetaData=@RMD(name="variant", type=ReferenceOrderedDatum.class))
public class SelectVariants extends RodWalker<Integer, Integer> {
    @Argument(fullName="sample", shortName="sn", doc="Sample(s) to include.  Can be a single sample, specified multiple times for many samples, or a regular expression to select many samples.", required=false)
    public Set<String> SAMPLE_EXPRESSIONS;

    @Argument(shortName="select", doc="One or more criteria to use when selecting the data.  Evaluated *after* the specified samples are extracted and the INFO-field annotations are updated.", required=false)
    public ArrayList<String> SELECT_EXPRESSIONS = new ArrayList<String>();

    @Argument(fullName="excludeNonVariants", shortName="env", doc="Don't include loci found to be non-variant after the subsetting procedure.", required=false)
    private boolean EXCLUDE_NON_VARIANTS = false;

    @Argument(fullName="excludeFiltered", shortName="ef", doc="Don't include filtered loci.", required=false)
    private boolean EXCLUDE_FILTERED = false;

    private ArrayList<String> selectNames = new ArrayList<String>();
    private List<VariantContextUtils.JexlVCMatchExp> jexls = null;

    private Set<String> samples = new HashSet<String>();
    private Set<String> possibleSampleRegexs = new HashSet<String>();
    private Set<String> sampleExpressionsThatDidNotWork = new HashSet<String>();

    private VCFWriter vcfWriter = null;

    /**
     * Set up the VCF writer, the sample expressions and regexs, and the JEXL matcher
     */
    public void initialize() {
        vcfWriter = new VCFWriter(out);

        ArrayList<String> rodNames = new ArrayList<String>();
        rodNames.add("variant");

        Map<String, VCFHeader> vcfRods = VCFUtils.getVCFHeadersFromRods(getToolkit(), rodNames);
        Set<String> vcfSamples = SampleUtils.getSampleList(vcfRods, VariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

        if (SAMPLE_EXPRESSIONS != null) {
            // Let's first assume that the values in SAMPLE_EXPRESSIONS are literal sample names and not regular
            // expressions.  Extract those samples specifically so we don't make the mistake of selecting more
            // than what the user really wants.
            for (String SAMPLE_EXPRESSION : SAMPLE_EXPRESSIONS) {
                if (vcfSamples.contains(SAMPLE_EXPRESSION)) {
                    samples.add(SAMPLE_EXPRESSION);
                } else {
                    possibleSampleRegexs.add(SAMPLE_EXPRESSION);
                }
            }

            // Now, check the expressions that weren't used in the previous step, and use them as if they're regular expressions
            for (String sampleRegex : possibleSampleRegexs) {
                Pattern p = Pattern.compile(sampleRegex);

                boolean patternWorked = false;

                for (String vcfSample : vcfSamples) {
                    Matcher m = p.matcher(vcfSample);
                    if (m.find()) {
                        samples.add(vcfSample);

                        patternWorked = true;
                    }
                }

                if (!patternWorked) {
                    sampleExpressionsThatDidNotWork.add(sampleRegex);
                }
            }

            // Finally, warn the user about any leftover sample expressions that had no effect
            if (sampleExpressionsThatDidNotWork.size() > 0) {
                for (String exp : sampleExpressionsThatDidNotWork) {
                    logger.warn("The sample expression '" + exp + "' had no effect (no matching sample or pattern match found).  Skipping.");
                }
            }
        } else {
            samples.addAll(vcfSamples);
        }

        for (String sample : samples) {
            logger.info("Including sample '" + sample + "'");
        }

        Set<VCFHeaderLine> headerLines = VCFUtils.smartMergeHeaders(vcfRods.values(), logger);
        headerLines.add(new VCFHeaderLine("source", "SelectVariants"));
        vcfWriter.writeHeader(new VCFHeader(headerLines, samples));

        for (int i = 0; i < SELECT_EXPRESSIONS.size(); i++) {
            // It's not necessary that the user supply select names for the JEXL expressions, since those
            // expressions will only be needed for omitting records.  Make up the select names here.
            selectNames.add(String.format("select-%d", i));
        }

        jexls = VariantContextUtils.initializeMatchExps(selectNames, SELECT_EXPRESSIONS);
    }

//    /**
//     * If JEXL expressions are supplied, include only records that satisfy the expression
//     *
//     * @param  tracker   the ROD tracker
//     * @param  ref       reference information
//     * @param  context   alignment info
//     * @return true if no JEXL expressions are supplied or if a record satisfies all JEXL criteria, false if otherwise
//     */
//    public boolean filter(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
//        VariantContext vc = tracker.getVariantContext(ref, "variant", null, context.getLocation(), true);
//
//        for ( VariantContextUtils.JexlVCMatchExp jexl : jexls ) {
//            if ( !VariantContextUtils.match(vc, jexl) ) {
//                return false;
//            }
//        }
//
//        return true;
//    }

    /**
     * Subset VC record if necessary and emit the modified record (provided it satisfies criteria for printing)
     *
     * @param  tracker   the ROD tracker
     * @param  ref       reference information
     * @param  context   alignment info
     * @return 1 if the record was printed to the output file, 0 if otherwise
     */
    @Override
    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return 0;

        VariantContext vc = tracker.getVariantContext(ref, "variant", null, context.getLocation(), true);
        VariantContext sub = subsetRecord(vc, samples);

        if ( (sub.isPolymorphic() || !EXCLUDE_NON_VARIANTS) && (!sub.isFiltered() || !EXCLUDE_FILTERED) ) {
            for ( VariantContextUtils.JexlVCMatchExp jexl : jexls ) {
                if ( !VariantContextUtils.match(sub, jexl) ) {
                    return 0;
                }
            }

            vcfWriter.add(sub, ref.getBase());
        }

        return 1;
    }

    /**
     * Helper method to subset a VC record, modifying some metadata stored in the INFO field (i.e. AN, AC, AF).
     *
     * @param vc       the VariantContext record to subset
     * @param samples  the samples to extract
     * @return the subsetted VariantContext
     */
    private VariantContext subsetRecord(VariantContext vc, Set<String> samples) {
        if ( samples == null || samples.isEmpty() )
            return vc;

        ArrayList<Genotype> genotypes = new ArrayList<Genotype>();
        for ( Map.Entry<String, Genotype> genotypePair : vc.getGenotypes().entrySet() ) {
            if ( samples.contains(genotypePair.getKey()) )
                genotypes.add(genotypePair.getValue());
        }
        
        VariantContext sub = vc.subContextFromGenotypes(genotypes);

        HashMap<String, Object> attributes = new HashMap<String, Object>(sub.getAttributes());

        int alleleCount = 0;
        int numberOfAlleles = 0;
        int depth = 0;
        for (String sample : sub.getSampleNames()) {
            Genotype g = sub.getGenotype(sample);

            if (g.isNotFiltered() && g.isCalled()) {
                numberOfAlleles += g.getPloidy();

                if (g.isHet()) { alleleCount++; }
                else if (g.isHomVar()) { alleleCount += 2; }
                
                String dp = (String) g.getAttribute("DP");
                if (dp != null) {
                    depth += Integer.valueOf(dp);
                }
            }
        }

        attributes.put("AC", alleleCount);
        attributes.put("AN", numberOfAlleles);
        if (numberOfAlleles == 0) {
            attributes.put("AF", 0.0);
        } else {
            attributes.put("AF", ((double) alleleCount) / ((double) numberOfAlleles));
        }
        attributes.put("DP", depth);

        sub = VariantContextUtils.modifyAttributes(sub, attributes);

        return sub;
    }

    @Override
    public Integer reduceInit() { return 0; }

    @Override
    public Integer reduce(Integer value, Integer sum) { return value + sum; }

    public void onTraversalDone(Integer result) { logger.info(result + " records processed."); }
}

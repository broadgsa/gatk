/*
 * Copyright (c) 2010.
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

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.vcf.*;
import org.broadinstitute.sting.utils.variantcontext.Allele;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.utils.variantcontext.VariantContextUtils;

import java.util.*;

/**
 * Annotates a validation (from Sequenom for example) VCF with QC metrics (HW-equilibrium, % failed probes)
 *
 * <p>
 * The Variant Validation Assessor is a tool for vetting/assessing validation data (containing genotypes).
 * The tool produces a VCF that is annotated with information pertaining to plate quality control and by
 * default is soft-filtered by high no-call rate or low Hardy-Weinberg probability.
 * If you have .ped files, please first convert them to VCF format
 * (see http://www.broadinstitute.org/gsa/wiki/index.php/Converting_ped_to_vcf).
 *
 * <h2>Input</h2>
 * <p>
 * A validation VCF to annotate.
 * </p>
 *
 * <h2>Output</h2>
 * <p>
 * An annotated VCF.  Additionally, a table like the following will be output:
 * <pre>
 *     Total number of samples assayed:                  185
 *     Total number of records processed:                152
 *     Number of Hardy-Weinberg violations:              34 (22%)
 *     Number of no-call violations:                     12 (7%)
 *     Number of homozygous variant violations:          0 (0%)
 *     Number of records passing all filters:            106 (69%)
 *     Number of passing records that are polymorphic:   98 (92%)
 * </pre>
 * </p>
 *
 * <h2>Examples</h2>
 * <pre>
 * java -Xmx2g -jar GenomeAnalysisTK.jar \
 *   -R ref.fasta \
 *   -T VariantValidationAssessor \
 *   --variant input.vcf \
 *   -o output.vcf
 * </pre>
 *
 */
@Reference(window=@Window(start=0,stop=40))
public class VariantValidationAssessor extends RodWalker<VariantContext,Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    @Output(doc="File to which variants should be written",required=true)
    protected VCFWriter vcfwriter = null;

    @Argument(fullName="maxHardy", doc="Maximum phred-scaled Hardy-Weinberg violation pvalue to consider an assay valid", required=false)
    protected double maxHardy = 20.0;

    /**
     * To disable, set to a value greater than 1.
     */
    @Argument(fullName="maxNoCall", doc="Maximum no-call rate (as a fraction) to consider an assay valid", required=false)
    protected double maxNoCall = 0.05;

    /**
     * To disable, set to a value greater than 1.
     */
    @Argument(fullName="maxHomVar", doc="Maximum homozygous variant rate (as a fraction) to consider an assay valid", required=false)
    protected double maxHomNonref = 1.1;

    //@Argument(fullName="populationFile", shortName="populations", doc="A tab-delimited file relating individuals to populations,"+
    //          "used for smart Hardy-Weinberg annotation",required = false)
    //private File popFile = null;

    // sample names
    private TreeSet<String> sampleNames = null;

    // variant context records
    private ArrayList<VariantContext> records = new ArrayList<VariantContext>();

    // statistics
    private int numRecords = 0;
    private int numHWViolations = 0;
    private int numNoCallViolations = 0;
    private int numHomVarViolations = 0;
    private int numTrueVariants = 0;

    //private HashMap<String,String> samplesToPopulation;

    public void initialize() {
        //if ( popFile != null ) {
        //    samplesToPopulation = parsePopulationFile(popFile);
        //}
    }

    public Integer reduceInit() {
        return 0;
    }

    public VariantContext map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        VariantContext vc = tracker.getFirstValue(variantCollection.variants, ref.getLocus());
        // ignore places where we don't have a variant
        if ( vc == null )
            return null;

        if ( sampleNames == null )
            sampleNames = new TreeSet<String>(vc.getSampleNames());        

        return addVariantInformationToCall(vc);
    }

    public Integer reduce(VariantContext call, Integer numVariants) {
        if ( call != null ) {
            numVariants++;
            records.add(call);
        }
        return numVariants;                        
    }

    public void onTraversalDone(Integer finalReduce) {
        final List<String> inputNames = Arrays.asList(variantCollection.variants.getName());

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), inputNames));

        // set up the info and filter headers
        hInfo.add(new VCFInfoHeaderLine("NoCallPct", 1, VCFHeaderLineType.Float, "Percent of no-calls"));
        hInfo.add(new VCFInfoHeaderLine("HomRefPct", 1, VCFHeaderLineType.Float, "Percent of homozygous reference genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HetPct", 1, VCFHeaderLineType.Float, "Percent of heterozygous genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HomVarPct", 1, VCFHeaderLineType.Float, "Percent homozygous variant genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HW", 1, VCFHeaderLineType.Float, "Phred-scaled Hardy-Weinberg violation p-value"));
        hInfo.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_COUNT_KEY, 1, VCFHeaderLineType.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
        hInfo.add(new VCFInfoHeaderLine(VCFConstants.ALLELE_NUMBER_KEY, 1, VCFHeaderLineType.Integer, "Total number of alleles in called genotypes"));
        hInfo.add(new VCFFilterHeaderLine("HardyWeinbergViolation", "The validation is in Hardy-Weinberg violation"));
        hInfo.add(new VCFFilterHeaderLine("HighNoCallRate", "The validation no-call rate is too high"));
        hInfo.add(new VCFFilterHeaderLine("TooManyHomVars", "The validation homozygous variant rate is too high"));

        // print out (and add to headers) the validation metrics
        System.out.println(String.format("Total number of samples assayed:\t\t\t%d", sampleNames.size()));
        hInfo.add(new VCFHeaderLine("ValidationMetrics_SamplesAssayed", String.format("%d", sampleNames.size())));
        System.out.println(String.format("Total number of records processed:\t\t\t%d", numRecords));
        hInfo.add(new VCFHeaderLine("ValidationMetrics_RecordsProcessed", String.format("%d", numRecords)));
        if ( numRecords > 0 ) {
            System.out.println(String.format("Number of Hardy-Weinberg violations:\t\t\t%d (%d%%)", numHWViolations, 100*numHWViolations/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_HardyWeinbergViolations", String.format("\"%d (%d%%)\"", numHWViolations, 100*numHWViolations/numRecords)));
            System.out.println(String.format("Number of no-call violations:\t\t\t\t%d (%d%%)", numNoCallViolations, 100*numNoCallViolations/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_NoCallViolations", String.format("\"%d (%d%%)\"", numNoCallViolations, 100*numNoCallViolations/numRecords)));
            System.out.println(String.format("Number of homozygous variant violations:\t\t%d (%d%%)", numHomVarViolations, 100*numHomVarViolations/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_HomVarViolations", String.format("\"%d (%d%%)\"", numHomVarViolations, 100*numHomVarViolations/numRecords)));
            int goodRecords = numRecords - numHWViolations - numNoCallViolations - numHomVarViolations;
            System.out.println(String.format("Number of records passing all filters:\t\t\t%d (%d%%)", goodRecords, 100*goodRecords/numRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_RecordsPassingFilters", String.format("\"%d (%d%%)\"", goodRecords, 100*goodRecords/numRecords)));
            if ( goodRecords > 0 ) {
                System.out.println(String.format("Number of passing records that are polymorphic:\t\t%d (%d%%)", numTrueVariants, 100*numTrueVariants/goodRecords));
                hInfo.add(new VCFHeaderLine("ValidationMetrics_PolymorphicPassingRecords", String.format("\"%d (%d%%)\"", numTrueVariants, 100*numTrueVariants/goodRecords)));
            }
        }
        
        vcfwriter.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));

        for ( VariantContext record : records )
            vcfwriter.add(record);
    }


    private VariantContext addVariantInformationToCall(VariantContext vContext) {

        // check possible filters
        double hwPvalue = hardyWeinbergCalculation(vContext);
        double hwScore = Math.abs(QualityUtils.phredScaleErrorRate(hwPvalue));
        double noCallProp = (double)vContext.getNoCallCount() / (double)vContext.getNSamples();
        double homRefProp = (double)vContext.getHomRefCount() / (double)vContext.getNSamples();
        double hetProp = (double)vContext.getHetCount() / (double)vContext.getNSamples();
        double homVarProp = (double)vContext.getHomVarCount() / (double)vContext.getNSamples();

        boolean isViolation = false;
        Set<String> filters = new HashSet<String>();
        if ( noCallProp > maxNoCall ) {
            filters.add("HighNoCallRate");
            numNoCallViolations++;
            isViolation = true;
        } else if ( hwScore > maxHardy ) {
            filters.add("HardyWeinbergViolation");
            numHWViolations++;
            isViolation = true;
        } else if ( homVarProp > maxHomNonref) {
            filters.add("TooManyHomVars");
            numHomVarViolations++;
            isViolation = true;
        }
        vContext = VariantContext.modifyFilters(vContext, filters);
        numRecords++;

        // add the info fields
        HashMap<String, Object> infoMap = new HashMap<String, Object>();
        infoMap.put("NoCallPct", String.format("%.1f", 100.0*noCallProp));
        infoMap.put("HomRefPct", String.format("%.1f", 100.0*homRefProp));
        infoMap.put("HomVarPct", String.format("%.1f", 100.0*homVarProp));
        infoMap.put("HetPct", String.format("%.1f", 100.0*hetProp));
        infoMap.put("HW", String.format("%.2f", hwScore));
        Collection<Allele> altAlleles = vContext.getAlternateAlleles();
        int altAlleleCount = altAlleles.size() == 0 ? 0 : vContext.getChromosomeCount(altAlleles.iterator().next());
        if ( !isViolation && altAlleleCount > 0 )
            numTrueVariants++;
        infoMap.put(VCFConstants.ALLELE_COUNT_KEY, String.format("%d", altAlleleCount));
        infoMap.put(VCFConstants.ALLELE_NUMBER_KEY, String.format("%d", vContext.getChromosomeCount()));

        return VariantContext.modifyAttributes(vContext, infoMap);
    }

    private double hardyWeinbergCalculation(VariantContext vc) {
        //if ( popFile != null ) {
        //    throw new StingException("We still need to implement this!");
        //} else {
        return VariantContextUtils.computeHardyWeinbergPvalue(vc);
        //}
    }

    // TODO -- REWRITE THIS TO WORK WITH VARIANT CONTEXT
    /******

    private String smartHardy(ReferenceContext ref, VCFRecord rec) {
        HashMap<String,ArrayList<Genotype>> genotypesByPopulation = new HashMap<String,ArrayList<Genotype>>(10);
        HashMap<String,String> hardyWeinbergByPopulation = new HashMap<String,String>(10);

        for ( String population : samplesToPopulation.values() ) {
            genotypesByPopulation.put(population,new ArrayList<Genotype>());
        }

        //for ( String name : sampleNames ) {
        //    String pop = samplesToPopulation.get(name);
        //    if ( rec.getGenotype(name) != null ) {
        //        genotypesByPopulation.get(pop).add(rec.getGenotype(name));
        //    }
        //}

        for ( String population : samplesToPopulation.values() ) {
            VCFVariationCall v = new VCFVariationCall(ref.getBase(),ref.getLocus(),VCFVariationCall.VARIANT_TYPE.SNP);
            v.setGenotypeCalls(genotypesByPopulation.get(population));
            hardyWeinbergByPopulation.put(population,HWCalc.annotate(null,ref,null,v));
        }

        return smartHardyString(hardyWeinbergByPopulation);
    }

    private String smartHardyString(HashMap<String,String> hwByPop) {
        // for now just return the maximum:
        int maxH = -100;
        for ( String pop : samplesToPopulation.values() ) {
            maxH = Integer.parseInt(hwByPop.get(pop)) > maxH ? Integer.parseInt(hwByPop.get(pop)) : maxH;
        }

        return String.format("%s",maxH);
    }

    *********/
}

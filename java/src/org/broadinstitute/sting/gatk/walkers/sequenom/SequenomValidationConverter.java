package org.broadinstitute.sting.gatk.walkers.sequenom;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.refdata.PlinkRod;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.util.*;

/**
 * Converts Sequenom files to a VCF annotated with QC metrics (HW-equilibrium, % failed probes)
 */
@Reference(window=@Window(start=0,stop=40))
@Requires(value={},referenceMetaData=@RMD(name="sequenom",type= ReferenceOrderedDatum.class))
public class SequenomValidationConverter extends RodWalker<VCFRecord,Integer> {
    @Argument(fullName="maxHardy", doc="Maximum phred-scaled Hardy-Weinberg violation pvalue to consider an assay valid [default:20]", required=false)
    protected double maxHardy = 20.0;
    @Argument(fullName="maxNoCall", doc="Maximum no-call rate (as a fraction) to consider an assay valid [default:0.05]", required=false)
    protected double maxNoCall = 0.05;
    @Argument(fullName="maxHomVar", doc="Maximum homozygous variant rate (as a fraction) to consider an assay valid [default:1.1, disabled]", required=false)
    protected double maxHomNonref = 1.1;

    //@Argument(fullName="populationFile", shortName="populations", doc="A tab-delimited file relating individuals to populations,"+
    //          "used for smart Hardy-Weinberg annotation",required = false)
    //private File popFile = null;

    // max allowable indel size (based on ref window)
    private static final int MAX_INDEL_SIZE = 40;

    // sample names
    private TreeSet<String> sampleNames = null;

    // vcf records
    private ArrayList<VCFRecord> records = new ArrayList<VCFRecord>();

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
        int numberOfVariantsProcessed = 0;
        return numberOfVariantsProcessed;
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        // get the sequenom rod at this locus if there is one
        List<Object> rods = tracker.getReferenceMetaData("sequenom");
        // ignore places where we don't have a variant
        if ( rods.size() == 0 )
            return null;

        Object rod = rods.get(0);

        // determine the reference allele
        Allele refAllele = determineRefAllele(rod, ref);

        VariantContext vc = VariantContextAdaptors.toVariantContext("sequenom", rod, refAllele);

        if ( sampleNames == null )
            sampleNames = new TreeSet<String>(vc.getSampleNames());        

        return addVariantInformationToCall(ref, vc, rod);
    }

    private Allele determineRefAllele(Object rod, ReferenceContext ref) {
        Allele refAllele;

        // ugly hack to get around the fact that the Plink rod needs
        // a very specific determination of the reference allele
        if ( rod instanceof PlinkRod ) {
            PlinkRod plink = (PlinkRod)rod;
            if ( !plink.isIndel() ) {
                refAllele = new Allele(Character.toString(ref.getBase()), true);
            } else if ( plink.isInsertion() ) {
                refAllele = new Allele(PlinkRod.SEQUENOM_NO_BASE, true);
            } else {
                if ( plink.getLength() > MAX_INDEL_SIZE )
                    throw new UnsupportedOperationException("PlinkToVCF currently can only handle indels up to length " + MAX_INDEL_SIZE);
                char[] deletion = new char[plink.getLength()];
                System.arraycopy(ref.getBases(), 1, deletion, 0, plink.getLength());
                refAllele = new Allele(new String(deletion), true);
            }
        } else {
            refAllele = new Allele(Character.toString(ref.getBase()), true);
        }

        return refAllele;
    }


    public Integer reduce(VCFRecord call, Integer numVariants) {
        if ( call != null ) {
            numVariants++;
            records.add(call);
        }
        return numVariants;                        
    }

    public void onTraversalDone(Integer finalReduce) {
        if ( sampleNames == null )
            sampleNames = new TreeSet<String>();

        VCFWriter vcfWriter = new VCFWriter(out);

        // set up the info and filter headers
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.add(new VCFHeaderLine("source", "SequenomValidationConverter"));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        hInfo.add(new VCFInfoHeaderLine("NoCallPct", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Percent of no-calls"));
        hInfo.add(new VCFInfoHeaderLine("HomRefPct", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Percent of homozygous reference genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HetPct", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Percent of heterozygous genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HomVarPct", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Percent homozygous variant genotypes"));
        hInfo.add(new VCFInfoHeaderLine("HW", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "Phred-scaled Hardy-Weinberg violation p-value"));
        hInfo.add(new VCFInfoHeaderLine(VCFRecord.ALLELE_COUNT_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Allele count in genotypes, for each ALT allele, in the same order as listed"));
        hInfo.add(new VCFInfoHeaderLine(VCFRecord.ALLELE_NUMBER_KEY, 1, VCFInfoHeaderLine.INFO_TYPE.Integer, "Total number of alleles in called genotypes"));
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
            System.out.println(String.format("Number of passing records that are polymorphic:\t\t%d (%d%%)", numTrueVariants, 100*numTrueVariants/goodRecords));
            hInfo.add(new VCFHeaderLine("ValidationMetrics_PolymorphicPassingRecords", String.format("\"%d (%d%%)\"", numTrueVariants, 100*numTrueVariants/goodRecords)));
        }
        
        VCFHeader header = new VCFHeader(hInfo, sampleNames);
        vcfWriter.writeHeader(header);

        for ( VCFRecord record : records )
            vcfWriter.addRecord(record);
        vcfWriter.close();
    }


    private VCFRecord addVariantInformationToCall(ReferenceContext ref, VariantContext vContext, Object rod) {

        VCFRecord record = VariantContextAdaptors.toVCF(vContext, ref.getBase());
        record.setGenotypeFormatString("GT");

        // check possible filters
        double hwPvalue = hardyWeinbergCalculation(vContext);
        double hwScore = Math.abs(QualityUtils.phredScaleErrorRate(hwPvalue));
        double noCallProp = (double)vContext.getNoCallCount() / (double)vContext.getNSamples();
        double homRefProp = (double)vContext.getHomRefCount() / (double)vContext.getNSamples();
        double hetProp = (double)vContext.getHetCount() / (double)vContext.getNSamples();
        double homVarProp = (double)vContext.getHomVarCount() / (double)vContext.getNSamples();

        boolean isViolation = false;
        if ( noCallProp > maxNoCall ) {
            record.setFilterString("HighNoCallRate");
            numNoCallViolations++;
            isViolation = true;
        } else if ( hwScore > maxHardy ) {
            record.setFilterString("HardyWeinbergViolation");
            numHWViolations++;
            isViolation = true;
        } else if ( homVarProp > maxHomNonref) {
            record.setFilterString("TooManyHomVars");
            numHomVarViolations++;
            isViolation = true;
        }
        numRecords++;

        // add the info fields
        HashMap<String, String> infoMap = new HashMap<String,String>(5);
        infoMap.put("NoCallPct", String.format("%.1f", 100.0*noCallProp));
        infoMap.put("HomRefPct", String.format("%.1f", 100.0*homRefProp));
        infoMap.put("HomVarPct", String.format("%.1f", 100.0*homVarProp));
        infoMap.put("HetPct", String.format("%.1f", 100.0*hetProp));
        infoMap.put("HW", String.format("%.2f", hwScore));
        Set<Allele> altAlleles = vContext.getAlternateAlleles();
        int altAlleleCount = altAlleles.size() == 0 ? 0 : vContext.getChromosomeCount(altAlleles.iterator().next());
        if ( !isViolation && altAlleleCount > 0 )
            numTrueVariants++;
        infoMap.put(VCFRecord.ALLELE_COUNT_KEY, String.format("%d", altAlleleCount));
        infoMap.put(VCFRecord.ALLELE_NUMBER_KEY, String.format("%d", vContext.getChromosomeCount()));
        record.addInfoFields(infoMap);

        // set the id if it's a plink rod
        if ( rod instanceof PlinkRod )
            record.setID(((PlinkRod)rod).getVariantName());

        return record;
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

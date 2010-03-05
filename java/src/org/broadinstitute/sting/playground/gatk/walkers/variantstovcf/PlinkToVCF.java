package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.QualityUtils;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.PlinkRod;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.Reference;
import org.broadinstitute.sting.gatk.walkers.Window;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Converts Sequenom files to a VCF annotated with QC metrics (HW-equilibrium, % failed probes)
 */
@Reference(window=@Window(start=0,stop=40))
public class PlinkToVCF extends RodWalker<VCFRecord,Integer> {
    @Argument(fullName="outputVCF", shortName="vcf", doc="The VCF file to write results", required=true)
    protected File vcfFile = null;
    @Argument(fullName="maxHardy", doc="Maximum phred-scaled Hardy-Weinberg violation pvalue to consider an assay valid [default:20]", required=false)
    protected double maxHardy = 20.0;
    @Argument(fullName="maxNoCall", doc="Maximum no-call rate (as a fraction) to consider an assay valid [default:0.05]", required=false)
    protected double maxNoCall = 0.05;
    @Argument(fullName="maxHomVar", doc="Maximum homozygous variant rate (as a fraction) to consider an assay valid [default:1.1, disabled]", required=false)
    protected double maxHomNonref = 1.1;

    @Argument(fullName="populationFile", shortName="populations", doc="A tab-delimited file relating individuals to populations,"+
              "used for smart Hardy-Weinberg annotation",required = false)
    public File popFile = null;

    // max allowable indel size (based on ref window)
    private static final int MAX_INDEL_SIZE = 40;

    // the VCF writer
    private VCFWriter vcfWriter = null;

    // statistics
    private int numRecords = 0;
    private int numHWViolations = 0;
    private int numNoCallViolations = 0;
    private int numHomVarViolations = 0;
    private int numTrueVariants = 0;

    private HashMap<String,String> samplesToPopulation;

    public void initialize() {
        if ( popFile != null ) {
            samplesToPopulation = parsePopulationFile(popFile);
        }
    }

    public Integer reduceInit() {
        int numberOfVariantsProcessed = 0;
        return numberOfVariantsProcessed;
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( tracker == null )
            return null;

        // get the Plink rod at this locus if there is one
        PlinkRod plinkRod = null;
        Iterator<ReferenceOrderedDatum> rods = tracker.getAllRods().iterator();
        while (rods.hasNext()) {
            ReferenceOrderedDatum rod = rods.next();
            if ( rod instanceof PlinkRod ) {
                plinkRod = (PlinkRod)rod;
                break;
            }
        }

        if ( plinkRod == null )
            return null;

        if ( vcfWriter == null )
            initializeWriter(plinkRod);

        return addVariantInformationToCall(ref, plinkRod);
    }

    private void initializeWriter(PlinkRod plinkRod) {
        vcfWriter = new VCFWriter(vcfFile);

        // set up the info and filter headers
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.add(new VCFHeaderLine("source", "PlinkToVCF"));
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

        VCFHeader header = new VCFHeader(hInfo, new TreeSet<String>(plinkRod.getSampleNames()));
        vcfWriter.writeHeader(header);
    }

    public Integer reduce(VCFRecord call, Integer numVariants) {
        if ( call != null ) {
            numVariants++;
            vcfWriter.addRecord(call);
        }
        return numVariants;                        
    }

    public void onTraversalDone(Integer finalReduce) {
        if ( vcfWriter != null )
            vcfWriter.close();
        System.out.println(String.format("Total number of records processed:\t\t\t%d", numRecords));
        System.out.println(String.format("Number of Hardy-Weinberg violations:\t\t\t%d (%d%%)", numHWViolations, 100*numHWViolations/numRecords));
        System.out.println(String.format("Number of no-call violations:\t\t\t\t%d (%d%%)", numNoCallViolations, 100*numNoCallViolations/numRecords));
        System.out.println(String.format("Number of homozygous variant violations:\t\t%d (%d%%)", numHomVarViolations, 100*numHomVarViolations/numRecords));
        int goodRecords = numRecords - numHWViolations - numNoCallViolations - numHomVarViolations;
        System.out.println(String.format("Number of records passing all filters:\t\t\t%d (%d%%)", goodRecords, 100*goodRecords/numRecords));
        System.out.println(String.format("Number of passing records that validated as true:\t%d (%d%%)", numTrueVariants, 100*numTrueVariants/goodRecords));
    }


    private VCFRecord addVariantInformationToCall(ReferenceContext ref, PlinkRod plinkRod) {

        // determine the reference allele
        Allele refAllele;
        if ( !plinkRod.isIndel() ) {
            refAllele = new Allele(Character.toString(ref.getBase()), true);
        } else if ( plinkRod.isInsertion() ) {
            refAllele = new Allele(PlinkRod.SEQUENOM_NO_BASE, true);
        } else {
            if ( plinkRod.getLength() > MAX_INDEL_SIZE )
                throw new UnsupportedOperationException("PlinkToVCF currently can only handle indels up to length " + MAX_INDEL_SIZE);
            char[] deletion = new char[plinkRod.getLength()];
            System.arraycopy(ref.getBases(), 1, deletion, 0, plinkRod.getLength());
            refAllele = new Allele(new String(deletion), true);
        }

        VariantContext vContext = VariantContextAdaptors.toVariantContext(plinkRod.getName(), plinkRod, refAllele);
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

        // add the id
        record.setID(plinkRod.getVariantName());

        return record;
    }

    private double hardyWeinbergCalculation(VariantContext vc) {
        if ( popFile != null ) {
            throw new StingException("We still need to implement this!");
        } else {
            return VariantContextUtils.computeHardyWeinbergPvalue(vc);
        }
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

    private HashMap<String,String> parsePopulationFile(File file) {
        HashMap<String,String> samplesToPopulation = new HashMap<String,String>();
        try {
            BufferedReader in = new BufferedReader( new FileReader( file ));
            String line = in.readLine();
            while ( line != null ) {
                String[] populationSamples = line.split("\t");
                String population = populationSamples[0];
                for ( int i = 1; i < populationSamples.length; i ++ ) {
                    samplesToPopulation.put(populationSamples[i], population);
                }
            }
        } catch ( IOException e) {
            throw new StingException("Error reading population file", e);
        }

        return samplesToPopulation;
    }
}

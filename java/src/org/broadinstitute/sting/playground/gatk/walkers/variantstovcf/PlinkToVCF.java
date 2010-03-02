package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.walkers.annotator.HardyWeinberg;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.PlinkRod;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
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
public class PlinkToVCF extends RodWalker<VCFRecord,Integer> {
    @Argument(fullName="outputVCF", shortName="vcf", doc="The VCF file to write to", required=true)
    public File vcfFile = null;
    @Argument(fullName="numberOfSamples", shortName="ns", doc="The number of samples sequenced", required=false)
    public int nSamples = 192; // The number of samples in the ped file I wrote this tool for
    @Argument(fullName="populationFile", shortName="populations", doc="A tab-delimited file relating individuals to populations,"+
              "used for smart Hardy-Weinberg annotation",required = false)
    public File popFile = null;
    @Argument(fullName="useb36ContigNames",shortName="b36contig",doc="Uses b36 contig names (1:1,000,000) rather than hg18 (chr1:1,000,000) for comparison with ref", required=false)
    public boolean useb36contigs=false;
    @Argument(fullName="maxHardy", doc="Maximum Hardy-Weinberg score to consider an assay valid", required=false)
    public double maxHardy = 10;
    @Argument(fullName="maxNoCall", doc="Maximum no-call rate (as a proportion) to consider an assay valid", required=false)
    public double maxNoCall = 0.05;
    @Argument(fullName="maxHomNonref", doc="Maximum homozygous-nonreference rate (as a proportion) to consider an assay valid", required = false)
    public double maxHomNonref = 1.1;

    private final Set<String> HEADER_FIELDS = new HashSet<String>(Arrays.asList("#Family ID","Individual ID","Sex","Paternal ID","Maternal ID","Phenotype",
                "FID","IID","PAT","MAT","SEX","PHENOTYPE"));
    private final int INIT_NUMBER_OF_POPULATIONS = 10;
    private final int DEFAULT_QUALITY = 20;
    private ArrayList<String> sampleNames = new ArrayList<String>(nSamples);
    private VCFGenotypeWriterAdapter vcfWriter;
    private final HardyWeinberg HWCalc = new HardyWeinberg();
    private final boolean useSmartHardy = popFile != null;
    private HashMap<String,String> samplesToPopulation;

    public void initialize() {
        vcfWriter = new VCFGenotypeWriterAdapter(vcfFile);
        if ( useSmartHardy ) {
            samplesToPopulation = parsePopulationFile(popFile);
        }
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "PlinkToVCF"));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        vcfWriter.writeHeader(new TreeSet<String>(sampleNames), hInfo);
        nSamples = sampleNames.size();
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

        return addVariantInformationToCall(ref, plinkRod);
    }

    public Integer reduce(VCFRecord call, Integer numVariants) {
        if ( call != null ) {
            numVariants++;
            printToVCF(call);
        }
        return numVariants;                        
    }

    public void onTraversalDone(Integer finalReduce) {
        logger.info("Variants processed="+finalReduce.toString());
        vcfWriter.close();
    }

    private void printToVCF(VCFRecord call) {
        try {
            vcfWriter.addRecord(call);
        } catch ( RuntimeException e ) {
            if ( e.getLocalizedMessage().equalsIgnoreCase("We have more genotype samples than the header specified")) {
                throw new StingException("We have more sample genotypes than sample names -- check that there are no duplicates in the .ped file",e);
            } else {
                throw new StingException("Error in VCF creation: "+e.getLocalizedMessage(),e);
            }
        }
    }

    private VCFRecord addVariantInformationToCall(ReferenceContext ref, PlinkRod plinkRod) {

        VariantContext vContext = plinkRod.getVariantContext();
        VCFRecord record = VariantContextAdaptors.toVCF(vContext);
        record.setGenotypeFormatString("GT");

        int numNoCalls = vContext.getNoCallCount();
        int numHomVarCalls = vContext.getHomVarCount();
        int numHetCalls = vContext.getHetCount();



        double noCallProp = ( (double) numNoCalls )/( (double) sampleNames.size());
        double homNonRProp = ( (double) numHomVarCalls )/( (double) sampleNames.size() - numNoCalls);

        record.setQual(DEFAULT_QUALITY);
        String hw = hardyWeinbergCalculation(ref,record);
        double hwScore = hw != null ? Double.valueOf(hw) : -0.0;
        // TODO -- record.addInfoFields(generateInfoField(record, numNoCalls,numHomVarCalls,numNonrefAlleles,ref, plinkRod, hwScore));
        if ( hwScore > maxHardy ) {
            record.setFilterString("Hardy-Weinberg");
        } else if ( noCallProp > maxNoCall ) {
            record.setFilterString("No-calls");
        } else if ( homNonRProp > maxHomNonref) {
            record.setFilterString("HomNonref-calls");
        }
        
        return record;

    }

    private String hardyWeinbergCalculation(ReferenceContext ref, VCFRecord rec) {
        if ( useSmartHardy ) {
            return smartHardy(ref, rec);
        } else {
            VCFVariationCall variant = new VCFVariationCall(ref.getBase(),ref.getLocus(),VCFVariationCall.VARIANT_TYPE.SNP);
            variant.setGenotypeCalls(rec.getGenotypes());
            return HWCalc.annotate(null, ref, null, variant);
        }
    }

    private Map<String,String> generateInfoField(VCFRecord rec, int nocall, int homnonref, int allnonref,
                                                 ReferenceContext ref, PlinkRod info, double hwScore) {
        double propNoCall = ( ( double ) nocall / (double) nSamples );
        double propHomNR = ( (double) homnonref / (double) nSamples );
        HashMap<String,String> infoMap = new HashMap<String,String>(1);
        putInfoStrings(infoMap,propNoCall,propHomNR,allnonref,hwScore,info.getName());

        return infoMap;
    }

    private void putInfoStrings(HashMap<String,String> infoMap, double pnc, double phnr, int nra, double hw, String nm) {

        infoMap.put("snpID",nm);
        infoMap.put("noCallPct",String.format("%.2f",100.0*pnc));
        infoMap.put("homNonrefPct",String.format("%.2f",100.0*phnr));
        infoMap.put("nonrefAlleles",String.format("%d",nra));
        infoMap.put("HW",String.format("%.2f",hw));

        //return String.format("snpID=%s;nocall=%f;homNonref=%4f;numNonrefAlleles=%d;HW=%s",nm,pnc,phnr,nra,hw);

    }

    private String smartHardy(ReferenceContext ref, VCFRecord rec) {
        HashMap<String,ArrayList<Genotype>> genotypesByPopulation = new HashMap<String,ArrayList<Genotype>>(INIT_NUMBER_OF_POPULATIONS);
        HashMap<String,String> hardyWeinbergByPopulation = new HashMap<String,String>(INIT_NUMBER_OF_POPULATIONS);

        for ( String population : samplesToPopulation.values() ) {
            genotypesByPopulation.put(population,new ArrayList<Genotype>());
        }

        for ( String name : sampleNames ) {
            String pop = samplesToPopulation.get(name);
            if ( rec.getGenotype(name) != null ) {
                // TODO -- genotypesByPopulation.get(pop).add(rec.getGenotype(name));
            }
        }

        for ( String population : samplesToPopulation.values() ) {
            VCFVariationCall v = new VCFVariationCall(ref.getBase(),ref.getLocus(),VCFVariationCall.VARIANT_TYPE.SNP);
            // TODO -- v.setGenotypeCalls(genotypesByPopulation.get(population));
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

    private HashMap<String,String> parsePopulationFile(File file) {
        HashMap<String,String> samplesToPopulation = new HashMap<String,String>();
        try {
            BufferedReader in = new BufferedReader( new FileReader( file ));
            String line = in.readLine();
            while ( line != null ) {
                String[] populationSamples = line.split("\t");
                String population = populationSamples[0];
                for ( int i = 1; i < populationSamples.length; i ++ ) {
                    samplesToPopulation.put(populationSamples[i],population);
                }
            }
        } catch ( IOException e) {
            throw new StingException("Error reading population file", e);
        }

        return samplesToPopulation;
    }
}

package org.broadinstitute.sting.playground.gatk.walkers.variantstovcf;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.walkers.annotator.HardyWeinberg;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.DiploidGenotype;
import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.Variation;
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Converts Sequenom files to a VCF annotated with QC metrics (HW-equilibrium, % failed probes)
 */
public class PlinkToVCF extends RefWalker<VCFRecord,Integer> {
    @Argument(fullName="sequenomePedFile", shortName="sPed", doc="The sequenome file from which to generate a VCF", required=true)
    public File seqFile = null;
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
    private HashMap<String, SequenomVariantInfo> sequenomResults = new HashMap<String,SequenomVariantInfo>();
    private ArrayList<String> sampleNames = new ArrayList<String>(nSamples);
    private VCFGenotypeWriterAdapter vcfWriter;
    private final HardyWeinberg HWCalc = new HardyWeinberg();
    private final boolean useSmartHardy = popFile != null;
    private HashMap<String,String> samplesToPopulation;

    public void initialize() {
        //System.out.println("Initializing... parse sequenom file");
        parseSequenomFile(seqFile);
        vcfWriter = new VCFGenotypeWriterAdapter(vcfFile);
        if ( useSmartHardy ) {
            samplesToPopulation = parsePopulationFile(popFile);
        }
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFHeaderLine("source", "Sequenom2VCF"));
        hInfo.add(new VCFHeaderLine("reference", getToolkit().getArguments().referenceFile.getName()));
        vcfWriter.writeHeader(new TreeSet<String>(sampleNames),hInfo);
        nSamples = sampleNames.size();
    }

    public Integer reduceInit() {
        int numberOfVariantsProcessed = 0;
        return numberOfVariantsProcessed;
    }

    public VCFRecord map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( sequenomResults.containsKey(context.getLocation().toString()) ) {
            SequenomVariantInfo varInfo = sequenomResults.remove(context.getLocation().toString());
            return addVariantInformationToCall(ref,varInfo);
        } else {
            return null;
        }
    }

    public Integer reduce(VCFRecord call, Integer numVariants) {
        if ( call == null ) {
            return numVariants;
        } else {
            printToVCF(call);
            return 1 + numVariants;
        }
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

    private VCFRecord addVariantInformationToCall(ReferenceContext ref, SequenomVariantInfo varInfo) {
        int numNoCalls = 0;
        int numHomNonrefCalls = 0;
        int numNonrefAlleles = 0;

        int sampleNumber = 0;
        VCFRecord record = new VCFRecord(ref.getBase(),ref.getLocus(),"GT");

        for ( String genTypeStr : varInfo.getGenotypes() ) {
            if ( genTypeStr.indexOf("0") == -1 ) {
                VCFGenotypeEncoding allele1 = new VCFGenotypeEncoding(genTypeStr.substring(0,1));
                VCFGenotypeEncoding allele2 = new VCFGenotypeEncoding(genTypeStr.substring(1));
                List<VCFGenotypeEncoding> alleles = new ArrayList<VCFGenotypeEncoding>(2);
                alleles.add(allele1);
                alleles.add(allele2);

                VCFGenotypeRecord genotype = new VCFGenotypeRecord(sampleNames.get(sampleNumber), alleles, VCFGenotypeRecord.PHASE.UNPHASED);
                genotype.setField("GQ",String.format("%d",DEFAULT_QUALITY));

                if ( genotype.isVariant(ref.getBase()) ) {
                    if ( genotype.isHom() ) {
                        numHomNonrefCalls++;
                        numNonrefAlleles+=2;
                        record.addAlternateBase(allele1);
                    } else {
                        numNonrefAlleles++;
                        record.addAlternateBase(allele1.getBases().equalsIgnoreCase(String.format("%c",ref.getBase())) ? allele2 : allele1);
                    }
                }

                record.addGenotypeRecord(genotype);

            } else {
                numNoCalls++;
            }
            sampleNumber++;
        }

        double noCallProp = ( (double) numNoCalls )/( (double) sampleNames.size());
        double homNonRProp = ( (double) numHomNonrefCalls )/( (double) sampleNames.size() - numNoCalls);

        record.setQual(DEFAULT_QUALITY);
        String hw = hardyWeinbergCalculation(ref,record);
        double hwScore = hw != null ? Double.valueOf(hw) : -0.0;
        record.addInfoFields(generateInfoField(record, numNoCalls,numHomNonrefCalls,numNonrefAlleles,ref, varInfo, hwScore));
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
                                                 ReferenceContext ref, SequenomVariantInfo info, double hwScore) {
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
                genotypesByPopulation.get(pop).add(rec.getGenotype(name));
            }
        }

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

    private void parseSequenomFile(File sequenomFile) {
        BufferedReader reader;
        try {
            reader = new BufferedReader(new FileReader(sequenomFile));
        } catch ( IOException e ) {
            throw new StingException("Sequenom file could not be opened",e);
        }
        String header;
        try {
            header = reader.readLine();
            //System.out.println("Read header line, it was "+header);
        } catch (IOException e) {
            throw new StingException(e.getMessage(),e);
        }
        HashMap<Integer,SequenomVariantInfo> sequenomVariants = parseSequenomHeader(header);
        try {
            String line;
            do {
                line = reader.readLine();
                //ystem.out.println("Read line, it was"+line);
                if ( line != null ) {
                    //System.out.println("Parsing line...");
                    parseSequenomLine(sequenomVariants,line);
                }
            } while ( line != null);

            reader.close();
            convertToLocusMap(sequenomVariants);

        } catch ( IOException e) {
            throw new StingException("Error reading sequenom file", e);
        }
    }

    private HashMap<Integer,SequenomVariantInfo> parseSequenomHeader(String header) {
        // the header will give us all the probe names and contigs
        //System.out.println("Parsing header");
        HashMap<Integer,SequenomVariantInfo> variantsByFileOffset = new HashMap<Integer,SequenomVariantInfo>();
        String[] fields = header.split("\t");
        int fieldOffset = 0;
        for ( String entry : fields ) {
            if ( ! HEADER_FIELDS.contains(entry) ) {
                //System.out.println(entry);
                // actually a SNP
                String snpName = entry.split("\\|")[1];
                //System.out.println("Entry: "+entry+" Name: "+snpName);
                variantsByFileOffset.put(fieldOffset,new SequenomVariantInfo(snpName,nSamples));
                //System.out.println("Adding entry for offset "+fieldOffset);
            }
            fieldOffset++;
        }

        return variantsByFileOffset;
    }

    private void parseSequenomLine(HashMap<Integer,SequenomVariantInfo> variants, String line) {
        if ( line == null ) {
            // early return
            //System.out.println("Line was null in file");
            return;
        }
        String[] entries = line.split("\t");
        int offset = 0;
        if ( entries[1].equalsIgnoreCase("empty")) {
            return;
        }
        sampleNames.add(entries[1]);
        for ( String entry : entries ) {
            if ( variants.containsKey(offset) ) { // actual SNP
                variants.get(offset).addGenotype(entry);
                //System.out.println("Added: "+entry+"To "+offset);
            }
            offset++;
        }
    }

    private void convertToLocusMap(HashMap<Integer,SequenomVariantInfo> variants) {
        for ( SequenomVariantInfo variant : variants.values() ) {
            //System.out.println("Variant name: "+variant.getName());
            String loc = genomeLocFromVariantName(variant.getName());
            //System.out.println(variant.getName());
            if ( loc == null ) {
                throw new StingException("Genome locus was null");
            }
            sequenomResults.put(loc,variant);
        }

        variants.clear();
        variants = null;
    }

    private String genomeLocFromVariantName(String name) {
        String[] nameInfo = name.split("_");
        //System.out.println(name);
        //System.out.println(nameInfo[0]);
        String chr = nameInfo[0].substring(1); // format: c3 becomes 3
        String pos = nameInfo[1].substring(1); // format: p9961104 becomes 9961104
        if ( useb36contigs ) {
            return chr+":"+pos;
        } else {
            return "chr"+chr+":"+pos;
        }

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

class SequenomVariantInfo {
    private ArrayList<String> genotypes;
    private String variantName;

    public SequenomVariantInfo(String name, int nPeople) {
        genotypes = new ArrayList<String>(nPeople);
        variantName = name;
    }

    public void addGenotype(String genotype) {
        String[] alleles = genotype.split(" ");
        genotypes.add(alleles[0]+alleles[1]);
    }

    public String getName() {
        return variantName;
    }

    public ArrayList<String> getGenotypes() {
        return genotypes;
    }
}

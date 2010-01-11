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
import org.broadinstitute.sting.utils.genotype.vcf.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Converts Sequenom files to a VCF annotated with QC metrics (HW-equilibrium, % failed probes)
 */
public class SequenomToVCF extends RefWalker<VCFVariationCall,Integer> {
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
        hInfo.add(new VCFHeaderLine("annotatorReference", getToolkit().getArguments().referenceFile.getName()));
        vcfWriter.writeHeader(new TreeSet<String>(sampleNames),hInfo);
        nSamples = sampleNames.size();
    }

    public Integer reduceInit() {
        int numberOfVariantsProcessed = 0;
        return numberOfVariantsProcessed;
    }

    public VCFVariationCall map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        if ( sequenomResults.containsKey(context.getLocation().toString()) ) {
            SequenomVariantInfo varInfo = sequenomResults.remove(context.getLocation().toString());
            return addVariantInformationToCall(ref,varInfo);
        } else {
            return null;
        }
    }

    public Integer reduce(VCFVariationCall call, Integer numVariants) {
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

    private void printToVCF(VCFVariationCall call) {
        vcfWriter.addMultiSampleCall(call.getGenotypes(),call);
    }

    private VCFVariationCall addVariantInformationToCall(ReferenceContext ref, SequenomVariantInfo varInfo) {
        int numNoCalls = 0;
        int numHomNonrefCalls = 0;
        int numNonrefAlleles = 0;

        int sampleNumber = 0;
        //System.out.println("Genotypes from varinfo:");
        //for ( String g : varInfo.getGenotypes()) {
            //System.out.println(g);
        //}
        ArrayList<VCFGenotypeCall> vcfGenotypeCalls = new ArrayList<VCFGenotypeCall>(nSamples);
        ArrayList<Genotype> genotypeCalls = new ArrayList<Genotype>(nSamples);
        VCFGenotypeCall vcfCall = new VCFGenotypeCall(ref.getBase(),ref.getLocus());
        boolean isSNP = false;

        for ( String genTypeStr : varInfo.getGenotypes() ) {
            if ( genTypeStr.indexOf("0") == -1 ) {
                vcfCall.setGenotype(DiploidGenotype.createDiploidGenotype(genTypeStr.charAt(0),genTypeStr.charAt(2)));
                vcfCall.setNegLog10PError((double) DEFAULT_QUALITY/10);
                vcfCall.setSampleName(sampleNames.get(sampleNumber));
                genotypeCalls.add( vcfCall.cloneCall() );
                vcfGenotypeCalls.add( vcfCall.cloneCall() );
                if ( vcfCall.isVariant(ref.getBase()) ) {
                    isSNP = true;
                    if ( vcfCall.isHom() ) {
                        numHomNonrefCalls++;
                        numNonrefAlleles+=2;
                    } else {
                        numNonrefAlleles++;
                    }
                }
            } else {
                numNoCalls++;
            }
            sampleNumber++;
        }

        VCFVariationCall variantCall = new VCFVariationCall(ref.getBase(),ref.getLocus(), isSNP ? VCFVariationCall.VARIANT_TYPE.SNP : VCFVariationCall.VARIANT_TYPE.REFERENCE);
        variantCall.setGenotypeCalls(genotypeCalls);
        variantCall.setConfidence((double) DEFAULT_QUALITY);
        variantCall.setFields(generateInfoField(numNoCalls,numHomNonrefCalls,numNonrefAlleles,variantCall,ref, varInfo, vcfGenotypeCalls));

        return variantCall;
    }

    private Map<String,String> generateInfoField(int nocall, int homnonref, int allnonref, VCFVariationCall call,
                                                 ReferenceContext ref, SequenomVariantInfo info, List<VCFGenotypeCall> vcfCalls) {
        double propNoCall = ( ( double ) nocall / (double) nSamples );
        double propHomNR = ( (double) homnonref / (double) nSamples );
        String hardy;
        if ( useSmartHardy ) {
            hardy = smartHardy(ref, call, info, vcfCalls);
        } else {
            hardy = HWCalc.annotate(null,ref, null, call);
        }
        HashMap<String,String> infoMap = new HashMap<String,String>(1);
        putInfoStrings(infoMap,propNoCall,propHomNR,allnonref,hardy,info.getName());

        return infoMap;
    }

    private void putInfoStrings(HashMap<String,String> infoMap, double pnc, double phnr, int nra, String hw, String nm) {

        infoMap.put("snpID",nm);
        infoMap.put("noCallPct",String.format("%.2f",100.0*pnc));
        infoMap.put("homNonrefPct",String.format("%.2f",100.0*phnr));
        infoMap.put("nonrefAlleles",String.format("%d",nra));
        infoMap.put("HW",hw);

        //return String.format("snpID=%s;nocall=%f;homNonref=%4f;numNonrefAlleles=%d;HW=%s",nm,pnc,phnr,nra,hw);

    }

    private String smartHardy(ReferenceContext ref, VCFVariationCall call, SequenomVariantInfo info, List<VCFGenotypeCall> vcfCalls) {
        HashMap<String,ArrayList<Genotype>> genotypesByPopulation = new HashMap<String,ArrayList<Genotype>>(INIT_NUMBER_OF_POPULATIONS);
        HashMap<String,String> hardyWeinbergByPopulation = new HashMap<String,String>(INIT_NUMBER_OF_POPULATIONS);

        for ( String population : samplesToPopulation.values() ) {
            genotypesByPopulation.put(population,new ArrayList<Genotype>());
        }

        for ( VCFGenotypeCall vgc : vcfCalls ) {
            genotypesByPopulation.get(vgc.getSampleName()).add(vgc);
        }

        for ( String population : samplesToPopulation.values() ) {
            VCFVariationCall v = new VCFVariationCall(ref.getBase(),ref.getLocus(),VCFVariationCall.VARIANT_TYPE.SNP);
            v.setGenotypeCalls(genotypesByPopulation.get(population));
            hardyWeinbergByPopulation.put(population,HWCalc.annotate(null,ref,null,v));
        }

        return smartHardyString(hardyWeinbergByPopulation,info);
    }

    private String smartHardyString(HashMap<String,String> hwByPop, SequenomVariantInfo varInfo) {
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
            if ( fieldOffset > 5 ) {
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
            if ( offset > 5 ) { // actual SNP
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
        genotypes.add(genotype);
    }

    public String getName() {
        return variantName;
    }

    public ArrayList<String> getGenotypes() {
        return genotypes;
    }
}

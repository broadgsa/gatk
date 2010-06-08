/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.gatk.walkers;

import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.StratifiedAlignmentContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.annotator.VariantAnnotatorEngine;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.genotype.vcf.VCFReader;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;

import java.io.*;
import java.util.*;
import static java.lang.Math.log10;


/**
 * Takes files produced by Beagle imputation engine and creates a vcf with modified annotations.
 */
@Requires(value={},referenceMetaData=@RMD(name=BeagleOutputToVCFWalker.INPUT_ROD_NAME,type= VCFRecord.class))
public class BeagleOutputToVCFWalker  extends RodWalker<Integer, Integer> {

    private VCFWriter vcfWriter;

    @Argument(fullName="input_prefix", shortName="input", doc="The prefix added to input Beagle files gprobs, r2, ...", required=true)
    private String INPUT_PREFIX = "beagle";

    @Argument(fullName="output_file", shortName="output", doc="VCF file to which output should be written", required=true)
    private String OUTPUT_FILE = null;

    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", required=false)
    protected String[] annotationsToUse = {"AlleleBalance"};

    @Argument(fullName="group", shortName="G", doc="One or more classes/groups of annotations to apply to variant calls", required=false)
    protected String[] annotationClassesToUse = {};

    @Argument(fullName="useAllAnnotations", shortName="all", doc="Use all possible annotations (not for the faint of heart)", required=false)
    protected Boolean USE_ALL_ANNOTATIONS = false;


    public static final String INPUT_ROD_NAME = "inputvcf";

    protected static BeagleFileReader gprobsReader = null;
    protected static BeagleFileReader phasedReader = null;
    protected static BeagleFileReader likeReader = null;
    protected static BeagleFileReader r2Reader = null;

    private VariantAnnotatorEngine engine;

    protected static String line = null;

    protected HashMap<String,BeagleSampleRecord> beagleSampleRecords;

    final TreeSet<String> samples = new TreeSet<String>();


    private final double MIN_PROB_ERROR = 0.000001;
    private final double MAX_GENOTYPE_QUALITY = 6.0;

    public void initialize() {
        if ( USE_ALL_ANNOTATIONS )
            engine = new VariantAnnotatorEngine(getToolkit());
        else
            engine = new VariantAnnotatorEngine(getToolkit(), annotationClassesToUse, annotationsToUse);

        // setup the header fields

        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFInfoHeaderLine("R2", 1, VCFInfoHeaderLine.INFO_TYPE.Float, "r2 Value reported by Beable on each site"));
        hInfo.add(new VCFHeaderLine("source", "BeagleImputation"));
        hInfo.addAll(engine.getVCFAnnotationDescriptions());

        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();

        // Open output file specified by output VCF ROD
        vcfWriter = new VCFWriter(new File(OUTPUT_FILE));

        for( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();
            if( rod.getRecordType().equals(VCFRecord.class) && rod.getName().equalsIgnoreCase(INPUT_ROD_NAME)) {
                final VCFReader reader = new VCFReader(rod.getFile());
                final Set<String> vcfSamples = reader.getHeader().getGenotypeSamples();
                samples.addAll(vcfSamples);
                reader.close();
            }

        }
        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);

        // Open Beagle files for processing
        gprobsReader = new BeagleFileReader(INPUT_PREFIX.concat(".gprobs"));
        likeReader = new BeagleFileReader(INPUT_PREFIX.concat(".like"));
        phasedReader = new BeagleFileReader(INPUT_PREFIX.concat(".phased"));
        r2Reader = new BeagleFileReader(INPUT_PREFIX.concat(".r2"));

        if (!CheckAllSamplesAreInList(likeReader.getHeaderString().split(" "), samples, 3, 3))
            throw new StingException("Inconsistent list of samples in File: " + likeReader.getFileName());

        if (!CheckAllSamplesAreInList(gprobsReader.getHeaderString().split(" "), samples, 3, 3))
            throw new StingException("Inconsistent list of samples in File: " + gprobsReader.getFileName());

        if (!CheckAllSamplesAreInList(phasedReader.getHeaderString().split(" "), samples, 2, 2))
            throw new StingException("Inconsistent list of samples in File: " + phasedReader.getFileName());


        beagleSampleRecords = new HashMap<String,BeagleSampleRecord>();

        // OK, now that we have all data in string array format inside readers, process each array individually,
        for(String[] gprobLine: gprobsReader.getDataSet())
        {
            int j;
            BeagleSampleRecord bglRecord = new BeagleSampleRecord();
            LikelihoodTrios ltrio;

            String markerKey = gprobLine[0];
            String alleleA = gprobLine[1];
            String alleleB = gprobLine[2];

            HashMap<String,LikelihoodTrios> genotypeProbabilities = new HashMap<String,LikelihoodTrios>();

            j = 3;
            for (String sample : samples) {
                ltrio = new LikelihoodTrios(new Double(gprobLine[j]), new Double(gprobLine[j+1]),
                        new Double(gprobLine[j+2]));
                j = j+3;
                genotypeProbabilities.put(sample,ltrio);
            }

            bglRecord.GenotypeProbabilities = genotypeProbabilities;
            bglRecord.AlleleA = Allele.create(alleleA);
            bglRecord.AlleleB = Allele.create(alleleB);

            beagleSampleRecords.put(markerKey,bglRecord);
        }

        // now fill in input likelihood info in the same way
        for (String[] likeLine: likeReader.getDataSet())
        {
            int j;

            LikelihoodTrios ltrio;

            String markerKey = likeLine[0];

            HashMap<String,LikelihoodTrios> genotypeLikelihoods = new HashMap<String,LikelihoodTrios>();

            j = 3;
            for (String sample : samples) {
                ltrio = new LikelihoodTrios(new Double(likeLine[j]), new Double(likeLine[j+1]),
                        new Double(likeLine[j+2]));
                j = j+3;
                genotypeLikelihoods.put(sample,ltrio);
            }

            BeagleSampleRecord bglRecord = beagleSampleRecords.get(markerKey);
            bglRecord.InputLikelihoods = genotypeLikelihoods;
            beagleSampleRecords.put(markerKey,bglRecord);

        }

        for (String[] phasedLine: phasedReader.getDataSet())
        {
            int j;

            HaplotypePair pair;

            String markerKey = phasedLine[1];

            HashMap<String,HaplotypePair> haplotypePairs = new HashMap<String,HaplotypePair>();

            j = 2;
            for (String sample : samples) {
                pair = new HaplotypePair(Allele.create(phasedLine[j]), Allele.create(phasedLine[j+1]));
                j = j+2;
                haplotypePairs.put(sample,pair);
            }

            BeagleSampleRecord bglRecord = beagleSampleRecords.get(markerKey);
            bglRecord.PhasedHaplotypes = haplotypePairs;
            beagleSampleRecords.put(markerKey,bglRecord);

        }

        for (String[] r2Line: r2Reader.getDataSet())
        {

            String[] vals = r2Line[0].split("\t");
            String markerKey = vals[0];

            BeagleSampleRecord bglRecord = beagleSampleRecords.get(markerKey);
            bglRecord.setR2Value(Double.valueOf(vals[1]));
            beagleSampleRecords.put(markerKey,bglRecord);

        }
    }

    private boolean CheckAllSamplesAreInList(String[] stringArray, TreeSet<String> samples, int valsPerSample, int initMarkers)
    {
        boolean allThere = true;

        if (stringArray.length != valsPerSample*samples.size()+initMarkers) {
            allThere = false;
        }

        for (int k=initMarkers; k < stringArray.length; k++)
            if (!samples.contains(stringArray[k]))
                allThere = false;

        return allThere;

    }


    private class BeagleFileReader {
        private String headerString;
        private BufferedReader reader;
        private String fileName;
        private HashSet<String[]> dataSet;

        public String getHeaderString() {
            return headerString;
        }
        public BufferedReader getReader() {
            return reader;
        }
        public String getFileName() {
            return fileName;
        }
        public HashSet<String[]> getDataSet() {
            return dataSet;
        }

        public BeagleFileReader(String fileName) {
            String curLine;

            try{
                reader = new BufferedReader(new FileReader(fileName));
            } catch ( FileNotFoundException e) {
                throw new StingException("Could not find required input file: " + fileName);
            }

            this.fileName = fileName;

            try {
                headerString = reader.readLine();

                dataSet = new HashSet<String[]>();
                while (( curLine = reader.readLine()) != null){
                    // Split current line along white spaces, ignore multiple white spaces
                    String lineTokens[] = curLine.split(" +");

                    dataSet.add(lineTokens);
                }
                reader.close();
            }
            catch (IOException e) {
                throw new StingException("Failed while reading data from input file: " + fileName);
            }

        }

    }

    private class LikelihoodTrios {
        public double HomRefLikelihood;
        public double HetLikelihood;
        public double HomVarLikelihood;

        public LikelihoodTrios(double hr, double het, double hv) {
            this.HomRefLikelihood = hr;
            this.HetLikelihood = het;
            this.HomVarLikelihood = hv;
        }

    }

    private class HaplotypePair {
        public Allele AlleleA;
        public Allele AlleleB;

        public HaplotypePair(Allele a, Allele b) {
            this.AlleleA = a;
            this.AlleleB = b;
        }
    }

    private class BeagleSampleRecord {
        // class containing, for each marker, all information necessary for Beagle input/output
        public Allele AlleleA;
        public Allele AlleleB;
        private double r2Value;
        public HashMap<String,LikelihoodTrios> GenotypeProbabilities;
        public HashMap<String,LikelihoodTrios> InputLikelihoods;
        public HashMap<String,HaplotypePair> PhasedHaplotypes;

        public Double getR2Value() {
            return r2Value;
        }

        public void setR2Value(Double x) {
            this.r2Value = x;
        }
    }


    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if ( tracker == null )
            return 0;

        EnumSet<VariantContext.Type> vc = EnumSet.of(VariantContext.Type.SNP);
        GenomeLoc loc = context.getLocation();
        VariantContext vc_input;


        try {
            vc_input = tracker.getVariantContext(ref,"inputvcf", vc, loc, true);
        } catch (java.util.NoSuchElementException e) {
            return 0;
        }



        if (beagleSampleRecords.containsKey(loc.toString())) {
            // Get record associated with this marker from output Beagle data
            BeagleSampleRecord bglRecord = beagleSampleRecords.get(loc.toString());

            // get reference base for current position
            byte refByte = ref.getBase();

            // make new Genotypes based on Beagle results
            Map<String, Genotype> genotypes;

            genotypes = new HashMap<String, Genotype>(vc_input.getGenotypes().size());

            // for each genotype, create a new object with Beagle information on it
            for ( Map.Entry<String, Genotype> genotype : vc_input.getGenotypes().entrySet() ) {

                Genotype g = genotype.getValue();
                Set<String> filters = new LinkedHashSet<String>(g.getFilters());

                // modify g here with Beagle info
                boolean genotypeIsPhased = true;
                String sample = g.getSampleName();

                LikelihoodTrios gtprobs = bglRecord.GenotypeProbabilities.get(sample);
//                LikelihoodTrios inplike = bglRecord.InputLikelihoods.get(sample);
                HaplotypePair hp = bglRecord.PhasedHaplotypes.get(sample);

                // We have phased genotype in hp. Need to set the isRef field in the allele.
                List<Allele> alleles = new ArrayList<Allele>();

                byte r[] = hp.AlleleA.getBases();
                byte rA = r[0];

                Boolean isRefA = (refByte  == rA);

                Allele refAllele = Allele.create(r, isRefA );
                alleles.add(refAllele);

                r = hp.AlleleB.getBases();
                byte rB = r[0];

                Boolean isRefB = (refByte  == rB);
                Allele altAllele = Allele.create(r,isRefB);
                alleles.add(altAllele);



                // Compute new GQ field = -10*log10Pr(Genotype call is wrong)
                // Beagle gives probability that genotype is AA, AB and BB.
                // Which, by definition, are prob of hom ref, het and hom var.
                Double probWrongGenotype, genotypeQuality;

                if (isRefA && isRefB) // HomRef call
                    probWrongGenotype = gtprobs.HetLikelihood + gtprobs.HomVarLikelihood;
                else if ((isRefB && !isRefA) || (isRefA && !isRefB))
                    probWrongGenotype = gtprobs.HomRefLikelihood + gtprobs.HomVarLikelihood;
                else // HomVar call
                    probWrongGenotype = gtprobs.HetLikelihood + gtprobs.HomRefLikelihood;


                if (probWrongGenotype < MIN_PROB_ERROR)
                    genotypeQuality = MAX_GENOTYPE_QUALITY;
                else
                    genotypeQuality = -log10(probWrongGenotype);

                Genotype imputedGenotype = new Genotype(genotype.getKey(), alleles, genotypeQuality, filters, g.getAttributes(), genotypeIsPhased);

                genotypes.put(genotype.getKey(), imputedGenotype);
            }


            VariantContext filteredVC = new VariantContext("outputvcf", vc_input.getLocation(), vc_input.getAlleles(), genotypes, vc_input.getNegLog10PError(), vc_input.getFilters(), vc_input.getAttributes());


            Set<Allele> altAlleles = filteredVC.getAlternateAlleles();
            StringBuffer altAlleleCountString = new StringBuffer();
            for ( Allele allele : altAlleles ) {
                if ( altAlleleCountString.length() > 0 )
                    altAlleleCountString.append(",");
                altAlleleCountString.append(filteredVC.getChromosomeCount(allele));
            }


            // if the reference base is not ambiguous, we can annotate
            Collection<VariantContext> annotatedVCs = Arrays.asList(filteredVC);
            Map<String, StratifiedAlignmentContext> stratifiedContexts;
            if ( BaseUtils.simpleBaseToBaseIndex(ref.getBase()) != -1 ) {
                if ( ! context.hasExtendedEventPileup() ) {
                    stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getBasePileup());
                } else {
                    stratifiedContexts = StratifiedAlignmentContext.splitContextBySample(context.getExtendedEventPileup());
                }
                if ( stratifiedContexts != null ) {
                    annotatedVCs = engine.annotateContext(tracker, ref, stratifiedContexts, filteredVC);
                }
            }


            for(VariantContext annotatedVC : annotatedVCs ) {
                VCFRecord vcf = VariantContextAdaptors.toVCF(filteredVC, ref.getBase());

                if ( annotatedVC.getChromosomeCount() > 0 ) {
                    vcf.addInfoField(VCFRecord.ALLELE_NUMBER_KEY, String.format("%d", annotatedVC.getChromosomeCount()));
                    if ( altAlleleCountString.length() > 0 )  {
                        vcf.addInfoField(VCFRecord.ALLELE_COUNT_KEY, altAlleleCountString.toString());
                        vcf.addInfoField(VCFRecord.ALLELE_FREQUENCY_KEY, String.format("%4.2f",
                                Double.valueOf(altAlleleCountString.toString())/(annotatedVC.getChromosomeCount())));
                    }
                }

                vcf.addInfoField("R2", (bglRecord.getR2Value()).toString() );
                vcfWriter.addRecord(vcf);
            }
        }

        return 1;

    }

    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    /**
     * Increment the number of loci processed.
     *
      * @param value result of the map.
      * @param sum   accumulator for the reduce.
      * @return the new number of loci processed.
      */
     public Integer reduce(Integer value, Integer sum) {
         return sum + value;
     }

     /**
      * Tell the user the number of loci processed and close out the new variants file.
      *
      * @param result  the number of loci seen.
      */
     public void onTraversalDone(Integer result) {
         out.printf("Processed %d loci.\n", result);

         vcfWriter.close();
     }
 }

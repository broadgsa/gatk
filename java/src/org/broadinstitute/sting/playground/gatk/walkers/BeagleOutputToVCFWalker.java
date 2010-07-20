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

import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.features.beagle.BeagleFeature;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFUtils;
import org.broadinstitute.sting.utils.genotype.vcf.VCFWriter;
import org.broad.tribble.vcf.*;

import java.io.*;
import java.util.*;
import static java.lang.Math.log10;


/**
 * Takes files produced by Beagle imputation engine and creates a vcf with modified annotations.
 */
@Requires(value={},referenceMetaData=@RMD(name=BeagleOutputToVCFWalker.INPUT_ROD_NAME, type=ReferenceOrderedDatum.class))

public class BeagleOutputToVCFWalker  extends RodWalker<Integer, Integer> {

    public static final String INPUT_ROD_NAME = "variant";
    public static final String COMP_ROD_NAME = "comp";
    public static final String R2_ROD_NAME = "beagleR2";
    public static final String PROBS_ROD_NAME = "beagleProbs";
    public static final String PHASED_ROD_NAME = "beaglePhased";

    @Argument(fullName="output_file", shortName="output", doc="VCF file to which output should be written", required=true)
    private String OUTPUT_FILE = null;

    @Argument(fullName="nocall_threshold", shortName="ncthr", doc="Threshold of confidence at which a genotype won't be called", required=false)
    private double noCallThreshold = 0.0;

    protected static String line = null;
    private VCFWriter vcfWriter;



    private final double MIN_PROB_ERROR = 0.000001;
    private final double MAX_GENOTYPE_QUALITY = 6.0;

    public void initialize() {

        // setup the header fields

        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFFormatHeaderLine("OG",1, VCFHeaderLineType.String, "Original Genotype input to Beagle"));
        hInfo.add(new VCFInfoHeaderLine("R2", 1, VCFHeaderLineType.Float, "r2 Value reported by Beagle on each site"));
        hInfo.add(new VCFInfoHeaderLine("NumGenotypesChanged", 1, VCFHeaderLineType.Integer, "r2 Value reported by Beagle on each site"));

        hInfo.add(new VCFHeaderLine("source", "BeagleImputation"));

        // Open output file specified by output VCF ROD
        vcfWriter = new VCFWriter(new File(OUTPUT_FILE));
        final List<ReferenceOrderedDataSource> dataSources = this.getToolkit().getRodDataSources();

        for( final ReferenceOrderedDataSource source : dataSources ) {
            final RMDTrack rod = source.getReferenceOrderedData();

            if (rod.getName().equals(COMP_ROD_NAME)) {
                hInfo.add(new VCFInfoHeaderLine("ACH", 1, VCFHeaderLineType.Integer, "Allele Count from Hapmap at this site"));
                hInfo.add(new VCFInfoHeaderLine("ANH", 1, VCFHeaderLineType.Integer, "Allele Frequency from Hapmap at this site"));
                hInfo.add(new VCFInfoHeaderLine("AFH", 1, VCFHeaderLineType.Float, "Allele Number from Hapmap at this site"));
                break;
            }

        }

        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(INPUT_ROD_NAME));

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if ( tracker == null )
            return 0;

        GenomeLoc loc = context.getLocation();
        VariantContext vc_input = tracker.getVariantContext(ref,INPUT_ROD_NAME, null, loc, false);

        VariantContext vc_comp = tracker.getVariantContext(ref,COMP_ROD_NAME, null, loc, false);

        if ( vc_input == null  )
            return 0;

        List<Object> r2rods = tracker.getReferenceMetaData(R2_ROD_NAME);

        // ignore places where we don't have a variant
        if ( r2rods.size() == 0 )
            return 0;

        BeagleFeature beagleR2Feature = (BeagleFeature)r2rods.get(0);

        List<Object> gProbsrods = tracker.getReferenceMetaData(PROBS_ROD_NAME);

        // ignore places where we don't have a variant
        if ( gProbsrods.size() == 0 )
            return 0;

        BeagleFeature beagleProbsFeature = (BeagleFeature)gProbsrods.get(0);

        List<Object> gPhasedrods = tracker.getReferenceMetaData(PHASED_ROD_NAME);

        // ignore places where we don't have a variant
        if ( gPhasedrods.size() == 0 )
            return 0;

        BeagleFeature beaglePhasedFeature = (BeagleFeature)gPhasedrods.get(0);

        // get reference base for current position
        byte refByte = ref.getBase();

        // make new Genotypes based on Beagle results
        Map<String, Genotype> genotypes = new HashMap<String, Genotype>(vc_input.getGenotypes().size());


        // for each genotype, create a new object with Beagle information on it

        int numGenotypesChangedByBeagle = 0;
        Integer alleleCountH = 0, chrCountH = 0;
        Double alleleFrequencyH = 0.0;

        Map<String,Genotype> hapmapGenotypes = null;

        if (vc_comp != null) {
            hapmapGenotypes = vc_comp.getGenotypes();
        }

        for ( Map.Entry<String, Genotype> originalGenotypes : vc_input.getGenotypes().entrySet() ) {

            Genotype g = originalGenotypes.getValue();
            Set<String> filters = new LinkedHashSet<String>(g.getFilters());

            boolean genotypeIsPhased = true;
            String sample = g.getSampleName();

            // If we have  a Hapmap (comp) ROD, compute Hapmap AC, AN and AF
            // use sample as key into genotypes structure
            if (vc_comp != null) {

                if (vc_input.getGenotypes().containsKey(sample) && hapmapGenotypes.containsKey(sample))  {

                    Genotype hapmapGenotype = hapmapGenotypes.get(sample);
                    if (hapmapGenotype.isCalled()){
                        chrCountH += 2;
                        if (hapmapGenotype.isHet()) {
                            alleleCountH += 1;
                        }    else if (hapmapGenotype.isHomVar()) {
                            alleleCountH += 2;
                        }
                    }
                }
            }

            ArrayList<String> beagleProbabilities = beagleProbsFeature.getProbLikelihoods().get(sample);
            ArrayList<String> beagleGenotypePairs = beaglePhasedFeature.getGenotypes().get(sample);

            Allele originalAlleleA = g.getAllele(0);
            Allele originalAlleleB = g.getAllele(1);

            // We have phased genotype in hp. Need to set the isRef field in the allele.
            List<Allele> alleles = new ArrayList<Allele>();

            String alleleA = beagleGenotypePairs.get(0);
            String alleleB = beagleGenotypePairs.get(1);

            byte[] r = alleleA.getBytes();
            byte rA = r[0];

            Boolean isRefA = (refByte  == rA);

            Allele refAllele = Allele.create(r, isRefA );
            alleles.add(refAllele);

            r = alleleB.getBytes();
            byte rB = r[0];

            Boolean isRefB = (refByte  == rB);
            Allele altAllele = Allele.create(r,isRefB);
            alleles.add(altAllele);

            // Compute new GQ field = -10*log10Pr(Genotype call is wrong)
            // Beagle gives probability that genotype is AA, AB and BB.
            // Which, by definition, are prob of hom ref, het and hom var.
            Double probWrongGenotype, genotypeQuality;
            Double homRefProbability = Double.valueOf(beagleProbabilities.get(0));
            Double hetProbability = Double.valueOf(beagleProbabilities.get(1));
            Double homVarProbability = Double.valueOf(beagleProbabilities.get(2));

            if (isRefA && isRefB) // HomRef call
                probWrongGenotype = hetProbability + homVarProbability;
            else if ((isRefB && !isRefA) || (isRefA && !isRefB))
                probWrongGenotype = homRefProbability + homVarProbability;
            else // HomVar call
                probWrongGenotype = hetProbability + homRefProbability;

            if (1-probWrongGenotype < noCallThreshold) {
                // quality is bad: don't call genotype
                alleles.clear();
                refAllele = originalAlleleA;
                altAllele = originalAlleleB;
                alleles.add(refAllele);
                alleles.add(altAllele);
                genotypeIsPhased = false;
            }

            if (probWrongGenotype < MIN_PROB_ERROR)
                genotypeQuality = MAX_GENOTYPE_QUALITY;
            else
                genotypeQuality = -log10(probWrongGenotype);

            HashMap<String,Object> originalAttributes = new HashMap<String,Object>(g.getAttributes());

            // get original encoding and add to keynotype attributes
            String a1, a2, og;
            if (originalAlleleA.isNoCall())
                a1 = ".";
            else if (originalAlleleA.isReference())
                a1 = "0";
            else
                a1 = "1";

            if (originalAlleleB.isNoCall())
                a2 = ".";
            else if (originalAlleleB.isReference())
                a2 = "0";
            else
                a2 = "1";

            og = a1+"/"+a2;

            // See if Beagle switched genotypes
            if (!((refAllele.equals(originalAlleleA) && altAllele.equals(originalAlleleB) ||
                    (refAllele.equals(originalAlleleB) && altAllele.equals(originalAlleleA))))){
                originalAttributes.put("OG",og);
                numGenotypesChangedByBeagle++;
            }
            else {
                originalAttributes.put("OG",".");
            }
            Genotype imputedGenotype = new Genotype(originalGenotypes.getKey(), alleles, genotypeQuality, filters,originalAttributes , genotypeIsPhased);


            genotypes.put(originalGenotypes.getKey(), imputedGenotype);

        }

        VariantContext filteredVC = new VariantContext("outputvcf", vc_input.getLocation(), vc_input.getAlleles(), genotypes, vc_input.getNegLog10PError(), vc_input.filtersWereApplied() ? vc_input.getFilters() : null, vc_input.getAttributes());

        Set<Allele> altAlleles = filteredVC.getAlternateAlleles();
        StringBuffer altAlleleCountString = new StringBuffer();
        for ( Allele allele : altAlleles ) {
            if ( altAlleleCountString.length() > 0 )
                altAlleleCountString.append(",");
            altAlleleCountString.append(filteredVC.getChromosomeCount(allele));
        }

        HashMap<String, Object> attributes = new HashMap<String, Object>(filteredVC.getAttributes());
        if ( filteredVC.getChromosomeCount() > 0 ) {
            attributes.put(VCFConstants.ALLELE_NUMBER_KEY, String.format("%d", filteredVC.getChromosomeCount()));
            if ( altAlleleCountString.length() > 0 )  {
                attributes.put(VCFConstants.ALLELE_COUNT_KEY, altAlleleCountString.toString());
                attributes.put(VCFConstants.ALLELE_FREQUENCY_KEY, String.format("%4.2f",
                        Double.valueOf(altAlleleCountString.toString())/(filteredVC.getChromosomeCount())));
            }
        }

        // Get Hapmap AC and AF
        if (vc_comp != null) {
            attributes.put("ACH", alleleCountH.toString() );
            attributes.put("ANH", chrCountH.toString() );
            attributes.put("AFH", String.format("%4.2f", (double)alleleCountH/chrCountH) );

        }

        attributes.put("NumGenotypesChanged", numGenotypesChangedByBeagle );
        attributes.put("R2", beagleR2Feature.getR2value().toString() );



        vcfWriter.add(VariantContextUtils.modifyAttributes(filteredVC, attributes), new byte[]{ref.getBase()});


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

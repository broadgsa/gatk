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
import org.broadinstitute.sting.gatk.contexts.variantcontext.*;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.features.beagle.BeagleFeature;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.VariantContextAdaptors;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrack;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
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

    @Argument(fullName="output_file", shortName="output", doc="VCF file to which output should be written", required=true)
    private String OUTPUT_FILE = null;


    public static final String INPUT_ROD_NAME = "inputvcf";

    protected static String line = null;

    // protected HashMap<String,BeagleSampleRecord> beagleSampleRecords;

    final TreeSet<String> samples = new TreeSet<String>();


    private final double MIN_PROB_ERROR = 0.000001;
    private final double MAX_GENOTYPE_QUALITY = 6.0;

    public void initialize() {

        // setup the header fields

        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFInfoHeaderLine("R2", 1, VCFHeaderLineType.Float, "r2 Value reported by Beable on each site"));
        hInfo.add(new VCFHeaderLine("source", "BeagleImputation"));

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

    }


     public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if ( tracker == null )
            return 0;

        GenomeLoc loc = context.getLocation();
        VariantContext vc_input = tracker.getVariantContext(ref,"inputvcf", null, loc, false);
        if ( vc_input == null  )
            return 0;


        List<Object> r2rods = tracker.getReferenceMetaData("beagleR2");

        // ignore places where we don't have a variant
        if ( r2rods.size() == 0 )
            return 0;

        BeagleFeature beagleR2Feature = (BeagleFeature)r2rods.get(0);

        List<Object> gProbsrods = tracker.getReferenceMetaData("beagleProbs");

        // ignore places where we don't have a variant
        if ( gProbsrods.size() == 0 )
            return 0;

        BeagleFeature beagleProbsFeature = (BeagleFeature)gProbsrods.get(0);

        List<Object> gPhasedrods = tracker.getReferenceMetaData("beaglePhased");

        // ignore places where we don't have a variant
        if ( gPhasedrods.size() == 0 )
            return 0;

        BeagleFeature beaglePhasedFeature = (BeagleFeature)gPhasedrods.get(0);

        // get reference base for current position
        byte refByte = ref.getBase();

        // make new Genotypes based on Beagle results
        Map<String, Genotype> genotypes = new HashMap<String, Genotype>(vc_input.getGenotypes().size());


        // for each genotype, create a new object with Beagle information on it
        for ( Map.Entry<String, Genotype> genotype : vc_input.getGenotypes().entrySet() ) {


            Genotype g = genotype.getValue();
            Set<String> filters = new LinkedHashSet<String>(g.getFilters());

            boolean genotypeIsPhased = true;
            String sample = g.getSampleName();

            ArrayList<String> beagleProbabilities = beagleProbsFeature.getProbLikelihoods().get(sample);
            ArrayList<String> beagleGenotypePairs = beaglePhasedFeature.getGenotypes().get(sample);


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

        VCFRecord vcf = VariantContextAdaptors.toVCF(filteredVC, ref.getBase());

        if ( filteredVC.getChromosomeCount() > 0 ) {
            vcf.addInfoField(VCFRecord.ALLELE_NUMBER_KEY, String.format("%d", filteredVC.getChromosomeCount()));
            if ( altAlleleCountString.length() > 0 )  {
                vcf.addInfoField(VCFRecord.ALLELE_COUNT_KEY, altAlleleCountString.toString());
                vcf.addInfoField(VCFRecord.ALLELE_FREQUENCY_KEY, String.format("%4.2f",
                        Double.valueOf(altAlleleCountString.toString())/(filteredVC.getChromosomeCount())));
            }
        }

        vcf.addInfoField("R2", beagleR2Feature.getR2value().toString() );
        vcfWriter.addRecord(vcf);


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

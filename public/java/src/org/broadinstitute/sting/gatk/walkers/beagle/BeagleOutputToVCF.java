/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.beagle;

import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.gatk.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.codecs.beagle.BeagleFeature;
import org.broadinstitute.sting.utils.variant.GATKVCFUtils;
import org.broadinstitute.variant.vcf.*;
import org.broadinstitute.sting.utils.help.DocumentedGATKFeature;
import org.broadinstitute.variant.variantcontext.*;
import org.broadinstitute.variant.variantcontext.writer.VariantContextWriter;

import java.util.*;

import static java.lang.Math.log10;


/**
 * Takes files produced by Beagle imputation engine and creates a vcf with modified annotations.
 *
 * <p>This walker is intended to be run after Beagle has successfully executed. The full calling sequence for using Beagle along with the GATK is:      </p>
 *
 * <p>1. Run ProduceBeagleInputWalker.  </p>
 * <p>2. Run Beagle</p>
 * <p>3. Uncompress output files</p>
 * <p>4. Run BeagleOutputToVCFWalker.</p>
 *
 *
 * Note that this walker requires all input files produced by Beagle.
 *
 *
 * <h2>Example</h2>
 * <pre>
 *     java -Xmx4000m -jar dist/GenomeAnalysisTK.jar \
 *      -R reffile.fasta -T BeagleOutputToVCF \
 *      -V input_vcf.vcf \
 *      -beagleR2:BEAGLE /myrun.beagle_output.r2 \
 *      -beaglePhased:BEAGLE /myrun.beagle_output.phased \
 *      -beagleProbs:BEAGLE /myrun.beagle_output.gprobs \
 *      -o output_vcf.vcf
 *      </pre>

 <p> Note that Beagle produces some of these files compressed as .gz, so gunzip must be run on them before walker is run in order to decompress them </p>

 */
@DocumentedGATKFeature( groupName = "Variant Discovery Tools", extraDocs = {CommandLineGATK.class} )
public class BeagleOutputToVCF extends RodWalker<Integer, Integer> {

    @ArgumentCollection
    protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

    /**
     * If this argument is present, the original allele frequencies and counts from this vcf are added as annotations ACH,AFH and ANH. at each record present in this vcf
     */
    @Input(fullName="comp", shortName = "comp", doc="Comparison VCF file", required=false)
    public RodBinding<VariantContext> comp;


    /**
     * This required argument is used to annotate each site in the vcf INFO field with R2 annotation. Will be NaN if Beagle determined there are no variant samples.
     */
    @Input(fullName="beagleR2", shortName = "beagleR2", doc="Beagle-produced .r2 file containing R^2 values for all markers", required=true)
    public RodBinding<BeagleFeature> beagleR2;

    /**
     * These values will populate the GL field for each sample and contain the posterior probability of each genotype given the data after phasing and imputation.
     */
    @Input(fullName="beagleProbs", shortName = "beagleProbs", doc="Beagle-produced .probs file containing posterior genotype probabilities", required=true)
    public RodBinding<BeagleFeature> beagleProbs;

    /**
     * By default, all genotypes will be marked in the VCF as "phased", using the "|" separator after Beagle.
     */
    @Input(fullName="beaglePhased", shortName = "beaglePhased", doc="Beagle-produced .phased file containing phased genotypes", required=true)
    public RodBinding<BeagleFeature> beaglePhased;

    @Output(doc="VCF File to which variants should be written",required=true)
    protected VariantContextWriter vcfWriter = null;

    /**
     * If this argument is absent, and if Beagle determines that there is no sample in a site that has a variant genotype, the site will be marked as filtered (Default behavior).
     * If the argument is present, the site won't be marked as filtered under this condition even if there are no variant genotypes.
     */
    @Argument(fullName="dont_mark_monomorphic_sites_as_filtered", shortName="keep_monomorphic", doc="If provided, we won't filter sites that beagle tags as monomorphic.  Useful for imputing a sample's genotypes from a reference panel" ,required=false)
    public boolean DONT_FILTER_MONOMORPHIC_SITES = false;

    /**
     * Value between 0 and 1. If the probability of getting a genotype correctly (based on the posterior genotype probabilities and the actual genotype) is below this threshold,
     * a genotype will be substitute by a no-call.
     */
    @Argument(fullName="no" +
            "call_threshold", shortName="ncthr", doc="Threshold of confidence at which a genotype won't be called", required=false)
    private double noCallThreshold = 0.0;

    protected static String line = null;

    private final double MIN_PROB_ERROR = 0.000001;
    private final double MAX_GENOTYPE_QUALITY = -6.0;

    public void initialize() {

        // setup the header fields

        final Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(GATKVCFUtils.getHeaderFields(getToolkit()));
        hInfo.add(new VCFFormatHeaderLine("OG",1, VCFHeaderLineType.String, "Original Genotype input to Beagle"));
        hInfo.add(new VCFInfoHeaderLine("R2", 1, VCFHeaderLineType.Float, "r2 Value reported by Beagle on each site"));
        hInfo.add(new VCFInfoHeaderLine("NumGenotypesChanged", 1, VCFHeaderLineType.Integer, "The number of genotypes changed by Beagle"));
        hInfo.add(new VCFFilterHeaderLine("BGL_RM_WAS_A", "This 'A' site was set to monomorphic by Beagle"));
        hInfo.add(new VCFFilterHeaderLine("BGL_RM_WAS_C", "This 'C' site was set to monomorphic by Beagle"));
        hInfo.add(new VCFFilterHeaderLine("BGL_RM_WAS_G", "This 'G' site was set to monomorphic by Beagle"));
        hInfo.add(new VCFFilterHeaderLine("BGL_RM_WAS_T", "This 'T' site was set to monomorphic by Beagle"));

        if ( comp.isBound() ) {
            hInfo.add(new VCFInfoHeaderLine("ACH", 1, VCFHeaderLineType.Integer, "Allele Count from Comparison ROD at this site"));
            hInfo.add(new VCFInfoHeaderLine("ANH", 1, VCFHeaderLineType.Integer, "Allele Frequency from Comparison ROD at this site"));
            hInfo.add(new VCFInfoHeaderLine("AFH", 1, VCFHeaderLineType.Float, "Allele Number from Comparison ROD at this site"));
        }

        Set<String> samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(variantCollection.variants.getName()));

        final VCFHeader vcfHeader = new VCFHeader(hInfo, samples);
        vcfWriter.writeHeader(vcfHeader);
    }

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if ( tracker == null )
            return 0;

        GenomeLoc loc = context.getLocation();
        VariantContext vc_input = tracker.getFirstValue(variantCollection.variants, loc);

        VariantContext vc_comp = tracker.getFirstValue(comp, loc);

        if ( vc_input == null  )
            return 0;

        if (vc_input.isFiltered()) {
            vcfWriter.add(vc_input);
            return 1;
        }

        BeagleFeature beagleR2Feature = tracker.getFirstValue(beagleR2);
        BeagleFeature beagleProbsFeature = tracker.getFirstValue(beagleProbs);
        BeagleFeature beaglePhasedFeature = tracker.getFirstValue(beaglePhased);

        // ignore places where we don't have a variant
        if ( beagleR2Feature == null || beagleProbsFeature == null ||  beaglePhasedFeature == null)
        {
            vcfWriter.add(vc_input);
            return 1;
        }


        // get reference base for current position
        byte refByte = ref.getBase();

        // make new Genotypes based on Beagle results
        GenotypesContext genotypes = GenotypesContext.create(vc_input.getGenotypes().size());

        // for each genotype, create a new object with Beagle information on it

        int numGenotypesChangedByBeagle = 0;
        Integer alleleCountH = 0, chrCountH = 0;
        Double alleleFrequencyH = 0.0;
        int beagleVarCounts = 0;

        GenotypesContext hapmapGenotypes = null;

        if (vc_comp != null) {
            hapmapGenotypes = vc_comp.getGenotypes();
        }

        for ( final Genotype g : vc_input.getGenotypes() ) {
            boolean genotypeIsPhased = true;
            String sample = g.getSampleName();

            // If we have  a Hapmap (comp) ROD, compute Hapmap AC, AN and AF
            // use sample as key into genotypes structure
            if (vc_comp != null) {

                if (vc_input.getGenotypes().containsSample(sample) && hapmapGenotypes.containsSample(sample))  {

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

            // original alleles at this genotype
            Allele originalAlleleA = g.getAllele(0);

            Allele originalAlleleB = (g.getAlleles().size() == 2) ? g.getAllele(1) : g.getAllele(0); // hack to deal with no-call genotypes


            // We have phased genotype in hp. Need to set the isRef field in the allele.
            List<Allele> alleles = new ArrayList<Allele>();

            String alleleA = beagleGenotypePairs.get(0);
            String alleleB = beagleGenotypePairs.get(1);

            if ( alleleA.equals("null") || alleleB.equals("null") ) {
                logger.warn("Beagle produced 'null' alleles at location "+ref.getLocus().toString()+". Ignoring.");
                return 0;
            }

            // Beagle always produces genotype strings based on the strings we input in the likelihood file.
            String refString = vc_input.getReference().getDisplayString();

            Allele bglAlleleA, bglAlleleB;

            if (alleleA.matches(refString))
                bglAlleleA = Allele.create(alleleA,true);
            else
                bglAlleleA = Allele.create(alleleA,false);

            if (alleleB.matches(refString))
                bglAlleleB = Allele.create(alleleB,true);
            else
                bglAlleleB = Allele.create(alleleB,false);


            alleles.add(bglAlleleA);
            alleles.add(bglAlleleB);

            // Compute new GQ field = -10*log10Pr(Genotype call is wrong)
            // Beagle gives probability that genotype is AA, AB and BB.
            // Which, by definition, are prob of hom ref, het and hom var.
            double probWrongGenotype, genotypeQuality;
            Double homRefProbability = Double.valueOf(beagleProbabilities.get(0));
            Double hetProbability = Double.valueOf(beagleProbabilities.get(1));
            Double homVarProbability = Double.valueOf(beagleProbabilities.get(2));

            if (bglAlleleA.isReference() && bglAlleleB.isReference()) // HomRef call
                probWrongGenotype = hetProbability + homVarProbability;
            else if ((bglAlleleB.isReference() && bglAlleleA.isNonReference()) || (bglAlleleA.isReference() && bglAlleleB.isNonReference()))
                probWrongGenotype = homRefProbability + homVarProbability;
            else // HomVar call
                probWrongGenotype = hetProbability + homRefProbability;

            // deal with numerical errors coming from limited formatting value on Beagle output files
            if (probWrongGenotype > 1 - MIN_PROB_ERROR)
                probWrongGenotype = 1 - MIN_PROB_ERROR;

            if (1-probWrongGenotype < noCallThreshold) {
                // quality is bad: don't call genotype
                alleles.clear();
                alleles.add(originalAlleleA);
                alleles.add(originalAlleleB);
                genotypeIsPhased = false;
            }

            if (probWrongGenotype < MIN_PROB_ERROR)
                genotypeQuality = MAX_GENOTYPE_QUALITY;
            else
                genotypeQuality = log10(probWrongGenotype);

            HashMap<String,Object> originalAttributes = new HashMap<String,Object>(g.getExtendedAttributes());

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
            if (! originalAlleleA.equals(Allele.NO_CALL) && beagleSwitchedGenotypes(bglAlleleA,originalAlleleA,bglAlleleB,originalAlleleB)){
                originalAttributes.put("OG",og);
                numGenotypesChangedByBeagle++;
            }
            else {
                originalAttributes.put("OG",".");
            }
            Genotype imputedGenotype = new GenotypeBuilder(g).alleles(alleles).log10PError(genotypeQuality).attributes(originalAttributes).phased(genotypeIsPhased).make();
            if ( imputedGenotype.isHet() || imputedGenotype.isHomVar() ) {
                beagleVarCounts++;
            }

            genotypes.add(imputedGenotype);
        }

        final VariantContextBuilder builder = new VariantContextBuilder(vc_input).source("outputvcf").genotypes(genotypes);
        if ( ! ( beagleVarCounts > 0 || DONT_FILTER_MONOMORPHIC_SITES ) ) {
            Set<String> removedFilters = vc_input.filtersWereApplied() ? new HashSet<String>(vc_input.getFilters()) : new HashSet<String>(1);
            removedFilters.add(String.format("BGL_RM_WAS_%s",vc_input.getAlternateAllele(0)));
            builder.alleles(new HashSet<Allele>(Arrays.asList(vc_input.getReference()))).filters(removedFilters);
        }

        // re-compute chromosome counts
        VariantContextUtils.calculateChromosomeCounts(builder, false);

        // Get Hapmap AC and AF
        if (vc_comp != null) {
            builder.attribute("ACH", alleleCountH.toString() );
            builder.attribute("ANH", chrCountH.toString() );
            builder.attribute("AFH", String.format("%4.2f", (double)alleleCountH/chrCountH) );

        }

        builder.attribute("NumGenotypesChanged", numGenotypesChangedByBeagle );
        if( !beagleR2Feature.getR2value().equals(Double.NaN) ) {
            builder.attribute("R2", beagleR2Feature.getR2value().toString() );
        }

        vcfWriter.add(builder.make());

        return 1;
    }

    private boolean beagleSwitchedGenotypes(Allele bglAlleleA, Allele originalAlleleA, Allele bglAlleleB, Allele originalAlleleB) {
       return !((bglAlleleA.equals(originalAlleleA) && bglAlleleB.equals(originalAlleleB) ||
                    (bglAlleleA.equals(originalAlleleB) && bglAlleleB.equals(originalAlleleA))));
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
        System.out.printf("Processed %d loci.\n", result);
    }
}

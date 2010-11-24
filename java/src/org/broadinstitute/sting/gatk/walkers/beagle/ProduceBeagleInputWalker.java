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

package org.broadinstitute.sting.gatk.walkers.beagle;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.*;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Hidden;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.vcf.VCFUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Produces an input file to Beagle imputation engine, listing genotype likelihoods for each sample in input variant file
 */
@Requires(value={},referenceMetaData=@RMD(name=ProduceBeagleInputWalker.ROD_NAME, type=VariantContext.class))
public class ProduceBeagleInputWalker extends RodWalker<Integer, Integer> {

    public static final String ROD_NAME = "variant";
    public static final String VALIDATION_ROD_NAME = "validation";

    @Output(doc="File to which BEAGLE input should be written",required=true)
    protected PrintStream  beagleWriter = null;

    @Argument(fullName = "genotypes_file", shortName = "genotypes", doc = "File to print reference genotypes for error analysis", required = false)
    public PrintStream beagleGenotypesWriter  = null;
    @Argument(fullName = "inserted_nocall_rate", shortName = "nc_rate", doc = "Rate (0-1) at which genotype no-calls will be randomly inserted, for testing", required = false)
    public double insertedNoCallRate  = 0;
    @Argument(fullName = "validation_genotype_ptrue", shortName = "valp", doc = "Flat probability to assign to validation genotypes. Will override GL field.", required = false)
    public double validationPrior = -1.0;
    @Argument(fullName = "validation_bootstrap", shortName = "bs", doc = "Proportion of records to be used in bootstrap set", required = false)
    public double bootstrap = 0.0;
    @Argument(fullName = "bootstrap_vcf",shortName = "bvcf", doc = "Output a VCF with the records used for bootstrapping filtered out", required = false)
    VCFWriter bootstrapVCFOutput = null;
    @Argument(fullName = "checkIsMaleOnChrX", shortName = "checkIsMaleOnChrX", doc = "Set to true when Beagle-ing chrX and want to ensure male samples don't have heterozygous calls.", required = false)
    public boolean CHECK_IS_MALE_ON_CHR_X = false;

    @Hidden
    @Argument(fullName = "variant_genotype_ptrue", shortName = "varp", doc = "Flat probability prior to assign to variant (not validation) genotypes. Does not override GL field.", required = false)
    public double variantPrior = 0.96;


    private Set<String> samples = null;
    private Set<String> BOOTSTRAP_FILTER = new HashSet<String>( Arrays.asList("bootstrap") );
    private Random generator;
    private int bootstrapCount = -1;
    private int bootstrapSetSize = 0;
    private int testSetSize = 0;

    public void initialize() {
        generator = new Random();

        samples = SampleUtils.getSampleListWithVCFHeader(getToolkit(), Arrays.asList(ROD_NAME));

        beagleWriter.print("marker alleleA alleleB");
        for ( String sample : samples )
            beagleWriter.print(String.format(" %s %s %s", sample, sample, sample));

        beagleWriter.println();
        if (beagleGenotypesWriter != null)
            beagleGenotypesWriter.println("dummy header");

        if ( bootstrapVCFOutput != null ) {
            initializeVcfWriter();
        }
    }

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker != null ) {
            GenomeLoc loc = context.getLocation();
            VariantContext variant_eval;
            VariantContext validation_eval;

            variant_eval = tracker.getVariantContext(ref, ROD_NAME, null, loc, true);
            validation_eval = tracker.getVariantContext(ref,VALIDATION_ROD_NAME,null,loc, true);
            if ( goodSite(variant_eval,validation_eval) ) {
                if ( useValidation(variant_eval,validation_eval, ref) ) {
                    writeBeagleOutput(validation_eval,variant_eval,true,validationPrior);
                } else {
                    if ( goodSite(variant_eval) ) {
                        writeBeagleOutput(variant_eval,validation_eval,false,variantPrior);
                        return 1;
                    } else { // todo -- if the variant site is bad, validation is good, but not in bootstrap set -- what do?
                        return 0;
                    }
                }
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }
    }

    public boolean goodSite(VariantContext a, VariantContext b) {
        return goodSite(a) || goodSite(b);
    }

    public boolean goodSite(VariantContext v) {
       return v != null && ! v.isFiltered() && v.isBiallelic() && v.hasGenotypes();
    }

    public boolean useValidation(VariantContext variant, VariantContext validation, ReferenceContext ref) {
        if( goodSite(validation) ) {
            // if using record keeps us below expected proportion, use it
            logger.debug(String.format("boot: %d, test: %d, total: %d", bootstrapSetSize, testSetSize, bootstrapSetSize+testSetSize+1));
            if ( (bootstrapSetSize+1.0)/(1.0+bootstrapSetSize+testSetSize) <= bootstrap ) {
                if ( bootstrapVCFOutput != null ) {
                    bootstrapVCFOutput.add(VariantContext.modifyFilters(validation, BOOTSTRAP_FILTER), ref.getBase() );
                }
                bootstrapSetSize ++;
                return true;
            } else {
                if ( bootstrapVCFOutput != null ) {
                    bootstrapVCFOutput.add(validation,ref.getBase());
                }
                testSetSize++;
                return false;
            }
        } else {
            if ( validation != null && bootstrapVCFOutput != null ) {
                bootstrapVCFOutput.add(validation,ref.getBase());
            }
            return false;
        }
    }

    public void writeBeagleOutput(VariantContext preferredVC, VariantContext otherVC, boolean isValidationSite, double prior) {
        GenomeLoc currentLoc = VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(),preferredVC);
        beagleWriter.print(String.format("%s:%d ",currentLoc.getContig(),currentLoc.getStart()));
        if ( beagleGenotypesWriter != null ) {
            beagleGenotypesWriter.print(String.format("%s ",VariantContextUtils.getLocation(getToolkit().getGenomeLocParser(),preferredVC).toString()));
        }

        for ( Allele allele : preferredVC.getAlleles() ) {
            String bglPrintString;
            if (allele.isNoCall() || allele.isNull())
                bglPrintString = "-";
            else
                bglPrintString = allele.getBaseString();  // get rid of * in case of reference allele

            beagleWriter.print(String.format("%s ", bglPrintString));
        }

        Map<String,Genotype> preferredGenotypes = preferredVC.getGenotypes();
        Map<String,Genotype> otherGenotypes = goodSite(otherVC) ? otherVC.getGenotypes() : null;
        boolean isValidation;
        for ( String sample : samples ) {
            boolean isMaleOnChrX = false;
            if( CHECK_IS_MALE_ON_CHR_X && getToolkit().getSampleById(sample).isMale() ) {
                isMaleOnChrX = true;
            }
            Genotype genotype;
            // use sample as key into genotypes structure
            if ( preferredGenotypes.keySet().contains(sample) ) {
                genotype = preferredGenotypes.get(sample);
                isValidation = isValidationSite;
            } else if ( otherGenotypes != null && otherGenotypes.keySet().contains(sample) ) {
                genotype = otherGenotypes.get(sample);
                isValidation = ! isValidationSite;
            } else {
                // there is magically no genotype for this sample.
                throw new StingException("Sample "+sample+" arose with no genotype in variant or validation VCF. This should never happen.");
            }
            /**
             * Use likelihoods if: is validation, prior is negative; or: is not validation, has genotype key
             */
            if ( (isValidation && prior < 0.0) || genotype.isCalled() && genotype.hasLikelihoods()) {
                double[] likeArray = genotype.getLikelihoods().getAsVector();
                if( isMaleOnChrX ) {
                    likeArray[1] = -255;
                }
                double[] normalizedLikelihoods = MathUtils.normalizeFromLog10(likeArray);
                // see if we need to randomly mask out genotype in this position.
                Double d = generator.nextDouble();
                if (d > insertedNoCallRate ) {
//                    System.out.format("%5.4f ", d);
                    for (Double likeVal: normalizedLikelihoods)
                        beagleWriter.print(String.format("%5.4f ",likeVal));
                }
                else {
                    // we are masking out this genotype
                    if( isMaleOnChrX ) {
                        beagleWriter.print("0.5 0.0 0.5 ");
                    } else {
                        beagleWriter.print("0.33 0.33 0.33 ");
                    }
                }

                if (beagleGenotypesWriter != null) {
                    char a = genotype.getAllele(0).toString().charAt(0);
                    char b = genotype.getAllele(0).toString().charAt(0);

                    beagleGenotypesWriter.format("%c %c ", a, b);
                }
            }
            /**
             * otherwise, use the prior uniformly
             */
            else if (! isValidation && genotype.isCalled() && ! genotype.hasLikelihoods() ) {
                // hack to deal with input VCFs with no genotype likelihoods.  Just assume the called genotype
                // is confident.  This is useful for Hapmap and 1KG release VCFs.
                double AA = (1.0-prior)/2.0;
                double AB = (1.0-prior)/2.0;
                double BB = (1.0-prior)/2.0;

                if (genotype.isHomRef()) { AA = prior; }
                else if (genotype.isHet()) { AB = prior; }
                else if (genotype.isHomVar()) { BB = prior; }

                if( isMaleOnChrX ) {
                    beagleWriter.printf("%.2f %.2f %.2f ", AA, 0.0, BB);
                } else {
                    beagleWriter.printf("%.2f %.2f %.2f ", AA, AB, BB);
                }

                if (beagleGenotypesWriter != null) {
                    char a = genotype.getAllele(0).toString().charAt(0);
                    char b = genotype.getAllele(0).toString().charAt(0);

                    beagleGenotypesWriter.format("%c %c ", a, b);
                }
            }
            else  {
                if( isMaleOnChrX ) {
                    beagleWriter.print("0.5 0.0 0.5 ");
                } else {
                    beagleWriter.print("0.33 0.33 0.33 ");
                } // write 1/3 likelihoods for uncalled genotypes.
                if (beagleGenotypesWriter != null)
                    beagleGenotypesWriter.print(". . ");
            }
        }

        beagleWriter.println();
        if (beagleGenotypesWriter != null)
            beagleGenotypesWriter.println();
    }


    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    public Integer reduce( Integer value, Integer sum ) {
        return 0; // Nothing to do here
    }

    public void onTraversalDone( Integer sum ) {

    }

    private void initializeVcfWriter() {

        final ArrayList<String> inputNames = new ArrayList<String>();
        inputNames.add( VALIDATION_ROD_NAME );

        // setup the header fields
        Set<VCFHeaderLine> hInfo = new HashSet<VCFHeaderLine>();
        hInfo.addAll(VCFUtils.getHeaderFields(getToolkit(), inputNames));
        hInfo.add(new VCFFilterHeaderLine("bootstrap","This site used for genotype bootstrapping with ProduceBeagleInputWalker"));

        bootstrapVCFOutput.writeHeader(new VCFHeader(hInfo, SampleUtils.getUniqueSamplesFromRods(getToolkit(), inputNames)));
    }

}

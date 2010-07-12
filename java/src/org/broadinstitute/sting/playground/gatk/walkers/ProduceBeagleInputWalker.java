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

import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Allele;
import org.broadinstitute.sting.gatk.contexts.variantcontext.Genotype;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.MathUtils;

import java.io.PrintStream;
import java.util.*;

/**
 * Produces an input file to Beagle imputation engine, listing genotype likelihoods for each sample in input variant file
 */
@Requires(value={},referenceMetaData=@RMD(name=ProduceBeagleInputWalker.ROD_NAME, type=ReferenceOrderedDatum.class))
public class ProduceBeagleInputWalker extends RodWalker<Integer, Integer> {

    public static final String ROD_NAME = "variant";

    @Argument(fullName = "beagle_file", shortName = "beagle", doc = "File to print BEAGLE-specific data for use with imputation", required = true)
    public PrintStream  beagleWriter = null;
    @Argument(fullName = "genotypes_file", shortName = "genotypes", doc = "File to print reference genotypes for error analysis", required = false)
    public PrintStream beagleGenotypesWriter  = null;
    @Argument(fullName = "inserted_nocall_rate", shortName = "nc_rate", doc = "Rate (0-1) at which genotype no-calls will be randomly inserted, for testing", required = false)
    public double insertedNoCallRate  = 0;

    private TreeSet<String> samples = null;
    private Random generator;

    private void initialize(Set<String> sampleNames) {
        generator = new Random();
        beagleWriter.print("marker alleleA alleleB");
        samples = new TreeSet<String>(sampleNames);
        for ( String sample : samples )
            beagleWriter.print(String.format(" %s %s %s", sample, sample, sample));

        beagleWriter.println();
        if (beagleGenotypesWriter != null)
            beagleGenotypesWriter.println("dummy header");
    }

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

         if( tracker != null ) {
            GenomeLoc loc = context.getLocation();
            VariantContext vc_eval;

            vc_eval = tracker.getVariantContext(ref, ROD_NAME, null, loc, false);
            if ( vc_eval == null || vc_eval.isFiltered() )
                return 0;

            if ( samples == null ) {
                initialize(vc_eval.getSampleNames());
            }

            // output marker ID to Beagle input file
            beagleWriter.print(String.format("%s ", vc_eval.getLocation().toString()));

            if (beagleGenotypesWriter != null)
                beagleGenotypesWriter.print(String.format("%s ", vc_eval.getLocation().toString()));

            for (Allele allele: vc_eval.getAlleles()) {
                // TODO -- check whether this is really needed by Beagle
                String bglPrintString;
                if (allele.isNoCall() || allele.isNull())
                    bglPrintString = "0";
                else
                    bglPrintString = allele.toString().substring(0,1);  // get rid of * in case of reference allele

                beagleWriter.print(String.format("%s ", bglPrintString));
            }

            if ( !vc_eval.hasGenotypes() )
                return 0;

            Map<String, Genotype> genotypes = vc_eval.getGenotypes();
            for ( String sample : samples ) {
                // use sample as key into genotypes structure
                Genotype genotype = genotypes.get(sample);
                if (genotype.isCalled() && genotype.hasAttribute(VCFConstants.GENOTYPE_LIKELIHOODS_KEY)) {
                    String[] glArray = genotype.getAttributeAsString(VCFConstants.GENOTYPE_LIKELIHOODS_KEY).split(",");

                    Double maxLikelihood = -100.0;
                    ArrayList<Double> likeArray = new ArrayList<Double>();

                    for (String gl : glArray) {
                        // need to normalize likelihoods to avoid precision loss. In worst case, if all 3 log-likelihoods are too
                        // small, we could end up with linear likelihoods of form 0.00 0.00 0.00 which will mess up imputation.
                        Double dg = Double.valueOf(gl);
                        if (dg> maxLikelihood)
                            maxLikelihood = dg;

                        likeArray.add(dg);
                    }

                    // see if we need to randomly mask out genotype in this position.
                    Double d = generator.nextDouble();
                    if (d > insertedNoCallRate ) {
//                    System.out.format("%5.4f ", d);
                        for (Double likeVal: likeArray)
                            beagleWriter.print(String.format("%5.4f ",Math.pow(10, likeVal-maxLikelihood)));
                    }
                    else {
                        // we are masking out this genotype
                        beagleWriter.print("0.33 0.33 0.33 ");
                    }

                    if (beagleGenotypesWriter != null) {
                        char a = genotype.getAllele(0).toString().charAt(0);
                        char b = genotype.getAllele(0).toString().charAt(0);

                        beagleGenotypesWriter.format("%c %c ", a, b);
                    }
                }
                else  {
                    beagleWriter.print("0.33 0.33 0.33 "); // write 1/3 likelihoods for uncalled genotypes.
                    if (beagleGenotypesWriter != null)
                         beagleGenotypesWriter.print(". . ");
                }
            }

            beagleWriter.println();
            if (beagleGenotypesWriter != null)
                 beagleGenotypesWriter.println();
        }
        return 1;
        
    }

    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    public Integer reduce( Integer value, Integer sum ) {
        return 0; // Nothing to do here
    }

    public void onTraversalDone( Integer sum ) {

   }
    
}

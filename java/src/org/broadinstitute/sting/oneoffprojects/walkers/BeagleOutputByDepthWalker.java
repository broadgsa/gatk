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

package org.broadinstitute.sting.oneoffprojects.walkers;

import org.broad.tribble.util.variantcontext.Allele;
import org.broad.tribble.util.variantcontext.Genotype;
import org.broad.tribble.util.variantcontext.VariantContext;
import org.broad.tribble.vcf.VCFConstants;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.commandline.Output;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.GenomeLoc;

import java.io.PrintStream;
import java.util.*;

/**
 * Produces an input file to Beagle imputation engine, listing genotype likelihoods for each sample in input VCF file
 * @help.summary Produces an input file to Beagle imputation engine, listing genotype likelihoods for each sample in input VCF file
 */
public class BeagleOutputByDepthWalker extends RodWalker<Integer, Integer> {


    public static final String POSTBEAGLE_EVAL_ROD_NAME = "postbeaglevcf";
    public static final String PREBEAGLE_EVAL_ROD_NAME = "prebeaglevcf";
    public static final String INPUT_HAPMAP_ROD_NAME = "hapmap";
    public static final String INPUT_COMP_ROD_NAME = "comp";

    @Output
    protected PrintStream outputWriter = null;


    public void initialize() {

    }


    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker != null ) {
            GenomeLoc loc = context.getLocation();

            VariantContext vc_postbgl = tracker.getVariantContext(ref,POSTBEAGLE_EVAL_ROD_NAME, null, loc, false);
            VariantContext vc_prebgl = tracker.getVariantContext(ref,PREBEAGLE_EVAL_ROD_NAME, null, loc, false);
            VariantContext vc_hapmap = tracker.getVariantContext(ref,INPUT_HAPMAP_ROD_NAME, null, loc, false);
            VariantContext vc_comp = tracker.getVariantContext(ref,INPUT_COMP_ROD_NAME, null, loc, false);
            if ( vc_postbgl == null || vc_prebgl == null || vc_comp == null)
                return 0;


            if (!vc_prebgl.hasGenotypes() || !vc_postbgl.hasGenotypes() )
                return 0;

            if (vc_postbgl.isFiltered())
                return 0;

            Map<String, Genotype> compGenotypes = vc_comp.getGenotypes();
            Integer alleleCountH = 0, chrCountH = 0, alleleCountEmp=0, chrCountEmp;

            // Get Hapmap AC and AF
            if (vc_hapmap != null) {
                Map<String, Genotype> hapmapGenotypes = vc_hapmap.getGenotypes();
                for ( String sample : vc_postbgl.getSampleNames() ) {
                    // use sample as key into genotypes structure
                    if (vc_postbgl.getGenotypes().containsKey(sample) && hapmapGenotypes.containsKey(sample))  {

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
            }
            else {
                alleleCountH = -1;
                chrCountH = -1;
            }

            chrCountEmp =  vc_postbgl.getChromosomeCount();
//System.out.println(chrCountH);

            if ( vc_postbgl.getAlternateAlleles().size() > 0 ) {
                for ( Allele allele : vc_postbgl.getAlternateAlleles() ) {
                    alleleCountEmp = alleleCountEmp+vc_postbgl.getChromosomeCount(allele);
                }
//System.out.println(alleleCountH);
            }




            for ( String sample : vc_postbgl.getSampleNames() ) {
                if (sample.compareToIgnoreCase("NA12878")!=0)
                    continue;

                // use sample as key into genotypes structure

                Genotype postbglGenotype = vc_postbgl.getGenotype(sample);
                Genotype prebglGenotype = vc_prebgl.getGenotype(sample);
                Genotype compGenotype = compGenotypes.get(sample);


                outputWriter.format("%d %d %d %d %d ", vc_postbgl.getStart(), alleleCountH, chrCountH,
                        alleleCountEmp, chrCountEmp);


                String dps = postbglGenotype.getAttributeAsString(VCFConstants.DEPTH_KEY);

                int dp;
                if (dps.compareTo(".")==0)
                    dp = -1;
                else
                    dp = Integer.valueOf(dps);

                int hg, bg, pg;
                if (compGenotype.isNoCall())
                    hg = -1;
                else if (compGenotype.isHomRef())
                    hg = 0;
                else if (compGenotype.isHet())
                    hg = 1;
                else if (compGenotype.isHomVar())
                    hg = 2;
                else
                    throw new ReviewedStingException("Bug! invalid genotype!");

                if (postbglGenotype.isNoCall())
                    bg = -1;
                else if (postbglGenotype.isHomRef())
                    bg = 0;
                else if (postbglGenotype.isHet())
                    bg = 1;
                else if (postbglGenotype.isHomVar())
                    bg = 2;
                else
                    throw new ReviewedStingException("Bug! invalid genotype!");


                if (prebglGenotype.isNoCall())
                    pg = -1;
                else if (prebglGenotype.isHomRef())
                    pg = 0;
                else if (prebglGenotype.isHet())
                    pg = 1;
                else if (prebglGenotype.isHomVar())
                    pg = 2;
                else
                    throw new ReviewedStingException("Bug! invalid genotype!");

                outputWriter.format("%d %d %d %d\n",dp, hg, pg, bg);


            }
            return 1;
        }
        return 0;
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


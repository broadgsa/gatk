package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.refdata.RodGLF;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.genotype.VariantBackedByGenotype;

import java.util.HashMap;

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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * Takes variant calls as .vcf files and creates plots of truth metrics as a function of the various annotations found in the INFO field.
 *
 * @author rpoplin
 * @since Jan 15, 2010
 * @help.summary Takes variant calls as .vcf files and creates plots of truth metrics as a function of the various annotations found in the INFO field. 
 */

public class AnalyzeAnnotationsWalker extends RodWalker<Integer, Integer> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName = "output_prefix", shortName = "output", doc = "The output path and name to prepend to all plots and intermediate data files", required = false)
    private String OUTPUT_PREFIX = "analyzeAnnotations/";
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is probably /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";
    @Argument(fullName = "min_variants_per_bin", shortName = "minBinSize", doc = "The minimum number of variants in a bin in order to calculate truth metrics.", required = false)
    private int MIN_VARIANTS_PER_BIN = 1000;
    @Argument(fullName = "max_variants_per_bin", shortName = "maxBinSize", doc = "The maximum number of variants in a bin.", required = false)
    private int MAX_VARIANTS_PER_BIN = 20000;
    @Argument(fullName = "sampleName", shortName = "sampleName", doc = "If supplied, only process variants found in this sample.", required = false)
    private String SAMPLE_NAME = null;
    @Argument(fullName = "name", shortName = "name", doc = "Labels for the annotations to make plots look nicer. Each name is a separate -name argument. For example, -name DP,Depth -name AB,AlleleBalance", required = false)
    private String[] ANNOTATION_NAMES = null;
    @Argument(fullName = "indicate_mean_num_vars", shortName = "meanNumVars", doc = "If supplied, plots will indicate the distribution of number of variants instead of distribution of value of annotation", required = false)
    private boolean INDICATE_MEAN_NUM_VARS = false;
        
    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private AnnotationDataManager dataManager;
    
    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {

        // Create a HashMap associating the names of the annotations to full Strings that can be used as labels on plots
        HashMap<String,String> nameMap = null;
        if( ANNOTATION_NAMES != null ) {
            nameMap = new HashMap<String,String>();
            for( String nameLine : ANNOTATION_NAMES ) {
                String[] vals = nameLine.split(",");
                nameMap.put(vals[0],vals[1]);
            }
        }
        dataManager = new AnnotationDataManager( nameMap, INDICATE_MEAN_NUM_VARS );

        if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        if( tracker != null ) {

            // First find out if this variant is in the truth sets
            boolean isInTruthSet = false;
            boolean isTrueVariant = false;
            for( final ReferenceOrderedDatum rod : tracker.getAllRods() ) {
                if( rod != null && rod.getName().toUpperCase().startsWith("TRUTH") ) {
                    isInTruthSet = true;

                    // Next see if the truth sets say this site is variant or reference
                    if( rod instanceof RodVCF ) {
                        if( ((RodVCF) rod).isSNP() ) {
                            isTrueVariant = true;
                        }
                    } else if( rod instanceof RodGLF ) {
                        if( ((RodGLF) rod).isSNP() ) {
                            isTrueVariant = true;
                        }
                    } else if( rod instanceof VariantBackedByGenotype ) {
                        if( ((VariantBackedByGenotype) rod).getCalledGenotype().isVariant(ref.getBase()) ) {
                            isTrueVariant = true;
                        }
                    } else {
                        throw new StingException( "Truth ROD is of unknown ROD type: " + rod.getName() );
                    }
                }
            }
            
            // Add each annotation in this VCF Record to the dataManager
            for( final ReferenceOrderedDatum rod : tracker.getAllRods() ) {
                if( rod != null && rod instanceof RodVCF && !rod.getName().toUpperCase().startsWith("TRUTH") ) {
                    final RodVCF variant = (RodVCF) rod;
                    if( variant.isSNP() ) {
                        dataManager.addAnnotations( variant, SAMPLE_NAME, isInTruthSet, isTrueVariant );
                    }
                }
            }
        }

        return 1; // This value isn't actually used for anything
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer reduceInit() {
        return 0; // Nothing to do here
    }

    public Integer reduce( Integer value, Integer sum ) {
        return 0; // Nothing to do here
    }

    public void onTraversalDone( Integer sum ) {

        // For each annotation, decide how to cut up the data, output intermediate cumulative p(true) tables, and call RScript to plot the tables
        dataManager.plotCumulativeTables(PATH_TO_RSCRIPT, PATH_TO_RESOURCES, OUTPUT_PREFIX, MIN_VARIANTS_PER_BIN, MAX_VARIANTS_PER_BIN);
    }
}
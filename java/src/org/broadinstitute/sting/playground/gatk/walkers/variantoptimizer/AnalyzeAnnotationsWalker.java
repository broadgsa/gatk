package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.IOException;

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
    @Argument(fullName = "output_dir", shortName = "outputDir", doc = "The directory in which to output all the plots and intermediate data files", required = false)
    private String OUTPUT_DIR = "analyzeAnnotations/";
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


    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final AnnotationDataManager dataManager = new AnnotationDataManager();
    
    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    /**
     * Create the output directory and setup the path variables
     */
    public void initialize() {

        // create the output directory where all the data tables and plots will go
        try {
            Process p = Runtime.getRuntime().exec("mkdir " + OUTPUT_DIR);
        } catch (IOException e) {
            throw new RuntimeException("Couldn't create directory: " + OUTPUT_DIR);
        }
        if( !OUTPUT_DIR.endsWith("/") ) { OUTPUT_DIR = OUTPUT_DIR + "/"; }
        if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public Integer map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        // Add each VCF record to each annotation's data structure
        if( tracker != null ) {
            for( ReferenceOrderedDatum rod : tracker.getAllRods() ) {
                if( rod != null && rod instanceof RodVCF ) {
                    RodVCF variant = (RodVCF) rod;
                    if( variant.isSNP() ) {
                        dataManager.addAnnotations( variant, SAMPLE_NAME );
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
        dataManager.plotCumulativeTables(PATH_TO_RSCRIPT, PATH_TO_RESOURCES, OUTPUT_DIR, MIN_VARIANTS_PER_BIN, MAX_VARIANTS_PER_BIN);
    }
}
package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.datasources.simpleDataSources.ReferenceOrderedDataSource;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.gatk.walkers.RodWalker;
import org.broadinstitute.sting.utils.ExpandingArrayList;
import org.broadinstitute.sting.utils.Pair;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;

import java.io.IOException;
import java.io.PrintStream;
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
 * Calculates variant concordance with the given truth sets and plots ROC curves for each input call set.
 * 
 * @author rpoplin
 * @since Mar 18, 2010
 *
 * @help.summary Calculates variant concordance with the given truth sets and plots ROC curves for each input call set.
 */

public class VariantConcordanceROCCurveWalker extends RodWalker<ExpandingArrayList<Pair<String,VariantDatum>>, HashMap<String,ExpandingArrayList<VariantDatum>>> {

    /////////////////////////////
    // Command Line Arguments
    /////////////////////////////
    @Argument(fullName="output_prefix", shortName="output", doc="The prefix added to output VCF file name and optimization curve pdf file name", required=false)
    private String OUTPUT_PREFIX = "optimizer";
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript. For Broad users this is probably /broad/tools/apps/R-2.6.0/bin/Rscript", required = false)
    private String PATH_TO_RSCRIPT = "/broad/tools/apps/R-2.6.0/bin/Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting R scripts.", required = false)
    private String PATH_TO_RESOURCES = "R/";

    /////////////////////////////
    // Private Member Variables
    /////////////////////////////
    private final ExpandingArrayList<String> inputRodNames = new ExpandingArrayList<String>();
    private int numCurves;
    private int[] trueNegGlobal;
    private int[] falseNegGlobal;
    private String sampleName = null;

    //---------------------------------------------------------------------------------------------------------------
    //
    // initialize
    //
    //---------------------------------------------------------------------------------------------------------------

    public void initialize() {
        if( !PATH_TO_RESOURCES.endsWith("/") ) { PATH_TO_RESOURCES = PATH_TO_RESOURCES + "/"; }
        
        for( ReferenceOrderedDataSource rod : this.getToolkit().getRodDataSources() ) {
            if( rod != null && !rod.getName().toUpperCase().startsWith("TRUTH") ) {
                if( rod.getReferenceOrderedData().iterator().next().get(0) instanceof RodVCF ) {
                    inputRodNames.add(rod.getName());
                    if( sampleName == null ) {
                        sampleName = ((RodVCF)rod.getReferenceOrderedData().iterator().next().get(0)).getSampleNames()[0]; // BUGBUG: single sample calls only for now
                    }
                }
            }
        }
        numCurves = inputRodNames.size();
        trueNegGlobal = new int[numCurves];
        falseNegGlobal = new int[numCurves];
        for( int kkk = 0; kkk < numCurves; kkk++ ) {
            trueNegGlobal[kkk] = 0;
            falseNegGlobal[kkk] = 0;
        }
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // map
    //
    //---------------------------------------------------------------------------------------------------------------

    public ExpandingArrayList<Pair<String,VariantDatum>> map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {

        final ExpandingArrayList<Pair<String,VariantDatum>> mapList = new ExpandingArrayList<Pair<String,VariantDatum>>();

        if( tracker == null ) { // For some reason RodWalkers get map calls with null trackers
            return mapList;
        }

        boolean isInTruthSet = false;
        boolean isTrueVariant = false;

        for( final VariantContext vc : tracker.getAllVariantContexts(null, context.getLocation(), false, false) ) {
            if( vc != null && vc.getName().toUpperCase().startsWith("TRUTH") ) {
                if( vc.isSNP() && !vc.isFiltered() ) {
                    if( !vc.getGenotype(sampleName).isNoCall() ) {
                        isInTruthSet = true;

                        if( !vc.getGenotype(sampleName).isHomRef() ) {
                            isTrueVariant = true;
                        }
                    }
                }
                //if( vc.isPolymorphic() ) { //BUGBUG: I don't think this is the right thing to do here, there are many polymorphic sites in the truth data because there are many samples
                //    isTrueVariant = true;
                //}
            }
        }

        if( !isInTruthSet ) {
            return mapList; // jump out if this location isn't in the truth set
        }

        int curveIndex = 0;
        for( final String inputName : inputRodNames ) {
            final VariantContext vc = tracker.getVariantContext( inputName, null, context.getLocation(), false ); // assuming single variant per track per location

            if( vc != null && vc.isPolymorphic() && !vc.isFiltered() ) {
                final VariantDatum variantDatum = new VariantDatum();
                variantDatum.qual = vc.getPhredScaledQual();
                variantDatum.isTrueVariant = isTrueVariant;
                mapList.add( new Pair<String,VariantDatum>(inputName, variantDatum) );
            } else { // Either not in the call set at all, is a monomorphic call, or was filtered out
                if( isTrueVariant ) { // ... but it is a true variant so this is a false negative call
                    falseNegGlobal[curveIndex]++;
                } else { // ... and it is not a variant site so this is a true negative call
                    trueNegGlobal[curveIndex]++;
                }
            }
            curveIndex++;
        }

        return mapList;
    }

    //---------------------------------------------------------------------------------------------------------------
    //
    // reduce
    //
    //---------------------------------------------------------------------------------------------------------------

    public HashMap<String,ExpandingArrayList<VariantDatum>> reduceInit() {
        final HashMap<String,ExpandingArrayList<VariantDatum>> init = new HashMap<String,ExpandingArrayList<VariantDatum>>();
        for( final String inputName : inputRodNames ) {
            init.put( inputName, new ExpandingArrayList<VariantDatum>() );
        }
        return init;
    }

    public HashMap<String,ExpandingArrayList<VariantDatum>> reduce( final ExpandingArrayList<Pair<String,VariantDatum>> mapValue, final HashMap<String,ExpandingArrayList<VariantDatum>> reduceSum ) {
        for( Pair<String,VariantDatum> value : mapValue ) {
            final ExpandingArrayList<VariantDatum> list = reduceSum.get(value.getFirst());
            list.add(value.getSecond());
            reduceSum.put(value.getFirst(),list);
        }
        return reduceSum;
    }

    public void onTraversalDone( HashMap<String,ExpandingArrayList<VariantDatum>> reduceSum ) {

        final int NUM_CURVES = numCurves;
        final HashMap<String, VariantDataManager> dataManagerMap = new HashMap<String, VariantDataManager>();
        for( final String inputName : inputRodNames ) {
            System.out.println("Creating data manager for: " + inputName);
            dataManagerMap.put(inputName, new VariantDataManager( reduceSum.get(inputName), null ));
        }
        reduceSum.clear(); // Don't need this ever again, clean up some memory

        final double[] minQual = new double[NUM_CURVES];
        final double[] maxQual = new double[NUM_CURVES];
        final double[] incrementQual = new double[NUM_CURVES];
        final double[] qualCut = new double[NUM_CURVES];

        final int NUM_STEPS = 200;

        int curveIndex = 0;
        for( final String inputName : inputRodNames ) {
            final VariantDataManager dataManager = dataManagerMap.get(inputName);
            minQual[curveIndex] =  dataManager.data[0].qual;
            maxQual[curveIndex] =  dataManager.data[0].qual;
            for( int iii = 1; iii < dataManager.data.length; iii++ ) {
                final double qual = dataManager.data[iii].qual;
                if( qual < minQual[curveIndex] ) { minQual[curveIndex] = qual; }
                else if( qual > maxQual[curveIndex] ) { maxQual[curveIndex] = qual; }
            }
            incrementQual[curveIndex] = (maxQual[curveIndex] - minQual[curveIndex]) / ((double)NUM_STEPS);
            qualCut[curveIndex] = minQual[curveIndex];
            curveIndex++;
        }

        final int[] truePos = new int[NUM_CURVES];
        final int[] falsePos = new int[NUM_CURVES];
        final int[] trueNeg = new int[NUM_CURVES];
        final int[] falseNeg = new int[NUM_CURVES];
     
        PrintStream outputFile;
        try {
            outputFile = new PrintStream( OUTPUT_PREFIX + ".dat" );
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + OUTPUT_PREFIX + ".dat" );
        }
        int jjj = 1;
        for( final String inputName : inputRodNames ) {
            outputFile.print(inputName + ",sensitivity" + jjj + ",specificity" + jjj + ",");
            jjj++;
        }
        outputFile.println("sentinel");

        for( int step = 0; step < NUM_STEPS; step++ ) {
            curveIndex = 0;
            for( final String inputName : inputRodNames ) {
                final VariantDataManager dataManager = dataManagerMap.get(inputName);
                truePos[curveIndex] = 0;
                falsePos[curveIndex] = 0;
                trueNeg[curveIndex] = 0;
                falseNeg[curveIndex] = 0;
                final int NUM_VARIANTS = dataManager.data.length;
                for( int iii  = 0; iii < NUM_VARIANTS; iii++ ) {
                    if( dataManager.data[iii].qual >= qualCut[curveIndex] ) { // this var is in this hypothetical call set
                        if( dataManager.data[iii].isTrueVariant ) {
                            truePos[curveIndex]++;
                        } else {
                            falsePos[curveIndex]++;
                        }
                    } else { // this var is out of this hypothetical call set
                        if( dataManager.data[iii].isTrueVariant ) {
                            falseNeg[curveIndex]++;
                        } else {
                            trueNeg[curveIndex]++;
                        }
                    }
                }
                final double sensitivity = ((double) truePos[curveIndex]) / ((double) truePos[curveIndex] + falseNegGlobal[curveIndex] + falseNeg[curveIndex]);
                final double specificity = ((double) trueNegGlobal[curveIndex] + trueNeg[curveIndex]) /
                                                ((double) falsePos[curveIndex] + trueNegGlobal[curveIndex] + trueNeg[curveIndex]);
                outputFile.print( String.format("%.8f,%.8f,%.8f,", qualCut[curveIndex], sensitivity, 1.0 - specificity) );
                qualCut[curveIndex] += incrementQual[curveIndex];
                curveIndex++;
            }
            outputFile.println("-1");
        }

        outputFile.close();
    
        // Execute Rscript command to plot the optimization curve
        // Print out the command line to make it clear to the user what is being executed and how one might modify it
        final String rScriptCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_variantROCCurve.R" + " " + OUTPUT_PREFIX + ".dat";
        System.out.println( rScriptCommandLine );

        // Execute the RScript command to plot the table of truth values
        try {
            Runtime.getRuntime().exec( rScriptCommandLine );
        } catch ( IOException e ) {
            throw new StingException( "Unable to execute RScript command: " + rScriptCommandLine );
        }
    }
}

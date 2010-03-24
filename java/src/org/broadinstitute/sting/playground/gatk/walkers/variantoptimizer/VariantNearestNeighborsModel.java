package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.utils.StingException;

import java.io.PrintStream;

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
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 1, 2010
 */

public final class VariantNearestNeighborsModel extends VariantOptimizationModel {

    private final int numKNN;

    public VariantNearestNeighborsModel( VariantDataManager _dataManager, final double _targetTITV, final int _numKNN ) {
        super( _targetTITV );
        //dataManager = _dataManager;
        numKNN = _numKNN;
    }
    
    public void run( final String outputPrefix ) {

        throw new StingException( "Nearest Neighbors model hasn't been updated yet." );
        /*
        final int numVariants = dataManager.numVariants;

        final double[] pTrueVariant = new double[numVariants];

        final VariantTree vTree = new VariantTree( numKNN );
        vTree.createTreeFromData( dataManager.data );

        System.out.println("Finished creating the kd-tree.");

        for(int iii = 0; iii < numVariants; iii++) {
            pTrueVariant[iii] = calcTruePositiveRateFromTITV( vTree.calcNeighborhoodTITV( dataManager.data[iii] ) );
        }

        PrintStream outputFile;
        try {
            outputFile = new PrintStream( outputPrefix + ".knn.optimize" );
        } catch (Exception e) {
            throw new StingException( "Unable to create output file: " + outputPrefix + ".knn.optimize" );
        }
        for(int iii = 0; iii < numVariants; iii++) {
            outputFile.print(String.format("%.4f",pTrueVariant[iii]) + ",");
            outputFile.println( (dataManager.data[iii].isTransition ? 1 : 0)
                    + "," + (dataManager.data[iii].isKnown? 1 : 0));
        }
        */
    }
}

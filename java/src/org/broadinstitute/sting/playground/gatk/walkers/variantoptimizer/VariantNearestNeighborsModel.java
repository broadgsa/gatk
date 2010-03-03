package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

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

public class VariantNearestNeighborsModel extends VariantOptimizationModel {

    public VariantNearestNeighborsModel( VariantDataManager _dataManager, final double _targetTITV ) {
       super( _dataManager, _targetTITV );
    }
    
    public double[] run() {

        final int numVariants = dataManager.numVariants;

        final double[] pTrueVariant = new double[numVariants];

        final VariantTree vTree = new VariantTree( 2000 );
        vTree.createTreeFromData( dataManager.data );

        System.out.println("Finished creating the kd-tree.");

        for(int iii = 0; iii < numVariants; iii++) {
            pTrueVariant[iii] = calcTruePositiveRateFromTITV( vTree.calcNeighborhoodTITV( dataManager.data[iii] ) );
        }

        return pTrueVariant;
    }
}

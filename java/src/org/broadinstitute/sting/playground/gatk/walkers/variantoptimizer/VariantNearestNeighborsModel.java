package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Mar 1, 2010
 */

public class VariantNearestNeighborsModel extends VariantOptimizationModel {

    public VariantNearestNeighborsModel( VariantDataManager _dataManager, final double _targetTITV ) {
       super( _dataManager, _targetTITV );
    }
    
    public double[][] run() {

        final int numVariants = dataManager.numVariants;

        final double[][] pTrueVariant = new double[1][numVariants];

        final VariantTree vTree = new VariantTree( 2000 );
        vTree.createTreeFromData( dataManager.data );

        System.out.println("Finished creating the kd-tree.");

        for(int iii = 0; iii < numVariants; iii++) {
            pTrueVariant[0][iii] = calcTruePositiveRateFromTITV( vTree.calcNeighborhoodTITV( dataManager.data[iii] ) );
        }

        return pTrueVariant;
    }
}

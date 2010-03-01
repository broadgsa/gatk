package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Feb 26, 2010
 */

public abstract class VariantOptimizationModel implements VariantOptimizationInterface {
    protected final VariantDataManager dataManager;
    protected final double targetTITV;

    public VariantOptimizationModel( VariantDataManager _dataManager, final double _targetTITV ) {
        dataManager = _dataManager;
        targetTITV = _targetTITV;
    }

    public double calcTruePositiveRateFromTITV( double titv ) {
        if( titv > targetTITV ) { titv -= 2.0f*(titv-targetTITV); }
        if( titv < 0.5f ) { titv = 0.5f; }

        return ( (titv - 0.5f) / (targetTITV - 0.5f) );
    }
}

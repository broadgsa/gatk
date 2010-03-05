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
 * Date: Feb 26, 2010
 */

public abstract class VariantOptimizationModel implements VariantOptimizationInterface {
    protected final VariantDataManager dataManager;
    protected final double targetTITV;

    public VariantOptimizationModel( VariantDataManager _dataManager, final double _targetTITV ) {
        dataManager = _dataManager;
        targetTITV = _targetTITV;
    }

    public final double calcTruePositiveRateFromTITV( double titv ) {
        if( titv > targetTITV ) { titv -= 2.0f*(titv-targetTITV); }
        if( titv < 0.5 ) { titv = 0.5; }
        return ( (titv - 0.5) / (targetTITV - 0.5) );
        //if( titv < 0.0 ) { titv = 0.0; }
        //return ( titv / targetTITV );
    }
}

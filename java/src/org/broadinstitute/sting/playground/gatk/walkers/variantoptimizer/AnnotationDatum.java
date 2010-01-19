package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.utils.StingException;

import java.util.Comparator;

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
 * Date: Jan 18, 2010
 */

public class AnnotationDatum implements Comparator<AnnotationDatum> {
    
    public final float value;
    public int numTransitions;
    public int numTransversions;

    public AnnotationDatum() {
        value = 0.0f;
        numTransitions = 0;
        numTransversions = 0;
    }

    public AnnotationDatum( float _value ) {
        value = _value;
        numTransitions = 0;
        numTransversions = 0;
    }

    final public void incrementTi() {
        numTransitions++;
    }

    final public void incrementTv() {
        numTransversions++;
    }

    final public float calcTiTv() {

        if( numTransitions < 0 || numTransversions < 0) {
            throw new StingException( "Integer overflow detected! There are too many variants piled up in one annotation bin." );
        }

        if( numTransversions == 0) {
            return 0.0f;
        }

        return ((float) numTransitions) / ((float) numTransversions);
    }

    public int compare( AnnotationDatum a1, AnnotationDatum a2 ) {
        if( a1.value < a2.value ) { return -1; }
        if( a1.value > a2.value ) { return 1; }
        return 0;
    }

    public int equals( AnnotationDatum that ) {
        if( this.value < that.value ) { return -1; }
        if( this.value > that.value ) { return 1; }
        return 0;
    }
}

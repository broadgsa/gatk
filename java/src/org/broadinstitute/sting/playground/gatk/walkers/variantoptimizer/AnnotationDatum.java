package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.utils.StingException;

import java.util.Comparator;

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

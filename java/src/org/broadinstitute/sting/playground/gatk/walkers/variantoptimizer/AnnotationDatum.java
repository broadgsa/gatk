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
    
    public float value;
    private final int[] ti;
    private final int[] tv;

    public static final int FULL_SET = 0;
    public static final int NOVEL_SET = 1;
    public static final int DBSNP_SET = 2;
    public static final int TRUTH_SET = 3;
    public static final int TRUE_POSITIVE = 4;
    private static final int NUM_SETS = 5;

    public AnnotationDatum() {

        value = 0.0f;
        ti = new int[NUM_SETS];
        tv = new int[NUM_SETS];
        for( int iii = 0; iii < NUM_SETS; iii++ ) {
            ti[iii] = 0;
            tv[iii] = 0;
        }
    }

    public AnnotationDatum( final float _value ) {

        value = _value;
        ti = new int[NUM_SETS];
        tv = new int[NUM_SETS];
        for( int iii = 0; iii < NUM_SETS; iii++ ) {
            ti[iii] = 0;
            tv[iii] = 0;
        }
    }

    final public void incrementTi( final boolean isNovelVariant, final boolean isInTruthSet, final boolean isTrueVariant ) {

        ti[FULL_SET]++;
        if( isNovelVariant ) {
            ti[NOVEL_SET]++;
        } else { // Is known, in DBsnp
            ti[DBSNP_SET]++;
        }
        if( isInTruthSet ) {
            ti[TRUTH_SET]++;
            if( isTrueVariant ) {
                ti[TRUE_POSITIVE]++;
            }
        }
    }

    final public void incrementTv( final boolean isNovelVariant, final boolean isInTruthSet, final boolean isTrueVariant ) {

        tv[FULL_SET]++;
        if( isNovelVariant ) {
            tv[NOVEL_SET]++;
        } else { // Is known, in DBsnp
            tv[DBSNP_SET]++;
        }
        if( isInTruthSet ) {
            tv[TRUTH_SET]++;
            if( isTrueVariant ) {
                tv[TRUE_POSITIVE]++;
            }
        }
    }

    final public void combine( final AnnotationDatum that ) {

        for( int iii = 0; iii < NUM_SETS; iii++ ) {
            this.ti[iii] += that.ti[iii];
            this.tv[iii] += that.tv[iii];
        }
        this.value = that.value; // Overwrite this bin's value
    }

    final public float calcTiTv( final int INDEX ) {

        if( ti[INDEX] < 0 || tv[INDEX] < 0 ) {
            throw new StingException( "Integer overflow detected! There are too many variants piled up in one annotation bin." );
        }

        if( tv[INDEX] == 0 ) { // Don't divide by zero
            return 0.0f;
        }

        return ((float) ti[INDEX]) / ((float) tv[INDEX]);
    }

    final public float calcDBsnpRate() {

        if( ti[FULL_SET] + tv[FULL_SET] == 0 ) { // Don't divide by zero
            return 0.0f;
        }

        return 100.0f * ((float) ti[DBSNP_SET] + tv[DBSNP_SET]) /
                ((float) ti[FULL_SET] + tv[FULL_SET]);
    }

    final public float calcTPrate() {

        if( ti[TRUTH_SET] + tv[TRUTH_SET] == 0 ) { // Don't divide by zero
            return 0.0f;
        }

        return 100.0f * ((float) ti[TRUE_POSITIVE] + tv[TRUE_POSITIVE]) /
                ((float) ti[TRUTH_SET] + tv[TRUTH_SET]);
    }

    final public int numVariants( final int INDEX ) {
        return ti[INDEX] + tv[INDEX];
    }

    final public void clearBin() {
        value = 0.0f;
        for( int iii = 0; iii < NUM_SETS; iii++ ) {
            ti[iii] = 0;
            tv[iii] = 0;
        }
    }

    public int compare( AnnotationDatum a1, AnnotationDatum a2 ) { // Function needed for this to be a Comparator
        if( a1.value < a2.value ) { return -1; }
        if( a1.value > a2.value ) { return 1; }
        return 0;
    }

    public int equals( AnnotationDatum that ) { // Function needed for this to be sorted correctly in a TreeSet
        if( this.value < that.value ) { return -1; }
        if( this.value > that.value ) { return 1; }
        return 0;
    }
}

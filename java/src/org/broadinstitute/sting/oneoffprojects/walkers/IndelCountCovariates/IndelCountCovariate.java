package org.broadinstitute.sting.oneoffprojects.walkers.IndelCountCovariates;


import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
//g import org.broadinstitute.sting.gatk.walkers.recalibration.ExperimentalCovariate;
//g import org.broadinstitute.sting.gatk.walkers.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.sting.utils.BaseUtils;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * Date: Nov 3, 2009
 *
 * The Reported Quality Score covariate.
 */

public class IndelCountCovariate implements ExperimentalCovariate {

    // Initialize any member variables using the command-line arguments passed to the walkers
    public void initialize( final RecalibrationArgumentCollection RAC ) {
    }


    // Used to pick out the covariate's value from attributes of the read
    public final Comparable getValue( final SAMRecord read ) {

        Cigar c = read.getCigar();

        int indelCount = 0;

        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
             CigarElement ce = c.getCigarElement(i);
             switch( ce.getOperator() ) {
                 case D:
                     indelCount++;
                     break;
                 case I:
                     indelCount++;
                     break;
                 case M:
                 case N:
                 case S:
                 default:
                     break;
             }
         }

        return indelCount;
    }

    public void getValues(SAMRecord read, Comparable[] comparable) {
/*        Comparable numIndels = getValue(read);
         for(int iii = 0; iii < read.getReadLength(); iii++) {
             comparable[iii] = numIndels; // BUGBUG: this can be optimized
         }    */
        int ind = 0;
        Cigar c = read.getCigar();
        System.out.println(c.toString());
        for ( int i = 0 ; i < c.numCigarElements() ; i++ ) {
              CigarElement ce = c.getCigarElement(i);

            switch( ce.getOperator() ) {
                case D:
                case I:
                    for (int k=0; k < ce.getLength(); k++)
                        comparable[ind++] = ce.getLength();
                    break;
                case M:
                case N:
                case S:
                    for (int k=0; k < ce.getLength(); k++)
                        comparable[ind++] = 0;
                default:
                    break;
            }
        }
     }


    // Used to get the covariate's value from input csv file in TableRecalibrationWalker
    public final Comparable getValue( final String str ) {
        return Integer.parseInt( str );
    }

}
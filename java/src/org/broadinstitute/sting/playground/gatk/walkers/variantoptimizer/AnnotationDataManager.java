package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;
import java.io.IOException;
import java.io.PrintStream;
import java.io.FileNotFoundException;

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

public class AnnotationDataManager {

    public final HashMap<String, TreeSet<AnnotationDatum>> data;

    public AnnotationDataManager() {
        data = new HashMap<String, TreeSet<AnnotationDatum>>();
    }

    public void addAnnotations( final RodVCF variant, final String sampleName ) {

        if( sampleName != null ) { // Only process variants that are found in the sample with this sampleName
            if( variant.getGenotype(sampleName).isNoCall() ) { // This variant isn't found in this sample so break out
                return;
            }
        } // else, process all samples

        // Loop over each annotation in the vcf record
        final Map<String,String> infoField = variant.getInfoValues();
        for( String annotationKey : infoField.keySet() ) {

            float value = 0.0f;
            try {
                value = Float.parseFloat( infoField.get( annotationKey ) );
            } catch( NumberFormatException e ) {
                continue; // Skip over annotations that aren't floats, like "DB"
            }

            TreeSet<AnnotationDatum> treeSet = data.get( annotationKey );
            if( treeSet == null ) { // This annotation hasn't been seen before
                treeSet = new TreeSet<AnnotationDatum>( new AnnotationDatum() ); // AnnotationDatum is a Comparator that orders variants by the value of the Annotation
                data.put( annotationKey, treeSet );
            }
            AnnotationDatum datum = new AnnotationDatum( value );
            if( treeSet.contains(datum) ) { // contains uses AnnotationDatum's equals function, so it only checks if the value field is already present
                datum = treeSet.tailSet(datum).first();
            } else {
                treeSet.add(datum);
            }

            final boolean isNovelVariant = variant.getID().equals(".");
            final boolean isTrueVariant = false; //BUGBUG: Check truth input file to see if this variant is in the truth set

            // Decide if the variant is a transition or transversion
            if( BaseUtils.isTransition( variant.getReferenceForSNP(), variant.getAlternativeBaseForSNP()) ) {
                datum.incrementTi( isNovelVariant, isTrueVariant );
            } else {
                datum.incrementTv( isNovelVariant, isTrueVariant );
            }
        }
    }

    public void plotCumulativeTables( final String PATH_TO_RSCRIPT, final String PATH_TO_RESOURCES, final String OUTPUT_PREFIX,
                                      final int MIN_VARIANTS_PER_BIN, final int MAX_VARIANTS_PER_BIN ) {

        final AnnotationDatum thisAnnotationBin = new AnnotationDatum();
        System.out.println( "\nExecuting RScript commands:" );

        // For each annotation we've seen
        for( String annotationKey : data.keySet() ) {

            PrintStream output;
            try {
                output = new PrintStream(OUTPUT_PREFIX + annotationKey + ".dat"); // Create the data file for this annotation
            } catch ( FileNotFoundException e ) {
                throw new StingException("Can't create intermediate output annotation data file. Does the output directory exist? " +
                                            OUTPUT_PREFIX + annotationKey + ".dat");
            }

            // Output a header line
            output.println("value\ttitv\tnumVariants\tcategory");

            // Bin SNPs and calculate truth metrics for each bin
            thisAnnotationBin.clearBin();
            for( AnnotationDatum datum : data.get( annotationKey ) ) {
                thisAnnotationBin.combine( datum );
                if( thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) >= MAX_VARIANTS_PER_BIN ) { // This annotation bin is full
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.FULL_SET ) + "\t" + thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) + "\tall");
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.NOVEL_SET ) + "\t" + thisAnnotationBin.numVariants( AnnotationDatum.NOVEL_SET ) + "\tnovel");
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.DBSNP_SET ) + "\t" + thisAnnotationBin.numVariants( AnnotationDatum.DBSNP_SET ) + "\tdbsnp");
                    thisAnnotationBin.clearBin();
                }
                // else, continue accumulating variants because this bin isn't full yet
            }

            if( thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) != 0 ) { // One final bin that may not have been dumped out
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.FULL_SET ) + "\t" + thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) + "\tall");
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.NOVEL_SET ) + "\t" + thisAnnotationBin.numVariants( AnnotationDatum.NOVEL_SET ) + "\tnovel");
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.DBSNP_SET ) + "\t" + thisAnnotationBin.numVariants( AnnotationDatum.DBSNP_SET ) + "\tdbsnp");

            }

            // Close the PrintStream
            output.close();

            // Print out the command line to make it clear to the user what is being executed and how one might modify it
            final String rScriptCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_BinnedTiTv.R" + " " +
                                   OUTPUT_PREFIX + annotationKey + ".dat" + " " + annotationKey + " " + MIN_VARIANTS_PER_BIN;
            System.out.println( rScriptCommandLine );

            // Execute the RScript command to plot the table of TiTv values
            try {
                Runtime.getRuntime().exec( rScriptCommandLine );
            } catch ( IOException e ) {
                throw new StingException( "Unable to execute RScript command: " + rScriptCommandLine );
            }
        }
    }
}

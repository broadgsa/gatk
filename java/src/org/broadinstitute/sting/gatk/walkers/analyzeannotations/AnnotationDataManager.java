package org.broadinstitute.sting.gatk.walkers.analyzeannotations;

import org.broad.tribble.util.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.contexts.variantcontext.VariantContextUtils;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
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

    private final HashMap<String, TreeSet<AnnotationDatum>> data;
    private final HashMap<String,String> nameMap;
    private final boolean INDICATE_MEAN_NUM_VARS;

    public AnnotationDataManager( HashMap<String,String> _nameMap, boolean _INDICATE_MEAN_NUM_VARS ) {
        data = new HashMap<String, TreeSet<AnnotationDatum>>();
        nameMap = _nameMap;
        INDICATE_MEAN_NUM_VARS = _INDICATE_MEAN_NUM_VARS;
    }

    public void addAnnotations( final VariantContext vc, final String sampleName, final boolean isInTruthSet, final boolean isTrueVariant ) {

        if( sampleName != null ) { // Only process variants that are found in the sample with this sampleName
            if( vc.getGenotype(sampleName).isNoCall() ) { // This variant isn't found in this sample so break out
                return;
            }
        } // else, process all samples

        // Loop over each annotation in the vcf record
        final Map<String,Object> infoField = new HashMap<String, Object>(vc.getAttributes());
        infoField.put("QUAL", ((Double)vc.getPhredScaledQual()).toString() ); // add QUAL field to annotations
        for( Map.Entry<String, Object> annotation : infoField.entrySet() ) {

            float value;
            try {
                value = Float.parseFloat( annotation.getValue().toString() );
            } catch( NumberFormatException e ) {
                continue; // Skip over annotations that aren't floats, like "DB"
            }

            TreeSet<AnnotationDatum> treeSet = data.get( annotation.getKey() );
            if( treeSet == null ) { // This annotation hasn't been seen before
                treeSet = new TreeSet<AnnotationDatum>( new AnnotationDatum() ); // AnnotationDatum is a Comparator that orders variants by the value of the Annotation
                data.put( annotation.getKey(), treeSet );
            }
            AnnotationDatum datum = new AnnotationDatum( value );
            if( treeSet.contains(datum) ) { // contains() uses AnnotationDatum's equals function, so it only checks if the value field is already present
                datum = treeSet.tailSet(datum).first();
            } else {
                treeSet.add(datum);
            }

            final boolean isNovelVariant = !vc.hasID() || !vc.getID().contains("rs");

            // Decide if the variant is a transition or transversion
            if( VariantContextUtils.getSNPSubstitutionType(vc).compareTo(BaseUtils.BaseSubstitutionType.TRANSITION) == 0 ) {
                datum.incrementTi( isNovelVariant, isInTruthSet, isTrueVariant );
            } else {
                datum.incrementTv( isNovelVariant, isInTruthSet, isTrueVariant );
            }
        }
    }

    public void plotCumulativeTables( final String PATH_TO_RSCRIPT, final String PATH_TO_RESOURCES, final String OUTPUT_PREFIX,
                                      final int MIN_VARIANTS_PER_BIN, final int MAX_VARIANTS_PER_BIN ) {

        final AnnotationDatum thisAnnotationBin = new AnnotationDatum();
        System.out.println( "\nFinished reading variants into memory. Executing RScript commands:" );

        // For each annotation we've seen
        for( final String annotationKey : data.keySet() ) {

            String filename = OUTPUT_PREFIX + annotationKey + ".dat";
            PrintStream output;
            try {
                output = new PrintStream(filename); // Create the intermediate data file for this annotation
            } catch ( FileNotFoundException e ) {
                throw new UserException.CouldNotCreateOutputFile(new File(filename), "Can't create intermediate output annotation data file. Does the output directory exist?", e);
            }

            // Output a header line
            output.println("value\ttitv\tdbsnp\ttruePositive\tnumVariants\tcategory");

            // Bin SNPs and calculate truth metrics for each bin
            thisAnnotationBin.clearBin();
            for( final AnnotationDatum datum : data.get( annotationKey ) ) {
                thisAnnotationBin.combine( datum );
                if( thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) >= MAX_VARIANTS_PER_BIN ) { // This annotation bin is full
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.FULL_SET ) + "\t" + thisAnnotationBin.calcDBsnpRate() + "\t" + thisAnnotationBin.calcTPrate() +
                            "\t" + thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) + "\tall");
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.NOVEL_SET ) + "\t-1\t-1\t" + thisAnnotationBin.numVariants( AnnotationDatum.NOVEL_SET ) + "\tnovel");
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.DBSNP_SET ) + "\t-1\t-1\t" + thisAnnotationBin.numVariants( AnnotationDatum.DBSNP_SET ) + "\tdbsnp");
                    output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.TRUTH_SET ) + "\t-1\t-1\t" + thisAnnotationBin.numVariants( AnnotationDatum.TRUTH_SET ) + "\ttruth");
                    thisAnnotationBin.clearBin();
                }
                // else, continue accumulating variants because this bin isn't full yet
            }

            // One final bin that may not have been dumped out
            if( thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) != 0 ) {
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.FULL_SET ) + "\t" + thisAnnotationBin.calcDBsnpRate() + "\t" + thisAnnotationBin.calcTPrate() +
                            "\t" + thisAnnotationBin.numVariants( AnnotationDatum.FULL_SET ) + "\tall");
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.NOVEL_SET ) + "\t-1\t-1\t" + thisAnnotationBin.numVariants( AnnotationDatum.NOVEL_SET ) + "\tnovel");
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.DBSNP_SET ) + "\t-1\t-1\t" + thisAnnotationBin.numVariants( AnnotationDatum.DBSNP_SET ) + "\tdbsnp");
                output.println( thisAnnotationBin.value + "\t" + thisAnnotationBin.calcTiTv( AnnotationDatum.TRUTH_SET ) + "\t-1\t-1\t" + thisAnnotationBin.numVariants( AnnotationDatum.TRUTH_SET ) + "\ttruth");
                thisAnnotationBin.clearBin();
            }

            // Close the PrintStream
            output.close();

            String annotationName = null;
            if( nameMap != null ) {
                annotationName = nameMap.get(annotationKey);
            }
            if( annotationName == null ) { // This annotation is not in the name map so use the key instead
                annotationName = annotationKey;
            }

            // Print out the command line to make it clear to the user what is being executed and how one might modify it
            final String rScriptCommandLine = PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_BinnedTruthMetrics.R" + " " +
                                   OUTPUT_PREFIX + annotationKey + ".dat" + " " + annotationName + " " + MIN_VARIANTS_PER_BIN + " " + INDICATE_MEAN_NUM_VARS;
            System.out.println( rScriptCommandLine );

            // Execute the RScript command to plot the table of truth values
            try {
                Runtime.getRuntime().exec( rScriptCommandLine );
            } catch ( IOException e ) {
                throw new UserException.CannotExecuteRScript( rScriptCommandLine, e );
            }
        }
    }
}

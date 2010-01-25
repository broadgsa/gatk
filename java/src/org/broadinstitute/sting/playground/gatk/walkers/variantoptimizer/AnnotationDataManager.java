package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;
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

    public final HashMap<String, TreeSet<AnnotationDatum>> dataFull;
    public final HashMap<String, TreeSet<AnnotationDatum>> dataNovel;
    public final HashMap<String, TreeSet<AnnotationDatum>> dataDBsnp;
    public final HashMap<String, TreeSet<AnnotationDatum>> dataTruthSet;

    public AnnotationDataManager() {
        dataFull = new HashMap<String, TreeSet<AnnotationDatum>>();
        dataNovel = new HashMap<String, TreeSet<AnnotationDatum>>();
        dataDBsnp = new HashMap<String, TreeSet<AnnotationDatum>>();
        dataTruthSet = new HashMap<String, TreeSet<AnnotationDatum>>();
    }

    public void addAnnotations( RodVCF variant, String sampleName ) {

        if( sampleName != null ) { // only process variants that are found in the sample with this sampleName
            if( variant.getGenotype(sampleName).isNoCall() ) { // this variant isn't found in this sample so break out
                return;
            }
        } // else, process all samples

        // Loop over each annotation in the vcf record
        final Map<String,String> infoField = variant.getInfoValues();
        for( String annotationKey : infoField.keySet() ) {

            float value = 0.0f;
            boolean skipThisAnnotation = false;
            try {
                value = Float.parseFloat( infoField.get( annotationKey ) );
            } catch( Exception e ) { // BUGBUG: Make this exception more specific. NumberFormatException??
                skipThisAnnotation = true; // skip over annotations that aren't floats, like "AC"?? and "DB"
            }
            if( skipThisAnnotation ) {
                continue; // skip over annotations that aren't floats, like "AC"?? and "DB"
            }

            if( variant.getID().equals(".") ) { // This is a novel variant
                TreeSet<AnnotationDatum> treeSet = dataNovel.get( annotationKey );
                if( treeSet == null ) { // This annotation hasn't been seen before
                    treeSet = new TreeSet<AnnotationDatum>( new AnnotationDatum() ); // AnnotationDatum is a Comparator that orders variants by the value of the Annotation
                    dataNovel.put( annotationKey, treeSet );
                }
                AnnotationDatum datum = new AnnotationDatum( value );
                if( treeSet.contains(datum) ) {
                    datum = treeSet.tailSet(datum).first();
                } else {
                    treeSet.add(datum);
                }

                // Decide if it was a Ti or a Tv
                if( BaseUtils.isTransition(variant.getReferenceForSNP(), variant.getAlternativeBaseForSNP()) ) {
                    datum.incrementTi();
                } else {
                    datum.incrementTv();
                }
            } else { // This is a DBsnp variant
                TreeSet<AnnotationDatum> treeSet = dataDBsnp.get( annotationKey );
                if( treeSet == null ) { // This annotation hasn't been seen before
                    treeSet = new TreeSet<AnnotationDatum>( new AnnotationDatum() ); // AnnotationDatum is a Comparator that orders variants by the value of the Annotation
                    dataDBsnp.put( annotationKey, treeSet );
                }
                AnnotationDatum datum = new AnnotationDatum( value );
                if( treeSet.contains(datum) ) {
                    datum = treeSet.tailSet(datum).first();
                } else {
                    treeSet.add(datum);
                }

                // Decide if it was a Ti or a Tv
                if( BaseUtils.isTransition(variant.getReferenceForSNP(), variant.getAlternativeBaseForSNP()) ) {
                    datum.incrementTi();
                } else {
                    datum.incrementTv();
                }
            }

            // the overall set, containing both novel and DBsnp variants
            TreeSet<AnnotationDatum> treeSet = dataFull.get( annotationKey );
            if( treeSet == null ) { // This annotation hasn't been seen before
                treeSet = new TreeSet<AnnotationDatum>( new AnnotationDatum() ); // AnnotationDatum is a Comparator that orders variants by the value of the Annotation
                dataFull.put( annotationKey, treeSet );
            }
            AnnotationDatum datum = new AnnotationDatum( value );
            if( treeSet.contains(datum) ) {
                datum = treeSet.tailSet(datum).first();
            } else {
                treeSet.add(datum);
            }

            // Decide if it was a Ti or a Tv
            if( BaseUtils.isTransition(variant.getReferenceForSNP(), variant.getAlternativeBaseForSNP()) ) {
                datum.incrementTi();
            } else {
                datum.incrementTv();
            }
        }
    }

    public void plotCumulativeTables( final String PATH_TO_RSCRIPT, final String PATH_TO_RESOURCES, final String OUTPUT_DIR,
                                      final int MIN_VARIANTS_PER_BIN, final int MAX_VARIANTS_PER_BIN ) {

        for( String annotationKey: dataFull.keySet() ) {

            /*
            PrintStream output;
            try {
                output = new PrintStream(OUTPUT_DIR + annotationKey + ".cumulative.dat");
            } catch (FileNotFoundException e) {
                throw new StingException("Can't create intermediate output annotation data file.");
            }
            // Output a header line
            output.println("value\tcumulativeTiTv\tnumVariants\tGT");

            // Filter SNPs greater than this annotation value
            int numTi = 0;
            int numTv = 0;
            for( AnnotationDatum datum : data.get( annotationKey ) ) {
                numTi += datum.numTransitions;
                numTv += datum.numTransversions;
                float titv;
                if( numTv == 0) { titv = 0.0f; }
                else { titv = ((float) numTi) / ((float) numTv); }
                output.println(datum.value + "\t" + titv + "\t" + (numTi+numTv) +"\t1");

            }

            // Filter SNPs less than this annotation value
            numTi = 0;
            numTv = 0;
            Iterator<AnnotationDatum> iter = data.get( annotationKey ).descendingIterator();
            while( iter.hasNext() ) {
                final AnnotationDatum datum = iter.next();
                numTi += datum.numTransitions;
                numTv += datum.numTransversions;
                float titv;
                if( numTv == 0) { titv = 0.0f; }
                else { titv = ((float) numTi) / ((float) numTv); }
                output.println(datum.value + "\t" + titv + "\t" + (numTi+numTv) +"\t0");

            }
            */

            PrintStream output;
            try {
                output = new PrintStream(OUTPUT_DIR + annotationKey + ".binned.dat");
            } catch (FileNotFoundException e) {
                throw new StingException("Can't create intermediate output annotation data file.");
            }

            // Output a header line
            output.println("value\ttitv\tnumVariants\tcategory");

            // Bin SNPs and calculate truth metrics for each bin, right now this is just TiTv
            int numTi = 0;
            int numTv = 0;
            AnnotationDatum lastDatum = null;
            for( AnnotationDatum datum : dataFull.get( annotationKey ) ) {
                numTi += datum.numTransitions;
                numTv += datum.numTransversions;
                lastDatum = datum;
                if( numTi + numTv >= MAX_VARIANTS_PER_BIN ) { // This annotation bin is full
                    float titv;
                    if( numTv == 0) { titv = 0.0f; }
                    else { titv = ((float) numTi) / ((float) numTv); }
                    output.println(datum.value + "\t" + titv + "\t" + (numTi+numTv) + "\t0");
                    numTi = 0;
                    numTv = 0;
                }
                // else, continue accumulating variants because this bin isn't full yet
            }
            if( numTi != 0 || numTv != 0 ) { // one final bin that may not have been dumped out
                float titv;
                if( numTv == 0) { titv = 0.0f; }
                else { titv = ((float) numTi) / ((float) numTv); }
                output.println(lastDatum.value + "\t" + titv + "\t" + (numTi+numTv)+ "\t0");
            }


            numTi = 0;
            numTv = 0;
            lastDatum = null;
            for( AnnotationDatum datum : dataNovel.get( annotationKey ) ) {
                numTi += datum.numTransitions;
                numTv += datum.numTransversions;
                lastDatum = datum;
                if( numTi + numTv >= MAX_VARIANTS_PER_BIN/2 ) { // This annotation bin is full
                    float titv;
                    if( numTv == 0) { titv = 0.0f; }
                    else { titv = ((float) numTi) / ((float) numTv); }
                    output.println(datum.value + "\t" + titv + "\t" + (numTi+numTv) + "\t1");
                    numTi = 0;
                    numTv = 0;
                }
                // else, continue accumulating variants because this bin isn't full yet
            }
            if( numTi != 0 || numTv != 0 ) { // one final bin that may not have been dumped out
                float titv;
                if( numTv == 0) { titv = 0.0f; }
                else { titv = ((float) numTi) / ((float) numTv); }
                output.println(lastDatum.value + "\t" + titv + "\t" + (numTi+numTv)+ "\t1");
            }

            numTi = 0;
            numTv = 0;
            lastDatum = null;
            for( AnnotationDatum datum : dataDBsnp.get( annotationKey ) ) {
                numTi += datum.numTransitions;
                numTv += datum.numTransversions;
                lastDatum = datum;
                if( numTi + numTv >= MAX_VARIANTS_PER_BIN ) { // This annotation bin is full
                    float titv;
                    if( numTv == 0) { titv = 0.0f; }
                    else { titv = ((float) numTi) / ((float) numTv); }
                    output.println(datum.value + "\t" + titv + "\t" + (numTi+numTv) + "\t2");
                    numTi = 0;
                    numTv = 0;
                }
                // else, continue accumulating variants because this bin isn't full yet
            }
            if( numTi != 0 || numTv != 0 ) { // one final bin that may not have been dumped out
                float titv;
                if( numTv == 0) { titv = 0.0f; }
                else { titv = ((float) numTi) / ((float) numTv); }
                output.println(lastDatum.value + "\t" + titv + "\t" + (numTi+numTv)+ "\t2");
            }


            /*
            System.out.println(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_CumulativeTiTv.R" + " " +
                                        OUTPUT_DIR + annotationKey + ".cumulative.dat" + " " + OUTPUT_DIR + " " + annotationKey);

            try {
                Process p = Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_CumulativeTiTv.R" + " " +
                                        OUTPUT_DIR + annotationKey + ".cumulative.dat"+  " " + OUTPUT_DIR + " " + annotationKey);
            } catch (Exception e) {
                throw new StingException("Unable to execute RScript command");
            }
            */

            System.out.println(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_BinnedTiTv.R" + " " +
                                        OUTPUT_DIR + annotationKey + ".binned.dat" + " " + OUTPUT_DIR + " " + annotationKey +
                                        " " + MIN_VARIANTS_PER_BIN);

            // Execute the RScript command to plot the table of TiTv values
            try {
                final Process p = Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_BinnedTiTv.R" + " " +
                                        OUTPUT_DIR + annotationKey + ".binned.dat" + " " + OUTPUT_DIR + " " + annotationKey +
                                        " " + MIN_VARIANTS_PER_BIN);
            } catch (Exception e) {
                throw new StingException("Unable to execute RScript command");
            }
        }
    }
}

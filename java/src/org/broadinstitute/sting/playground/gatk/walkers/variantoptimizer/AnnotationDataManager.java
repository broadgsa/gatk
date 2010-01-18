package org.broadinstitute.sting.playground.gatk.walkers.variantoptimizer;

import org.broadinstitute.sting.gatk.refdata.RodVCF;
import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.StingException;

import java.util.HashMap;
import java.util.TreeSet;
import java.util.Map;
import java.util.Iterator;
import java.io.PrintStream;
import java.io.FileNotFoundException;

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

    public void addAnnotations( RodVCF variant ) {

        // Loop over each annotation in the vcf record
        final Map<String,String> infoField = variant.getInfoValues();
        for( String annotationKey : infoField.keySet() ) {
            final float value = Float.parseFloat( infoField.get( annotationKey ) );

            TreeSet<AnnotationDatum> treeSet = data.get( annotationKey );
            if( treeSet == null ) { // This annotation hasn't been seen before
                treeSet = new TreeSet<AnnotationDatum>( new AnnotationDatum() );
                data.put( annotationKey, treeSet );
            }
            AnnotationDatum datum = new AnnotationDatum( value );
            if( treeSet.contains(datum) ) {
                datum = treeSet.tailSet(datum).first();
            } else {
                treeSet.add(datum);
            }

            if( BaseUtils.isTransition(variant.getReferenceForSNP(), variant.getAlternativeBaseForSNP()) ) {
                datum.incrementTi();
            } else {
                datum.incrementTv();
            }
        }
    }

    public void plotCumulativeTables(final String PATH_TO_RSCRIPT, final String PATH_TO_RESOURCES, final String OUTPUT_DIR) {

        for( String annotationKey: data.keySet() ) {
            PrintStream output;
            try {
                output = new PrintStream(OUTPUT_DIR + annotationKey + ".dat");
            } catch (FileNotFoundException e) {
                throw new StingException("Can't create intermediate output annotation data file.");
            }
            // Output a header line
            output.println("value\tcumulativeTiTv\tnumVariants\tGT");

            int numTi = 0;
            int numTv = 0;
            for( AnnotationDatum datum : data.get( annotationKey )) {
                numTi += datum.numTransitions;
                numTv += datum.numTransversions;
                float titv;
                if( numTv == 0) { titv = 0.0f; }
                else { titv = ((float) numTi) / ((float) numTv); }
                output.println(datum.value + "\t" + titv + "\t" + (numTi+numTv) +"\t1");

            }

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

            System.out.println(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_CumulativeTiTv.R" + " " +
                                        OUTPUT_DIR + annotationKey + ".dat" + " " + annotationKey);

            try {
                Process p = Runtime.getRuntime().exec(PATH_TO_RSCRIPT + " " + PATH_TO_RESOURCES + "plot_Annotations_CumulativeTiTv.R" + " " +
                                        OUTPUT_DIR + annotationKey + ".dat" + " " + annotationKey);
            } catch (Exception e) {
                throw new StingException("Unable to run RScript command");
            }
        }
    }
}

package org.broadinstitute.sting.playground.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.Tags;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/12/11
 */

public class TrainingSet {

    public String name;
    public boolean isKnown = false;
    public boolean isTraining = false;
    public boolean isTruth = false;
    public double prior = 3.0;

    protected final static Logger logger = Logger.getLogger(TrainingSet.class);

    public TrainingSet( final String name, final Tags tags ) {
        this.name = name;
        if( tags != null ) {
            isKnown = tags.containsKey("known") && tags.getValue("known").equals("true");
            isTraining = tags.containsKey("training") && tags.getValue("training").equals("true");
            isTruth = tags.containsKey("truth") && tags.getValue("truth").equals("true");
            prior = ( tags.containsKey("known") ? Double.parseDouble(tags.getValue("prior")) : prior );
        }
        logger.info(String.format( "Found %s track: \tKnown = %s \tTraining = %s \tTruth = %s \tPrior = Q%.1f", this.name, isKnown, isTraining, isTruth, prior) );
    }
}

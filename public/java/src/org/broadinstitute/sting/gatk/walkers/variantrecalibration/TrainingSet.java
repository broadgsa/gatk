package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

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
    public boolean isAntiTraining = false;
    public boolean isTruth = false;
    public boolean isConsensus = false;
    public double prior = 0.0;

    protected final static Logger logger = Logger.getLogger(TrainingSet.class);

    public TrainingSet( final String name, final Tags tags ) {
        this.name = name;

        // Parse the tags to decide which tracks have which properties
        if( tags != null ) {
            isKnown = tags.containsKey("known") && tags.getValue("known").equals("true");
            isTraining = tags.containsKey("training") && tags.getValue("training").equals("true");
            isAntiTraining = tags.containsKey("bad") && tags.getValue("bad").equals("true");
            isTruth = tags.containsKey("truth") && tags.getValue("truth").equals("true");
            isConsensus = tags.containsKey("consensus") && tags.getValue("consensus").equals("true");
            prior = ( tags.containsKey("prior") ? Double.parseDouble(tags.getValue("prior")) : prior );
        }

        // Report back to the user which tracks were found and the properties that were detected
        if( !isConsensus && !isAntiTraining ) {
            logger.info( String.format( "Found %s track: \tKnown = %s \tTraining = %s \tTruth = %s \tPrior = Q%.1f", this.name, isKnown, isTraining, isTruth, prior) );
        } else if( isConsensus ) {
            logger.info( String.format( "Found consensus track: %s", this.name) );
        } else {
            logger.info( String.format( "Found bad sites training track: %s", this.name) );
        }
    }
}

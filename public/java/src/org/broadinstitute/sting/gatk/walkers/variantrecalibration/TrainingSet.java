/*
 * Copyright (c) 2011 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.gatk.walkers.variantrecalibration;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.RodBinding;
import org.broadinstitute.sting.commandline.Tags;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 3/12/11
 */

public class TrainingSet {

    public RodBinding<VariantContext> rodBinding;
    public boolean isKnown = false;
    public boolean isTraining = false;
    public boolean isAntiTraining = false;
    public boolean isTruth = false;
    public boolean isConsensus = false;
    public double prior = 0.0;

    protected final static Logger logger = Logger.getLogger(TrainingSet.class);

    public TrainingSet( final RodBinding<VariantContext> rodBinding) {
        this.rodBinding = rodBinding;

        final Tags tags = rodBinding.getTags();
        final String name = rodBinding.getName();

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
            logger.info( String.format( "Found %s track: \tKnown = %s \tTraining = %s \tTruth = %s \tPrior = Q%.1f", name, isKnown, isTraining, isTruth, prior) );
        } else if( isConsensus ) {
            logger.info( String.format( "Found consensus track: %s", name) );
        } else {
            logger.info( String.format( "Found bad sites training track: %s", name) );
        }
    }
}

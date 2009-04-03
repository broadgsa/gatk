package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.gatk.refdata.AllelicVariant;

import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: andrewk
 * Date: Apr 2, 2009
 * Time: 9:01:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class AlleleMetricsWalker {
    // Class that will walk over various metrics in a reference ordered way
    // This class walks over the genome in reference order and calls AlleleMetrics on each class
    // Hapmap and dbSNP tracks are taken from the command line
    // At first pass, this will at least be able to walk over a GFF file and compare to the hapmap and dbsnp
    // tracks specified on the command line and handed in via the LocusContext

    public void map(List<AllelicVariant> avdata) {

    }
    

}

package org.broadinstitute.sting.gatk.dataSources.providers;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 8, 2009
 * Time: 5:01:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class ReferenceProvider {
    private ReferenceIterator reference;

    public ReferenceProvider( ReferenceIterator reference ) {
        this.reference = reference;
    }

    public ReferenceIterator getReferenceSequence( GenomeLoc genomeLoc ) {
        if( (genomeLoc.getStop() - genomeLoc.getStart()) > 0 )
            throw new RuntimeException( "Internal error :LocusContextProviders currently require 1-base genomeLocs.");

        // jump to the first reference site
        return reference.seekForward(genomeLoc);                    
    }
}

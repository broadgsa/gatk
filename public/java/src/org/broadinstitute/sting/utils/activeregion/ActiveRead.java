package org.broadinstitute.sting.utils.activeregion;

import org.broadinstitute.sting.utils.sam.GATKSAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: 1/4/12
 */

public class ActiveRead {
    final public GATKSAMRecord read;
    final public boolean isPrimaryRegion;

    ActiveRead( final GATKSAMRecord read, final boolean isPrimaryRegion ) {
        this.read = read;
        this.isPrimaryRegion = isPrimaryRegion;
    }
}

package org.broadinstitute.sting.gatk.datasources.providers;

import org.broadinstitute.sting.utils.fasta.IndexedFastaSequenceFile;

import java.util.Collections;
import java.util.Collection;
/**
 * User: hanna
 * Date: May 22, 2009
 * Time: 12:19:17 PM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * A view into the reference backing this shard.
 */
public class ReferenceView implements View {
    /**
     * The source of reference data.
     */
    protected IndexedFastaSequenceFile reference = null;

    /**
     * Create a new ReferenceView.
     * @param provider
     */
    public ReferenceView( ShardDataProvider provider ) {
        this.reference = provider.getReference();
    }

    /**
     * Reference views don't conflict with anything else.
     * @return Empty list.
     */
    public Collection<Class<? extends View>> getConflictingViews() { return Collections.emptyList(); }

    /**
     * Deinitialize pointers for fast fail.  Someone else will handle file management.
     */
    public void close() {
        reference = null;
    }
}

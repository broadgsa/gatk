package org.broadinstitute.sting.gatk.iterators;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.Iterator;
/**
 * User: hanna
 * Date: May 12, 2009
 * Time: 10:48:56 AM
 * BROAD INSTITUTE SOFTWARE COPYRIGHT NOTICE AND AGREEMENT
 * Software and documentation are copyright 2005 by the Broad Institute.
 * All rights are reserved.
 *
 * Users acknowledge that this software is supplied without any warranty or support.
 * The Broad Institute is not responsible for its use, misuse, or
 * functionality.
 */

/**
 * An iterator for iterating through loci on a genome.
 */

public interface LocusIterator extends Iterator<GenomeLoc> {
}

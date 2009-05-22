package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;

/** A simple wrapper interface that allows for keeping any of the "main" (point mutation) or "associated" (indel)
 * genotypes, or both, at any given reference position. This is the abstraction of "complete" genotyping information:
 * one might have or might have not genotyped independently a base at the given position (point) or an indel right after
 * that position. The query and getter methods allow clients to figure out what pieces of the genotyping information
 * are actually available and to retrieve them. Behavior of the getters when corresponding queries return false is
 * not bound by a contract but it is suggested that they return null or a "trivial" (reference) genotype. 
 * 
 * @author asivache
 *
 */

public interface GenotypeList {
	/** Returns true if point mutation genotyping information is available at the
	 * current position.
	 * @return
	 */
	boolean hasPointGenotype();

	/** Returns true if indel genotyping information is available at the
	 * current position.
	 * @return
	 */
	boolean hasIndelGenotype();
	
	/** Returns point mutation genotyping information if available;
	 * behavior is not specified in the case when there is no such information.
	 * @return
	 */
	Genotype getPointGenotype();

	/** Returns indel genotyping information if available;
	 * behavior is not specified in the case when there is no such information.
	 * @return
	 */
	Genotype getIndelGenotype();
	
	/** Returns position on the reference the stored genotyping data are associated with. 
	 * 
	 * @return
	 */
	GenomeLoc getLocation();
}

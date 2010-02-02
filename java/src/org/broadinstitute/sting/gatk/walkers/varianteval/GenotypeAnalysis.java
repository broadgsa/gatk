package org.broadinstitute.sting.gatk.walkers.varianteval;

/**
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 * If an analysis implements this interface, it asserts that it performs a genotype based analysis, as
 * opposed a straight variant analysis.  The difference here is that variants are not asserted to be
 * the actual genotype of a particular person, but are really just variation "out-there" in a population.
 * A genotype analysis would be something like covered bases, confidently called bases, genotyping
 * concordance, etc.
 *
 */
public interface GenotypeAnalysis {
    
}
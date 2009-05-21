package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.Utils;

/**
 * This class is a trivial wrapper for keeping together and passing around a few simple statistics relevant for genotyping:
 * e.g. number of covered bases (have any observation at all), number of "assessed" bases (e.g. those with confidence level
 * above some cutoff, so that a call was actually made), number of ref/variant calls etc.
 * @author asivache
 *
 */
public class GenotypingCallStats {
	public long covered = 0; // number of loci covered in all 3 individuals (not necessarily confidently called)
	public long assessed = 0; // number of loci with all 3 genotypes available at or above the specified cutoff
	public long ref = 0; // number of assessed loci, where all 3 people have homogeneous reference allele
	public int variant = 0; // number of assessed loci where a variant is observed in the individual
	public int consistent_variant = 0; // variants that are consistent in any (application-specific) sense, e.g. variant matches variants in other members of the family trio
	public int consistent_ref = 0; // reference calls that are consistent in any (app-specific) sense, e.g. consistent with other members of the family trio
	public int non_biallelic_variant = 0; // number of variant calls that are not biallelic
	
	public GenotypingCallStats add(GenotypingCallStats other) {
		this.covered += other.covered;
		this.assessed += other.assessed;
		this.ref += other.ref;
		this.variant += other.variant;
		this.consistent_variant += other.consistent_variant;
		this.consistent_ref += other.consistent_ref;
		this.non_biallelic_variant += other.non_biallelic_variant;
		return this;
	}
	
//	public int totalVariants() { return consistent_variant + inconsistent_variant + non_biallelic_variant; }
	
	public String toString() {
		StringBuilder b = new StringBuilder();
		
		
		b.append( String.format("      covered: %d%n      assessed: %d (%3.2f%%)%n         ref: %d (%3.2f%%)%n         variants: %d (%3.2f%%)%n            multiallelic: %d (%3.2f%%)%n", 
				covered, assessed, Utils.percentage(assessed, covered), 
				ref, Utils.percentage(ref,assessed),
				variant, Utils.percentage(variant,assessed),
				non_biallelic_variant, Utils.percentage(non_biallelic_variant, variant))
			);
		
		return b.toString();
	}
}

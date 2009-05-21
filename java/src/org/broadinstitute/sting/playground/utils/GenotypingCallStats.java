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
	public long covered = 0; // number of loci covered in an individual (not necessarily confidently called)
	public long assessed = 0; // number of loci with confident calls
	public long ref = 0; // number of assessed loci where the reference is called
	public int variant = 0; // number of assessed loci where a variant is observed in the individual
	// NOTE: consistent_ref + inconsistent_ref is the total number of ref calls assessed for consistency (by some external application).
	// this  number does not have to be equal to 'ref' ( total number of ref calls in this individual - we migh be unable to assess
	// the consistency for all of them!); same applies to (in)consistent_variant.
	public int consistent_variant = 0; // variants that are consistent in any (application-specific) sense, e.g. variant matches variants in other members of the family trio
	public int consistent_ref = 0; // reference calls that are consistent in any (app-specific) sense, e.g. consistent with other members of the family trio
	public int inconsistent_variant = 0; // variants that are inconsistent in any (application-specific) sense, e.g. variant does not match variants in other members of the family trio
	public int inconsistent_ref = 0; // reference calls that are inconsistent in any (app-specific) sense, e.g. inconsistent with other members of the family trio
	public int non_biallelic_variant = 0; // number of variant calls that are not biallelic
	
	public GenotypingCallStats add(GenotypingCallStats other) {
		this.covered += other.covered;
		this.assessed += other.assessed;
		this.ref += other.ref;
		this.variant += other.variant;
		this.consistent_variant += other.consistent_variant;
		this.consistent_ref += other.consistent_ref;
		this.inconsistent_variant += other.consistent_variant;
		this.inconsistent_ref += other.consistent_ref;
		this.non_biallelic_variant += other.non_biallelic_variant;
		return this;
	}
	
//	public int totalVariants() { return consistent_variant + inconsistent_variant + non_biallelic_variant; }
	
	public String toString() {
		StringBuilder b = new StringBuilder();
		
		
		b.append( String.format("      covered: %d%n      assessed: %d (%3.2f%% covered)%n", 
				covered, assessed, Utils.percentage(assessed, covered) )
			);
		
		b.append( String.format("         ref: %d (%3.2f%% assessed)%n", 
				ref, Utils.percentage(ref,assessed))
			);

		int z = consistent_ref+inconsistent_ref;
		b.append( String.format("            ref assessed for consistency: %d (%3.2f%% ref)%n", 
				z, Utils.percentage(z, ref) )
			);
		b.append( String.format("               consistent ref: %d (%3.2f%% consistency-assessed ref)%n", 
				consistent_ref, Utils.percentage(consistent_ref,z) )
			);

		b.append( String.format("         variants: %d (%3.2f%% assessed)%n", 
				variant, Utils.percentage(variant,assessed))
			);
		b.append( String.format("            multiallelic: %d (%3.2f%% variants)%n", 
				non_biallelic_variant, Utils.percentage(non_biallelic_variant, variant))
			);
		
		z = consistent_variant+inconsistent_variant;
		b.append( String.format("            variants assessed for consistency: %d (%3.2f%% variants)%n", 
				z, Utils.percentage(z, variant) )
			);
		b.append( String.format("               consistent variants: %d (%3.2f%% consistency-assessed variants)%n", 
				consistent_ref, Utils.percentage(consistent_variant,z) )
			);
		return b.toString();
	}
}

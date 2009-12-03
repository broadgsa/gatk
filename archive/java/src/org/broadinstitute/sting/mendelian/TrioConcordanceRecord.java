package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.Utils;

/**
 * This class is a trivial wrapper for keeping together and passing around counts of different possible outcomes of 
 * comparisons in a trio.  Mendelian walker uses this class to classify/accumulate events such as consistent snp, inconsistent snp 
 * (e.g. only in a kid, but not in parents), loci with no calls etc. 
 * @author asivache
 *
 */
public class TrioConcordanceRecord {
	
	public GenotypingCallStats mom;
	public GenotypingCallStats dad;
	public GenotypingCallStats kid;
	
	public GenotypingCallStats trio;

//	public long mom_assessed_ref; // number of ref calls in mother on positions assessed *in all 3 individuals*
//	public long dad_assessed_ref; // ditto
//	public long kid_assessed_ref; 
//	public int mom_assessed_variant; // number of variant calls in mother on  positions assessed *in all 3 individuals*
//	public int dad_assessed_variant; // ditto
//	public int kid_assessed_variant; 
	public int missing_variant_in_kid;
	public int nonmatching_variant_in_kid;
	public int missing_variant_in_parents;
	public int mom_passed_variant;
	public int dad_passed_variant;

	//	public long consistent_ref = 0; // number of assessed loci, where all 3 people have homogeneous reference allele
//	public int consistent_variant = 0; // number of assessed loci where a variant is observed in at least one individual and genotyping calls are consistent between the trio members
//	public int inconsistent_variant = 0; // number of assessed loci where a variant is observed in at least one individual and genotyping calls are inconsistent
//	public int missing_variant_in_parents = 0; // number of inconsistent variants (see above), where parent(s) have a variant but the kid does not while she should
//	public int missing_variant_in_kid = 0; // number of inconsistent variants (see above), where kid has a snp but the parents do not while they should
//	public int consistent_variant_passed = 0; // variants that are consistent and *passed* (i.e. present in kid and one of the parents)
//	public int non_biallelic_variant = 0; // number of variant calls that are not biallelic
//	public long unclassified_events = 0;
	
	public TrioConcordanceRecord() {
		mom = new GenotypingCallStats();
		dad = new GenotypingCallStats();
		kid = new GenotypingCallStats();
		trio = new GenotypingCallStats();
	}
	
	public TrioConcordanceRecord add(TrioConcordanceRecord other) {
		
		this.mom.add(other.mom);
		this.dad.add(other.dad);
		this.kid.add(other.kid);

		this.trio.add(other.trio);

//		this.mom_assessed_ref += other.mom_assessed_ref;
//		this.dad_assessed_ref += other.dad_assessed_ref;
//		this.kid_assessed_ref += other.kid_assessed_ref;
//		this.mom_assessed_variant += other.mom_assessed_variant;
//		this.dad_assessed_variant += other.dad_assessed_ref;
//		this.kid_assessed_variant += other.kid_assessed_variant;
		this.missing_variant_in_kid += other.missing_variant_in_kid ;
		this.nonmatching_variant_in_kid += other.nonmatching_variant_in_kid ;
		this.missing_variant_in_parents += other.missing_variant_in_parents ;
		this.mom_passed_variant += other.mom_passed_variant;
		this.dad_passed_variant += other.dad_passed_variant;

		//		this.consistent_ref += other.consistent_ref;
//		this.consistent_variant += other.consistent_variant;
//		this.inconsistent_variant += other.inconsistent_variant;
//		this.missing_variant_in_parents += other.missing_variant_in_parents;
//		this.missing_variant_in_kid += other.missing_variant_in_kid;
//		this.consistent_variant_passed += other.consistent_variant_passed;
//		this.non_biallelic_variant += other.non_biallelic_variant;
//		this.unclassified_events += other.unclassified_events;
		return this;
	}
	
	public int totalVariants() { return trio.consistent_variant + trio.inconsistent_variant + trio.non_biallelic_variant; }
	
	public String toString() {
		StringBuilder b = new StringBuilder();
		
		b.append(String.format("%ncovered in trio: %d%n", trio.covered ) );

		b.append(String.format("assessed in trio: %d (%3.2f%% covered)%n", 
				trio.assessed, Utils.percentage(trio.assessed,trio.covered )) );
		
		b.append(String.format("   reference in all samples: %d (%3.2f%% assessed)%n", 
				trio.ref, Utils.percentage(trio.ref,trio.assessed )) );
		
		b.append(String.format("   variant sites: %d (%3.2f%% assessed, or 1 per %3.2f kB)%n", 
				totalVariants(), Utils.percentage(totalVariants(), trio.assessed), ((double)trio.assessed/totalVariants())/1000.0 
		));
		
		b.append(String.format("      consistent variants: %d (%3.2f%% variants)%n", 
				trio.consistent_variant, Utils.percentage(trio.consistent_variant,totalVariants()) 
		));

//		b.append(String.format("         passed (in daughter and parent(s)): %d%n         lost (in parent(s) but not in daughter): %d%n",
//    		   consistent_variant_passed, consistent_variant - consistent_variant_passed));
       
       b.append(String.format("      multiallelic variant: %d (%3.2f%% variants)%n",
				trio.non_biallelic_variant, Utils.percentage(trio.non_biallelic_variant, totalVariants())
       ));

       b.append(String.format("      inconsistent variant: %d (%3.2f%% variants)%n",
				trio.inconsistent_variant, Utils.percentage(trio.inconsistent_variant, totalVariants())
      ));

       b.append(String.format("         missing from daughter: %d (%3.2f%% inconsistent variants)%n",
				missing_variant_in_kid, Utils.percentage(missing_variant_in_kid, trio.inconsistent_variant)
      ));

       b.append(String.format("         missing from both parents: %d (%3.2f%% inconsistent variants)%n",
				missing_variant_in_parents, Utils.percentage(missing_variant_in_parents, trio.inconsistent_variant)		
       ));

       b.append(String.format("         non-matching in daughter: %d (%3.2f%% inconsistent variants)%n",
				nonmatching_variant_in_kid, Utils.percentage(nonmatching_variant_in_kid, trio.inconsistent_variant)		
       ));
       
		b.append("per trio individual:\n");
		b.append("   mother:\n");
		b.append(mom.toString());
		b.append("   father:\n");
		b.append(dad.toString());
		b.append("   daughter:\n");
		b.append(kid.toString());
		
		return b.toString();
	}
}

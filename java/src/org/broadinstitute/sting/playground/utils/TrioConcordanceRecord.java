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
	public long assessed_loci = 0; // number of loci with all 3 genotypes available at or above the specified cutoff
	public long consistent_ref = 0; // number of assessed loci, where all 3 people have homogeneous reference allele
	public int consistent_snp = 0; // number of assessed loci where a SNP is observed in at least one individual and genotyping calls are consistent between the trio members
	public int inconsistent_snp = 0; // number of assessed loci where a SNP is observed in at least one individual and genotyping calls are inconsistent
	public int missing_snp_in_parents = 0; // number of inconsistent snps (see above), where parent(s) have a snp but the kid does not while she should
	public int missing_snp_in_kid = 0; // number of inconsistent snps (see above), where kid has a snp but the parents do not while they should
	public int consistent_indels = 0; // *_indels are same as *_snps, see above
	public int consistent_indels_in_mother = 0; // *_indels are same as *_snps, see above
	public int consistent_indels_in_father = 0; // *_indels are same as *_snps, see above
	public int inconsistent_indels  = 0 ;
	public int missing_indels_in_parents  = 0 ;
	public int missing_indels_in_kid  = 0 ;
	public int non_biallelic_snp = 0; // number of variant calls that are not biallelic
	public int non_biallelic_indel = 0; // number of variant calls that are not biallelic
	public long mom_assessed = 0; // number of assessed loci for mother (i.e. passing confidence threshold filter)
	public long dad_assessed = 0;
	public long kid_assessed = 0;
	public long mom_ref = 0; // number of reference calls (out of total assessed)
	public long dad_ref = 0;
	public long kid_ref = 0;
	public long mom_snp = 0; // number of snp calls (out of total assessed)
	public long dad_snp = 0;
	public long kid_snp = 0;
	public long mom_indel = 0; // number of snp calls (out of total assessed)
	public long dad_indel = 0;
	public long kid_indel = 0;
	public long unclassified_events = 0;
	
	public TrioConcordanceRecord add(TrioConcordanceRecord other) {
		this.assessed_loci += other.assessed_loci;
		this.consistent_ref += other.consistent_ref;
		this.consistent_snp += other.consistent_snp;
		this.inconsistent_snp += other.inconsistent_snp;
		this.missing_snp_in_parents += other.missing_snp_in_parents;
		this.missing_snp_in_kid += other.missing_snp_in_kid;
		this.consistent_indels += other.consistent_indels;
		this.consistent_indels_in_mother += other.consistent_indels_in_mother;
		this.consistent_indels_in_father += other.consistent_indels_in_father;
		this.inconsistent_indels += other.inconsistent_indels;
		this.missing_indels_in_parents += other.missing_indels_in_parents;
		this.missing_indels_in_kid += other.missing_indels_in_kid;
		this.non_biallelic_snp += other.non_biallelic_snp;
		this.non_biallelic_indel += other.non_biallelic_indel;
		this.mom_assessed += other.mom_assessed;
		this.dad_assessed += other.dad_assessed;
		this.kid_assessed += other.kid_assessed;
		this.mom_ref += other.mom_ref;
		this.dad_ref += other.dad_ref;
		this.kid_ref += other.kid_ref;
		this.mom_snp += other.mom_snp;
		this.dad_snp += other.dad_snp;
		this.kid_snp += other.kid_snp;
		this.mom_indel += other.mom_indel;
		this.dad_indel += other.dad_indel;
		this.kid_indel += other.kid_indel;
		this.unclassified_events += other.unclassified_events;
		return this;
	}
	
	public int totalSNP() { return consistent_snp + inconsistent_snp + non_biallelic_snp; }
	public int totalIndels() { return consistent_indels + inconsistent_indels + non_biallelic_indel; }
	
	public String toString() {
		StringBuilder b = new StringBuilder();
		
		b.append(String.format("%ntotal assessed in trio: %d%n   reference: %d (%3.2f)%n", assessed_loci, consistent_ref, Utils.percentage(consistent_ref,assessed_loci )) );
		b.append(String.format("   total snp sites: %d%n      consistent snp: %d (%3.2f)%n      multiallelic snp: %d (%3.2f)%n      inconsistent snp: %d (%3.2f)%n", 
				totalSNP(), consistent_snp, Utils.percentage(consistent_snp,totalSNP() ), non_biallelic_snp, Utils.percentage(non_biallelic_snp,totalSNP()), 
				inconsistent_snp, Utils.percentage(inconsistent_snp, totalSNP())  ) );
		
		b.append(String.format("   total indel sites: %d%n      consistent indel: %d (%3.2f)%n      multiallelic indel: %d (%3.2f)%n      inconsistent indels: %d (%3.2f)%n         missing from daughter: %d (%3.2f)%n         missing from both parents: %d (%3.2f)%n", 
				totalIndels(), consistent_indels, Utils.percentage(consistent_indels,totalIndels()), 
				non_biallelic_indel, Utils.percentage(non_biallelic_indel, totalIndels()),
				inconsistent_indels, Utils.percentage(inconsistent_indels, totalIndels()),
				missing_indels_in_kid, Utils.percentage(missing_indels_in_kid, inconsistent_indels),
				missing_indels_in_parents, Utils.percentage(missing_indels_in_parents, inconsistent_indels)
				));

		b.append(String.format("   unclassified (snp+indel): %d%n", unclassified_events));
		
		b.append( String.format("per trio individual:%n   mother:%n      assessed: %d%n      ref: %d%n      snps: %d%n      indels: %d%n", mom_assessed, mom_ref, mom_snp,mom_indel) );
		b.append( String.format("   father:%n      assessed: %d%n      ref: %d%n      snps: %d%n      indels: %d%n", dad_assessed, dad_ref, dad_snp,dad_indel) );
		b.append( String.format("   daughter:%n      assessed: %d%n      ref: %d%n      snps: %d%n      indels: %d%n", kid_assessed, kid_ref, kid_snp,kid_indel) );
		
		return b.toString();
	}
}

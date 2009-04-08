package org.broadinstitute.sting.playground.utils;

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
	public int inconsistent_snp_in_parent = 0; // number of inconsistent snps (see above), where parent(s) have a snp but the kid does not while she should
	public int inconsistent_snp_in_kid = 0; // number of inconsistent snps (see above), where kid has a snp but the parents do not while they should
	public int consistent_indels = 0; // *_indels are same as *_snps, see above
	public int inconsistent_indels  = 0 ;
	public int inconsistent_indels_in_parent  = 0 ;
	public int inconsistent_indels_in_kid  = 0 ;
	public int non_biallelic = 0; // number of variant calls that are not biallelic
	
	public TrioConcordanceRecord add(TrioConcordanceRecord other) {
		this.assessed_loci += other.assessed_loci;
		this.consistent_ref += other.consistent_ref;
		this.consistent_snp += other.consistent_snp;
		this.inconsistent_snp += other.inconsistent_snp;
		this.inconsistent_snp_in_parent += other.inconsistent_snp_in_parent;
		this.inconsistent_snp_in_kid += other.inconsistent_snp_in_kid;
		this.consistent_indels += other.consistent_indels;
		this.inconsistent_indels += other.inconsistent_indels;
		this.inconsistent_indels_in_parent += other.inconsistent_indels_in_parent;
		this.inconsistent_indels_in_kid += other.inconsistent_indels_in_kid;
		this.non_biallelic += other.non_biallelic;
		
		return this;
	}
	
	public int totalSNP() { return consistent_snp + inconsistent_snp + non_biallelic; }
	
	public String toString() {
		return String.format("assessed: %d; reference: %d (%3.2f);  total snp: %d;  consistent snp: %d (%3.2f); multiallelic: %d (%3.2f); " , 
				assessed_loci, consistent_ref, ((double)consistent_ref*100.00)/assessed_loci,totalSNP(), consistent_snp, ((double)consistent_snp*100.0)/totalSNP(),
				non_biallelic, ((double)non_biallelic*100.0)/totalSNP());
	}
}

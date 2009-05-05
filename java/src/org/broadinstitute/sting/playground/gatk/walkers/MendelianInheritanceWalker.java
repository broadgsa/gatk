package org.broadinstitute.sting.playground.gatk.walkers;


import java.util.List;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.rodSAMPileup;
import org.broadinstitute.sting.gatk.refdata.Genotype;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.playground.utils.TrioConcordanceRecord;
import org.broadinstitute.sting.utils.cmdLine.Argument;

public class MendelianInheritanceWalker  extends RefWalker<TrioConcordanceRecord, TrioConcordanceRecord> {

	@Argument(fullName="consensus_cutoff", shortName="XC", doc="consensus cutoff", required=true ) public Double CONS_CUTOFF;
	@Argument(fullName="snp_cutoff", shortName="XS", doc="snp cutoff", required=true ) public Double SNP_CUTOFF;
	@Argument(fullName="indel_cutoff", shortName="XI", doc="indel cutoff", required=true ) public Double INDEL_CUTOFF;
	@Argument(fullName="log_concordant", shortName="LC",doc="If set, all trio-concordant sites will be logged at level INFO") public boolean LOG_CONCORDANT; 
	@Argument(fullName="log_discordant", shortName="LD",doc="If set, all trio-discordant sites will be logged at level INFO") public boolean LOG_DISCORDANT; 
	
	private static Logger logger = Logger.getLogger(MendelianInheritanceWalker.class);	
	private final static String star = new String("*");
	
	@Override
	public TrioConcordanceRecord map(RefMetaDataTracker rodData, char ref, LocusContext context) {
		
		TrioConcordanceRecord t = new TrioConcordanceRecord();
		
//		String outLine = new String(context.getLocation() + " REF: "+ref + " RODS:" + rodData.getAllRods().size());
		
		Genotype mom = (rodSAMPileup)rodData.lookup("mother", null);
		Genotype dad = (rodSAMPileup)rodData.lookup("father", null);
		Genotype kid = (rodSAMPileup)rodData.lookup("daughter", null);
		
		if ( hasCall(mom)) {
			t.mom_assessed = 1;
			if ( mom.isReference() ) t.mom_ref = 1;
			else if ( mom.isSNP() ) t.mom_snp = 1;
			else if ( mom.isIndel() ) t.mom_indel = 1;
		}

		if ( hasCall(dad)) {
			t.dad_assessed = 1;
			if ( dad.isReference() ) t.dad_ref = 1;
			else if ( dad.isSNP() ) t.dad_snp = 1;
			else if ( dad.isIndel() ) t.dad_indel = 1;
		}
		
		if ( hasCall(kid)) {
			t.kid_assessed = 1;
			if ( kid.isReference() ) t.kid_ref = 1;
			else if ( kid.isSNP() ) t.kid_snp = 1;
			else if ( kid.isIndel() ) t.kid_indel = 1;
		}
		

		//		if ( hasCall(mom) && mom.isIndel() ) System.out.println("GOT INDEL: "+mom.toString());
		
		if (( t.mom_assessed + t.dad_assessed + t.kid_assessed) != 3 ) return t; // at least one person is not called; nothing to do, bail out
			
		// proceed only if we got confident calls
			
        t.assessed_loci = 1;
			
		
		// ok, we got a call and it is biallelic (at most)
		
		if ( mom.isReference() && dad.isReference() && kid.isReference() ) { // got all refs, we are done
			t.consistent_ref = 1;
			return t;
		}

		// by now we know that there's a SNP or an indel: at least one of the individuals is non-ref
		
		String kid_allele_1 = kid.getFWDAlleles().get(0);
		String kid_allele_2 = kid.getFWDAlleles().get(1);
		List<String> mom_alleles = mom.getFWDAlleles();
		List<String> dad_alleles = dad.getFWDAlleles();
		
		if ( ! mom.isIndel() && ! dad.isIndel() && ! kid.isIndel() ) {
			// individuals have SNPs (any member of the trio is ok, we'll check for consistency in a sec), but no indels at the site:
		
			// warning: no special processing for X/Y chromosomes yet; not an issue for daughter

			if ( ! mom.isBiallelic() || ! dad.isBiallelic() || ! kid.isBiallelic() ) {
				t.non_biallelic_snp = 1; // we currently ignore non-biallelic sites
				return t;
			}
			
			
			if ( mom_alleles.contains(kid_allele_1) && dad_alleles.contains(kid_allele_2) ||
					mom_alleles.contains(kid_allele_2) && dad_alleles.contains(kid_allele_1) ) {
				t.consistent_snp = 1;

				if ( LOG_CONCORDANT ) {
					logger.info("consistent SNP at "+context.getLocation() + 
							"("+ref+") " + mom_alleles.get(0)+"/" +mom_alleles.get(1) + "  " +
							dad_alleles.get(0)+"/" +dad_alleles.get(1) + "  " +
							kid_allele_1+"/" +kid_allele_2  
							);
				}
			} else {
				// we are inconsistent. let's see what happened:
				if ( kid.isSNP() && ! mom.isSNP() && ! dad.isSNP() ) t.missing_snp_in_parents = 1;
				if ( ! kid.isSNP() && ( mom.isSNP() && mom.isHom() || dad.isSNP() && dad.isHom() ) ) t.missing_snp_in_kid = 1;
				
				t.inconsistent_snp = 1;
				if ( LOG_DISCORDANT ) {
					logger.info("INconsistent SNP at "+context.getLocation() + 
							"("+ref+")  mom:" + mom_alleles.get(0)+"/" +mom_alleles.get(1) + " dad: " +
							dad_alleles.get(0)+"/" +dad_alleles.get(1) + " kid: " +
							kid_allele_1+"/" +kid_allele_2  
							);
				}
			}
			return t;
		}

		if ( ! mom.isSNP() && ! dad.isSNP() && ! kid.isSNP() ) {
			// individuals have indels (any member of the trio is ok, we'll check for consistency in a seq), but no indels at the site:
			
			// warning: no special processing for X/Y chromosomes yet; not an issue for daughter

			if ( ! mom.isBiallelic() || ! dad.isBiallelic() || ! kid.isBiallelic() ) {
				t.non_biallelic_indel = 1; // we currently ignore non-biallelic sites
				return t;
			}

			if ( mom_alleles.contains(kid_allele_1) && dad_alleles.contains(kid_allele_2) ||
						mom_alleles.contains(kid_allele_2) && dad_alleles.contains(kid_allele_1) ) {
					t.consistent_indels = 1;
				
					if ( LOG_CONCORDANT ) {
						logger.info("consistent INDEL at "+context.getLocation() + 
								"("+ref+")  mom:" + genotypeString(mom)+ " dad: " +
								genotypeString(dad) + " kid: " +
								genotypeString(kid)  
								);
					}
			} else {
				  if ( kid.isIndel() && ! mom.isIndel() && ! dad.isIndel() ) t.missing_indels_in_parents = 1;
				  if ( ! kid.isIndel() && ( mom.isIndel() && mom.isHom() || dad.isIndel() && dad.isHom() ) ) t.missing_indels_in_kid = 1;

				  t.inconsistent_indels = 1;
				  
				  if ( LOG_DISCORDANT ) {
					  logger.info("INconsistent INDEL at "+context.getLocation() + 
							  "("+ref+")  mom:" + genotypeString(mom)+ " dad: " +
							  genotypeString(dad) + " kid: " +
							  genotypeString(kid)  
								);
				  }
			}
			return t;
		}
		
		t.unclassified_events = 1;
		
		return t;
		
	}

	protected String alleleString(Genotype g, int n) {
		if ( g.getFWDAlleles().get(n).length() == 0 ) return star;
		return g.getFWDAlleles().get(n);
	}
	
	protected String genotypeString(Genotype g) {
		return alleleString(g, 0) +"/"+alleleString(g, 1);
	}

	@Override
	public TrioConcordanceRecord reduce(TrioConcordanceRecord value, TrioConcordanceRecord sum) {
		
		return sum.add(value);
	}

	@Override
	public TrioConcordanceRecord reduceInit() {
		
		return new TrioConcordanceRecord();
	}

	boolean hasCall(Genotype g) {
		if ( g == null ) return false; // there's no call if there's no rod data available, duh!
		if ( g.isReference() ) return g.getConsensusConfidence() >= CONS_CUTOFF ;
		else {
			if ( g.isSNP() ) return g.getVariantConfidence() >= SNP_CUTOFF;
			else return Math.max(g.getConsensusConfidence(),g.getVariantConfidence()) >= INDEL_CUTOFF;
		}
		
	}
	
	/*
	protected String shortLine(Genotype av) {
		
		if ( av == null )  return new String( "<NULL>");
		
		if ( av.isReference() ) return new String ("*");
		
		List<String> alleles = av.getFWDAlleles();
		
		if ( av.isSNP() ) {
            if  ( alleles.get(0).charAt(0) == av.getRef() ) return alleles.get(1);
            else return alleles.get(0);
		}
		if ( av.isInsertion() ) {
            if ( alleles.get(0).length() == 0 ) return new String('+'+alleles.get(1));
            else return new String('+'+alleles.get(0));
		}
		if ( av.isDeletion() ) {
            if ( alleles.get(0).length() == 0 ) return new String ('-'+alleles.get(0));
            else return new String('-'+alleles.get(1));
		}
		if ( av.isIndel() ) {
			return new String('?'+alleles.get(0));
		}
		return new String("<missing data!!!>");
	}
	*/
}

package org.broadinstitute.sting.playground.gatk.walkers;


import java.util.List;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.rodSAMPileup;
import org.broadinstitute.sting.gatk.refdata.AllelicVariant;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.playground.utils.TrioConcordanceRecord;
import org.broadinstitute.sting.utils.cmdLine.Argument;

public class MendelianInheritanceWalker  extends RefWalker<TrioConcordanceRecord, TrioConcordanceRecord> {

	@Argument(fullName="conf_cutoff", shortName="X",required=true ) public Double CONF_CUTOFF;
	
	private static Logger logger = Logger.getLogger(MendelianInheritanceWalker.class);	
	
	
	@Override
	public TrioConcordanceRecord map(RefMetaDataTracker rodData, char ref, LocusContext context) {
		
		TrioConcordanceRecord t = new TrioConcordanceRecord();
		
//		String outLine = new String(context.getLocation() + " REF: "+ref + " RODS:" + rodData.getAllRods().size());
		
		AllelicVariant mom = (rodSAMPileup)rodData.lookup("mother", null);
		AllelicVariant dad = (rodSAMPileup)rodData.lookup("father", null);
		AllelicVariant kid = (rodSAMPileup)rodData.lookup("daughter", null);
		
		if ( ! hasCall(mom)  || ! hasCall(dad) || ! hasCall(kid) ) return t; // at least one person is not called; nothing to do, bail out
			
		// proceed only if we got confident calls
			
        t.assessed_loci = 1;
			
		if ( ! mom.isBiallelic() || ! dad.isBiallelic() || ! kid.isBiallelic() ) {
			t.non_biallelic = 1; // we currently ignore non-biallelic sites
			return t;
		}
		
		// ok, we got a call and it is biallelic (at most)
		
		if ( mom.isReference() && dad.isReference() && kid.isReference() ) { // got all refs, we are done
			t.consistent_ref = 1;
			return t;
		}

		// by now we know that there's a SNP
		
		String kid_allele_1 = kid.getGenotype().get(0);
		String kid_allele_2 = kid.getGenotype().get(1);
		List<String> mom_alleles = mom.getGenotype();
		List<String> dad_alleles = dad.getGenotype();
		
		// kid must have one allele from mom and one allele from dad if it's not X chromosome!
		
		if ( mom_alleles.contains(kid_allele_1) && dad_alleles.contains(kid_allele_2) ||
				mom_alleles.contains(kid_allele_2) && dad_alleles.contains(kid_allele_1) ) {
			t.consistent_snp = 1;
			
			logger.info("consistent SNP at "+context.getLocation() + 
					"("+ref+") " + mom_alleles.get(0)+"/" +mom_alleles.get(1) + "  " +
					dad_alleles.get(0)+"/" +dad_alleles.get(1) + "  " +
					kid_allele_1+"/" +kid_allele_2  
					);
			return t;
	    }
		
		// snp is inconsistent
		
		t.inconsistent_snp = 1;
		logger.info("INconsistent SNP at "+context.getLocation() + 
				"("+ref+") " + mom_alleles.get(0)+"/" +mom_alleles.get(1) + "  " +
				dad_alleles.get(0)+"/" +dad_alleles.get(1) + "  " +
				kid_allele_1+"/" +kid_allele_2  
				);
		
		return t;
		
	}

	@Override
	public TrioConcordanceRecord reduce(TrioConcordanceRecord value, TrioConcordanceRecord sum) {
		
		return sum.add(value);
	}

	@Override
	public TrioConcordanceRecord reduceInit() {
		
		return new TrioConcordanceRecord();
	}

	boolean hasCall(AllelicVariant av) {
		if ( av == null ) return false; // there's no call if there's no rod data available, duh!
		if ( av.isReference() ) return av.getConsensusConfidence() >= CONF_CUTOFF ;
		else return av.getVariationConfidence() >= CONF_CUTOFF;
		
	}
	
	protected String shortLine(AllelicVariant av) {
		
		if ( av == null )  return new String( "<NULL>");
		
		if ( av.isReference() ) return new String ("*");
		
		if ( av.isSNP() ) return av.getAltBasesFWD();
		if ( av.isInsertion() ) return new String('+'+av.getAltBasesFWD());
		if ( av.isDeletion() ) return new String('-'+av.getRefBasesFWD());
		if ( av.isIndel() ) return new String('?'+av.getAltBasesFWD());
		return new String("<missing data!!!>");
	}
	
}

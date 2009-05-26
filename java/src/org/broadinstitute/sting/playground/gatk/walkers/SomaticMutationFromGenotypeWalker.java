package org.broadinstitute.sting.playground.gatk.walkers;


import java.io.BufferedOutputStream;
import java.io.FileWriter;
import java.io.OutputStream;

import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.LocusContext;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.gatk.refdata.Genotype;
import org.broadinstitute.sting.gatk.walkers.RefWalker;
import org.broadinstitute.sting.gatk.walkers.Requires;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.playground.utils.GenotypingCallStats;
import org.broadinstitute.sting.utils.GenotypeUtils;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.cmdLine.Argument;


//@Requires(value=DataSource.REFERENCE,referenceMetaData={@RMD(name="mother",type=rodSAMPileup.class),
//                                                        @RMD(name="father",type=rodSAMPileup.class),
//                                                        @RMD(name="daughter",type=rodSAMPileup.class)})
public class SomaticMutationFromGenotypeWalker  extends RefWalker<SomaticMutationRecord, SomaticMutationRecord> {

	@Argument(fullName="point_consensus_cutoff", shortName="XPC", doc="confidence cutoff for consensus in point genotype", required=true ) public Double POINT_CONS_CUTOFF;
	@Argument(fullName="point_variant_cutoff", shortName="XPV", doc="confidence cutoff for variant (snp) in point genotype", required=true ) public Double POINT_VAR_CUTOFF;
	@Argument(fullName="indel_consensus_cutoff", shortName="XIC", doc="confidence cutoff for consensus in indel genotype", required=true ) public Double INDEL_CONS_CUTOFF;
	@Argument(fullName="indel_variant_cutoff", shortName="XIV", doc="confidence cutoff for variant (snp) in indel genotype", required=true ) public Double INDEL_VAR_CUTOFF;
	@Argument(fullName="default_reference_calls",shortName="DRC",
			doc="If set,  any position where the requested genotype subtype (variant_type) is NOT explicitly specified, while the other is provided, is considered to be an implicit confident 'reference' (no-indel or no-snp) call") 
			public boolean defCalls;
	@Argument(fullName="variant_type",
		      shortName="VT",
		      doc="Look for variants of specified type (POINT mutations or INDELs)",
		      required=true) 
		      public String VTYPE_STR; 
	@Argument(fullName="bed_out", shortName="BED", doc="Write somatic variants into the specified file in BED format",required=false) public java.io.File BED_OUT;
	
	private java.io.Writer bed_stream = null;
	private static Logger logger = Logger.getLogger(MendelianInheritanceWalker.class);	
	private final static String star = new String("*");
	private GenotypeUtils.VariantType VARIANT_TYPE; 
	
	@Override
	public SomaticMutationRecord map(RefMetaDataTracker rodData, char ref, LocusContext context) {
				
//		String outLine = new String(context.getLocation() + " REF: "+ref + " RODS:" + rodData.getAllRods().size());
		
		ReferenceOrderedDatum rodNormal = rodData.lookup("normal", null);
		ReferenceOrderedDatum rodTumor = rodData.lookup("tumor", null);
		ReferenceOrderedDatum rodRefseq = rodData.lookup("refseq", null);
		
		Genotype normal = GenotypeUtils.extractGenotype(rodNormal,VARIANT_TYPE, defCalls);
		Genotype tumor = GenotypeUtils.extractGenotype(rodTumor,VARIANT_TYPE,defCalls);

		SomaticMutationRecord r = new SomaticMutationRecord();
		if ( rodRefseq == null ) return r;
		out.println(context.getLocation() + " " + rodRefseq.toString() );
		if ( true ) return r;
			

		
		assessGenotype(normal,r.normal);
		assessGenotype(tumor, r.tumor);

		
		
		
		if ( r.normal.covered == 1 && r.tumor.covered == 1 ) r.mut.covered = 1;
		if ( r.normal.assessed == 1 && r.tumor.assessed == 1 ) r.mut.assessed = 1;
		else {
			if ( r.tumor.variant == 1 ) {
				if ( r.normal.covered == 0 ) r.tumor_variant_without_normal = 1; // normal is not even covered
				else if ( r.normal.assessed == 0 ) r.tumor_variant_without_confident_normal = 1; // normal is covered but can not make confident call
			}
			if ( r.normal.variant == 1 ) {
				if ( r.tumor.covered == 0 ) r.normal_variant_without_tumor = 1; // normal is not even covered
				else if ( r.tumor.assessed == 0 ) r.normal_variant_without_confident_tumor = 1; // normal is covered but can not make confident call
			}
			return r;
		}
		
		// both tumor and normal are assessed
		
		if ( r.normal.variant == 0 && r.tumor.variant == 0 ) {
			r.mut.ref = r.mut.consistent_ref = 1;
			return r;
		}
		
		// at least one sample has a confident variant call
		
		String n_allele1 = normal.getFWDAlleles().get(0);
		String n_allele2 = normal.getFWDAlleles().get(1);
		String t_allele1 = tumor.getFWDAlleles().get(0);
		String t_allele2 = tumor.getFWDAlleles().get(1);
		
		if ( r.normal.variant == 1 && r.tumor.variant == 1 ) {
			// germline variant
			r.germline = 1;
			if ( !  ( n_allele1.equals( t_allele1 ) && n_allele2.equals( t_allele2 ) || 
					n_allele1.equals( t_allele2 ) && n_allele2.equals( t_allele1 ) ) ) r.non_matching_germline = 1;
			return r;
		}
		
		if ( r.normal.variant == 1 && r.tumor.variant == 0 ) {
			// variant in normal, but not in tumor, what the hell?
			r.lost_variant = 1;
			return r;
		}
		
		// last possibility: variant in tumor but not in normal (both are assessed!)
		
		r.mut.variant = 1;
		
		if ( bed_stream != null ) {
			long stop = normal.getLocation().getStop();
			if ( stop < normal.getLocation().getStart() + 1 ) stop = normal.getLocation().getStart()+1;
			try {
				bed_stream.append(normal.getLocation().getContig()+"\t"+
												  normal.getLocation().getStart() + "\t" +
												  stop +"\t"+
												  genotypeString(tumor));
			} catch ( Exception e ) {
				throw new StingException("Failed to write into BED output "+BED_OUT+": " + e.getMessage() );
			}
		}
		
		return r;
	}
	
	

	/** Takes a single genotype object and returns properly filled new assessment object (covered/assessed/ref/variant set to 0/1 
	 * according to what the genotype says)  
	 * @param g
	 * @return
	 */
	protected GenotypingCallStats assessGenotype(Genotype g, GenotypingCallStats stats) {
		
		if ( g != null ) stats.covered = 1;
		
		if ( hasCall(g)) {
			stats.assessed = 1;
			if ( g.isReference() ) stats.ref = 1;
			else {
				stats.variant = 1;
				if ( ! g.isBiallelic() ) stats.non_biallelic_variant = 1;
			}
		}
		return stats;
	}
		
	
	protected String alleleString(Genotype g, int n) {
		if ( g.getFWDAlleles().get(n).length() == 0 ) return star;
		return g.getFWDAlleles().get(n);
	}
	
	protected String genotypeString(Genotype g) {
		return alleleString(g, 0) +"/"+alleleString(g, 1);
	}

	@Override
	public SomaticMutationRecord reduce(SomaticMutationRecord value, SomaticMutationRecord sum) {
		return sum.add(value);
	}

	@Override
	public void initialize() {
		super.initialize();
		if ( BED_OUT != null ) {
			try {
				bed_stream = new java.io.BufferedWriter( new FileWriter(BED_OUT) );
			} catch ( Exception e ) {
				throw new StingException("Failed to initialize BED output "+BED_OUT+": " + e.getMessage() );
			}
		}
		VARIANT_TYPE = GenotypeUtils.VariantType.valueOf(VTYPE_STR.toUpperCase());

	}

	
	@Override
	public SomaticMutationRecord reduceInit() {
		
		return new SomaticMutationRecord();
	}

	boolean hasCall(Genotype g) {
		if ( g == null ) return false; // there's no call if there's no rod data available, duh!

		if ( g.isReference() ) {
			if ( g.isPointGenotype() ) return g.getConsensusConfidence() >= POINT_CONS_CUTOFF ;
			else return g.getConsensusConfidence() >= INDEL_CONS_CUTOFF ;
		}
		else { // it's a variant
			if ( g.isPointGenotype() ) return g.getVariantConfidence() >= POINT_VAR_CUTOFF ;
			else return g.getVariantConfidence() >= INDEL_VAR_CUTOFF ;
		}
		
	}
	
	public void onTraversalDone(SomaticMutationRecord result) {
		if ( bed_stream != null ) {
			try {
				bed_stream.close();
			} catch ( Exception e ) {
				throw new StingException("Failed to close BED output file "+BED_OUT+": " + e.getMessage() );
			} 
		}
		super.onTraversalDone(result);
	}
	
}

class SomaticMutationRecord {
	GenotypingCallStats normal;
	GenotypingCallStats tumor;
	GenotypingCallStats mut;
	int germline = 0;
	int non_matching_germline = 0;
	int tumor_variant_without_normal = 0;
	int tumor_variant_without_confident_normal = 0;
	int normal_variant_without_tumor = 0;
	int normal_variant_without_confident_tumor = 0;
	int lost_variant = 0;
	
	public SomaticMutationRecord() { 
		normal = new GenotypingCallStats();
		tumor = new GenotypingCallStats();
		mut = new GenotypingCallStats();
	}
	
	public SomaticMutationRecord add(SomaticMutationRecord other) {
		normal.add(other.normal);
		tumor.add(other.tumor);
		mut.add(other.mut);
		
		germline += other.germline;
		non_matching_germline += other.non_matching_germline;
		tumor_variant_without_normal += other.tumor_variant_without_normal;
		tumor_variant_without_confident_normal += other.tumor_variant_without_confident_normal;
		normal_variant_without_tumor += other.normal_variant_without_tumor;
		normal_variant_without_confident_tumor += other.normal_variant_without_confident_tumor;
		lost_variant += other.lost_variant;
		
		return this; 
	}
}


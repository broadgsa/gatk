package org.broadinstitute.sting.utils;

import java.util.Arrays;
import java.util.List;

import org.broadinstitute.sting.gatk.refdata.Genotype;
import org.broadinstitute.sting.gatk.refdata.GenotypeList;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;

/** Holds useful utility methods and auxiliary default classes for working with Genotype objects
 * 
 * @author asivache
 *
 */
public class GenotypeUtils {
	
	public static enum VariantType {
		POINT, INDEL;
		VariantType() {}
	}

	/** This method accepts rods that implement either Genotype or GenotypeList interface (all others will result in an exception). Variant 
	 * (Genotype object) of the specified type (point mutation or indel) will be extracted from GenotypeList rod if such variant exists, or the rod itself 
	 * will be typecasted and returned back if it implements Genotype and represents the specified variant type. If the last argument is false, then 
	 * null will be returned in all other cases. If the last argument is true and either a) rod is a GenotypeList that lacks a call of the specified type, but call
	 * of the other type was made, or b) rod is a Genotype and its type is difference from variant_type, then it will be assumed that the implicit ref/ref call
	 * of the specified type also exists at this genomic position, and new object representing such default call will be returned. If rod argument
	 * is null, then this method safely (and silently) returns null.
	 *  
	 * @param rod GenotypeList or Genotype to extract requested call from/upgrade to requested call
	 * @param variant_type type of the variant to extract (POINT mutations or INDEL)
	 * @param assume_implicit_ref_calls true if presence of only a call of different type means that ref/ref call of the requested type is implicitly present 
	 * @return Genotyping call of the requested type or null if no explicit or implicit call was made.
	 */
	public static Genotype extractGenotype(ReferenceOrderedDatum rod, VariantType variant_type, boolean assume_implicit_ref_calls) {
		
		if ( rod == null ) return null;

		if ( rod instanceof GenotypeList ) {

			GenotypeList rod_gl = (GenotypeList)rod;
			switch ( variant_type ) {
			case POINT :
				if ( rod_gl.hasPointGenotype() ) return rod_gl.getPointGenotype();
				else {
					if ( assume_implicit_ref_calls && rod_gl.hasIndelGenotype() ) throw new StingException("Default (reference) implicit POINT genotype is not implemented yet");
					else return null;
				}
			case INDEL:				
				if ( rod_gl.hasIndelGenotype() ) return rod_gl.getIndelGenotype();
				else {
					if ( assume_implicit_ref_calls && rod_gl.hasPointGenotype() ) return new DefaultIndelGenotype(rod_gl.getLocation());
					else return null;
				}
			default: throw new StingException("Unrecognized variant type: "+variant_type);
			}
			
		}	else if ( rod instanceof Genotype ) {
			
			Genotype rod_g = (Genotype)rod;
			switch ( variant_type ) {
			case POINT:
				if ( rod_g.isIndelGenotype() ) {
					if ( assume_implicit_ref_calls  ) throw new StingException("Default (reference) implicit POINT genotype is not implemented yet");
					else return null;
				}	else return rod_g;
			case INDEL: 
				if ( rod_g.isPointGenotype() ) {
					if ( assume_implicit_ref_calls  ) return new DefaultIndelGenotype(rod_g.getLocation());
					else return null;
				}	else return rod_g;
			default: throw new StingException("Unrecognized variant type: "+variant_type);
			} 
			
		}
		else throw new StingException("track "+rod.getName()+" is not a Genotype or GenotypeList");
	}

	
	/** This class represents a "default" indel-type genotype with homozygous reference (i.e. confidently no indel)
	 * call. All the interface methods are implemented and return consistent values. Use this class when working with
	 * genotyping data where absence of explicit indel call actually means that no evidence for an indel was observed 
	 *  (SAM pileup is one such example).
	 * @author asivache
	 *
	 */
	public static class DefaultIndelGenotype implements Genotype {
		private GenomeLoc location;
		private static final List<String> alleles = Arrays.asList("","") ;
		private int confidence ;

		/** Creates default indel genotype (ref/ref = no indel) at the specified position
		 * with default consensus confidence of 1000 (variant confidence is 0).
		 * @param l reference position to associate the genotyping call with
		 */
		public DefaultIndelGenotype(GenomeLoc l) {
			this(l,1000);
		}
		
		/** Creates ref/ref (i.e. absense of indel) indel genotype at the specified position
		 * with the specified consensus confidence (variant confidence is 0).
		 * @param l reference position to associate the genotyping call with
		 */
		public DefaultIndelGenotype(GenomeLoc l, int confidence) {
			location = l;
			this.confidence = confidence;
		}

		@Override
		public double getConsensusConfidence() {
			return confidence;
		}

		@Override
		public List<String> getFWDAlleles() {
			return alleles;
		}

		@Override
		public String getFWDRefBases() {
			return alleles.get(0);
		}

		@Override
		public GenomeLoc getLocation() {
			return location;
		}

		@Override
		public char getRef() {
			return '*';
		}

		@Override
		public double getVariantConfidence() {
			return 0.0;
		}

		@Override
		public boolean isBiallelic() {
			return false;
		}

		@Override
		public boolean isDeletion() {
			return false;
		}

		@Override
		public boolean isHet() {
			return false;
		}

		@Override
		public boolean isHom() {
			return true;
		}

		@Override
		public boolean isIndel() {
			return false;
		}

		@Override
		public boolean isIndelGenotype() {
			return true;
		}

		@Override
		public boolean isInsertion() {
			return false;
		}

		@Override
		public boolean isPointGenotype() {
			return false;
		}

		@Override
		public boolean isReference() {
			return true;
		}

		@Override
		public boolean isSNP() {
			return false;
		}

		@Override
		public int compareTo(ReferenceOrderedDatum o) {
			return location.compareTo(o.getLocation());
		}
		
	}

}

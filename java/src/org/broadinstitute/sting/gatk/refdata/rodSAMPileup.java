package org.broadinstitute.sting.gatk.refdata;


import java.util.*;

import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedDatum;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.gatk.refdata.AllelicVariant;

/**
 *  This class wraps Maq/samtools allele calls from pileup format and presents them as a ROD.<br>
 *  
 * Example format:<br>
 *   for SNP:<br>
 *      [chr]       [pos] [ref]  [consensus allele(s)]    [consensus confidence]  [snp confidence]   [max mapping qual]  [num reads in the pile] <br>
 *     chrX       466     T                     Y                                   170                                 170                          88                              32 ... (piles of read bases  and quals follow) <br>
 *     <br>
 *   for indel: <br>
 *    [chr]         [pos]     [always *]  [alleles]         [?]        [?]     [?]     [num reads in the pile] ... <br>     
 *     chrX       141444  *                +CA/+CA     32      468     255     25                                        +CA     *       5       2       12      6 <br>
 * User: asivache
 * Date: Apr 03, 2009
 * Time: 2:58:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class rodSAMPileup extends ReferenceOrderedDatum implements AllelicVariant {
	private static final int NO_VARIANT = -1;
	private static final int SNP_VARIANT = 0;
	private static final int INSERTION_VARIANT = 1;
	private static final int DELETION_VARIANT = 2;
	private static final int INDEL_VARIANT = 3;
	
	// allocate once and don't ever bother creating them again:
	private static final String baseA = new String("A");
	private static final String baseC = new String("C");
	private static final String baseG = new String("G");
	private static final String baseT = new String("T");
	private static final String emptyStr = new String(); // we will use this for "reference" allele in insertions

    protected GenomeLoc loc;       // genome location of SNP
                                // Reference sequence chromosome or scaffold
                                // Start and stop positions in chrom


    protected char refBaseChar; // what we have set for the reference base (is set to a '*' for indel!)
    protected String refBases;        // the reference base according to NCBI, in the dbSNP file
    protected String observedString; // store the actual string representation of observed alleles

    protected List<String> observedAlleles = null;    // The sequences of the observed alleles from rs-fasta files
    protected int varType = NO_VARIANT;
    protected int ploidy = 2; // how many allelic variants we get?
    protected int nNonref = 0; // number of non-reference alleles
    protected int eventLength = 0;
    
    protected double consensusScore;
    protected double variantScore;
    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
    public rodSAMPileup(final String name) {
        super(name);
    }

    @Override
    public void parseLine(final String[] parts) {
//    	          0          1             2         3                  4         5            6         7
//    	*     chrX     466           T         Y                170      170       88      32 ... (piles of read bases  and quals follow)
//    	 *    chrX    141444     *     +CA/+CA       32       468     255     25     +CA     *       5       2       12      6
       try {
    	   
            String contig = parts[0];
            long start = Long.parseLong(parts[1]) ;
            
            consensusScore = Double.parseDouble(parts[4]);
            variantScore = Double.parseDouble(parts[5]);

            refBaseChar = Character.toUpperCase(parts[2].charAt(0));
            
            parts[3] = parts[3].toUpperCase();
            observedString = parts[3];
            
            observedAlleles = new ArrayList<String>(2);
            
            if ( refBaseChar == '*' ) {
            	parseIndels(parts[3]) ;
            	if ( varType == DELETION_VARIANT ) loc = new GenomeLoc(contig, start, start+eventLength);
            	else loc = new GenomeLoc(contig, start, start); // if it's not a deletion and we are biallelic, this got to be an insertion; otherwise the state is inconsistent!!!! 
            }
            else {
            	// if the variant is a SNP or a reference base (i.e. no variant at all)
            	assert parts[3].length() == 1 : "non-indel genotype is expected to be represented by a single letter";
            	refBases = parts[2].toUpperCase();
            	eventLength = 1;
            	loc = new GenomeLoc(contig, start, start+1);

            	char ch = parts[3].charAt(0);
            	
            	switch ( ch ) {
            	case 'A':  observedAlleles.add(baseA); observedAlleles.add(baseA); break;
            	case 'C':  observedAlleles.add(baseC); observedAlleles.add(baseC); break;
            	case 'G':  observedAlleles.add(baseG); observedAlleles.add(baseG); break;
            	case 'T':   observedAlleles.add(baseT); observedAlleles.add(baseT); break;
            	case 'M': observedAlleles.add(baseA); observedAlleles.add(baseC); break;
            	case 'R': observedAlleles.add(baseA); observedAlleles.add(baseG); break;
            	case 'W': observedAlleles.add(baseA); observedAlleles.add(baseT); break;
            	case 'S': observedAlleles.add(baseC); observedAlleles.add(baseG); break;
            	case 'Y': observedAlleles.add(baseC); observedAlleles.add(baseT); break;
            	case 'K': observedAlleles.add(baseG); observedAlleles.add(baseT); break;
            	}
            	if ( observedAlleles.get(0).charAt(0) == refBaseChar && observedAlleles.get(1).charAt(0) == refBaseChar ) varType = NO_VARIANT;
            	else {
            		// we know that at least one allele is non-ref;
            		// if one is ref and the other is non-ref, or if both are non ref but they are the same (i.e.
            		// homozygous non-ref), we still have 2 allelic variants at the site (e.g. one ref and one nonref)
            		varType = SNP_VARIANT;
                	if ( observedAlleles.get(0).charAt(0) == refBaseChar || 
                			observedAlleles.get(1).charAt(0) == refBaseChar ||
                			observedAlleles.get(0) == observedAlleles.get(1)
                	   ) nNonref = 1;
                	else nNonref = 2; // if both observations differ from ref and they are not equal to one another, then we get multiallelic site...
            	}
            }
        } catch ( RuntimeException e ) {
            System.out.printf("  Exception caught during parsing Pileup line:  %s%n", Utils.join(" <=> ", parts));
            throw e;
        }
        if ( nNonref > 1 ) System.out.println("SAM pileup: WARNING: multi-allelic variant at "+ loc.toString());
    }

    
    private void parseIndels(String genotype) {
    	String [] obs = genotype.split("/"); // get observations, now need to tinker with them a bit
    	
    	// if reference allele is among the observed alleles, we will need to take special care of it since we do not have direct access to the reference;
    	// if we have an insertion, the "reference" allele is going to be empty; if it it is a deletion, we will deduce the "reference allele" bases
    	// from what we have recorded for the deletion allele (e.g. "-CAC")
    	boolean hasRefAllele = false; 
    	
    	for ( int i = 0 ; i < obs.length ; i++ ) {
    		if ( obs[i].length() == 1 && obs[i].charAt(0) == '*'  ) {
    			hasRefAllele = true;
    			continue; // we will fill reference allele later
    		}
    		
    		String varBases = obs[i].substring(1).toUpperCase(); 

    		switch ( obs[i].charAt(0) )  {
    		case '+': 
    			if ( varType != NO_VARIANT && varType != INSERTION_VARIANT )  varType = INDEL_VARIANT;
    			else varType = INSERTION_VARIANT;
    			refBases = emptyStr;
    			break;
    		case '-' :
    			if ( varType != -1 && varType != DELETION_VARIANT )  varType = INDEL_VARIANT;
    			else varType = DELETION_VARIANT;
    			refBases = varBases; // remember what was deleted, this will be saved as "reference allele"
    			break;
    		default: throw new RuntimeException("Can not interpret observed indel allele record: "+genotype);
    		}
    		observedAlleles.add(varBases);
			eventLength = obs[i].length() - 1; // inconsistent for non-biallelic indels!!
    	}
    	if ( hasRefAllele ) {
    		// we got at least one ref. allele (out of two recorded)
    		if ( varType == NO_VARIANT ) { // both top theories are actually ref allele; 
    			observedAlleles.add(emptyStr);
    			observedAlleles.add(emptyStr); // no inserted/deleted bases at the site!
    			nNonref = 0; // no observations of non-reference allele at all
    		} else {
    			nNonref = 1; // hasRefAllele = true, so one allele was definitely ref, hence there is only one left
    			
    			 // whether we have insertion or deletion, one allele (ref for insertion, or alt for deletion) is empty;
    			// the one that contain actual bases (alt for insertion or ref for deletion) was already filled above:
    			observedAlleles.add(emptyStr);
    		}
    	} else {
    		// we observe two non-ref alleles; they better be the same variant, otherwise the site is not bi-allelic and at the moment we
    		// fail to set data in a consistent way.. (the check for INDEL_VARIANT ensures that recorded variants are indeed both insertions
    		// or both deletions as compared to +ACC/-ACC which would still have the same bases (no matter how crazy and improbable
    		// such event would be)
    		if ( observedAlleles.get(0).equals(observedAlleles.get(1)) && varType != INDEL_VARIANT ) nNonref = 1;
    		else nNonref = 2;
    	}
    	// DONE with indels

    }
    
    
    public GenomeLoc getLocation() { return loc; }

    /** Returns bases in the reference allele as a String. String can be empty (as in insertion into
     * the reference), can contain a single character (as in SNP or one-base deletion), or multiple characters
     * (for longer indels).
     *
     * @return reference allele, forward strand
     */
    public String getRefBasesFWD() {
           return refBases;
    }

    /**
     * Returns reference (major) allele base for a SNP variant as a character; should throw IllegalStateException
     * if variant is not a SNP.
     *
     * @return reference base on the forward strand
     */
    public char getRefSnpFWD() throws IllegalStateException {
        if ( isIndel() ) throw new IllegalStateException("Variant is not a SNP");
        return refBaseChar;
    }

    @Override
    public List<String> getGenotype()  {
        return observedAlleles;
    }

    // ----------------------------------------------------------------------
    //
    // What kind of variant are we?
    //
    // ----------------------------------------------------------------------
    public boolean isSNP() { return varType == SNP_VARIANT ; }
    public boolean isInsertion() { return varType == INSERTION_VARIANT; }
    public boolean isDeletion() { return varType == DELETION_VARIANT ; }
    public boolean isIndel() { return isInsertion() || isDeletion() || varType == INDEL_VARIANT; }
    public boolean isReference()  { return varType == NO_VARIANT; }

    // ----------------------------------------------------------------------
    //
    // formatting
    //
    // ----------------------------------------------------------------------
    public String toString() {
        return String.format("%s\t%d\t%d\t%s\t%s\t%s",
                getContig(), getStart(), getStop(), name, refBases, observedString);
    }

    public String toSimpleString() {
        return String.format("%s:%s", name, observedString);
    }

    public String toMediumString() {
        String s = String.format("%s:%s:%s", getLocation().toString(), name, observedString);
        if ( isSNP() ) s += ": SNP";
        else {
        	if ( isInsertion() ) s += ": Insertoin";
        	else {
        		if ( isDeletion() ) s+= ": Deletion";
        		else {
        			if ( isIndel() ) s+=": Indel";
        			else s+=": Reference";
        		}
        	}
        }
        return s;        
    }

    public String repl() {
        return String.format("REPL not implemented yet");
    }

 
	@Override
	public String getAltBasesFWD() {
		if ( ! isSNP() && ! isIndel() ) return emptyStr;

		if ( isSNP() ) {
			if  ( observedAlleles.get(0).charAt(0) == refBaseChar ) return observedAlleles.get(1);
			else return observedAlleles.get(0);
		}
		
		if ( isInsertion() ) {
			if ( observedAlleles.get(0) == emptyStr ) return observedAlleles.get(1);
			else return observedAlleles.get(0);
		}
		
		if ( isDeletion() ) {
			if ( observedAlleles.get(0) == emptyStr ) return observedAlleles.get(0);
			else return observedAlleles.get(1);
			
		}
		
		System.out.printf("WARNING: unexpected variant type in pileup %s at %s%n",name,getLocation().toString());
		return null;
	}

	@Override
	public char getAltSnpFWD() throws IllegalStateException {
		if ( ! isSNP() ) throw new IllegalStateException("Variant is not a SNP");
		if  ( observedAlleles.get(0).charAt(0) == refBaseChar ) return observedAlleles.get(1).charAt(0);
		else return observedAlleles.get(0).charAt(0);

	}

	@Override
	public double getConsensusConfidence() {
		return consensusScore;
	}


	@Override
	public double getMAF() {
		if ( nNonref > 1 ) System.out.println("SAM pileup: WARNING: can not determine minor allele freq for multiallelic site");
		if ( isSNP() || isIndel() ) {
			if ( observedAlleles.get(0).equals(observedAlleles.get(1)) ) return 1.0;
			else return 0.5;
		}
		return 0;
	}

	@Override
	public int getPloidy() throws IllegalStateException {
		return 2; // ???
	}

	@Override
	public double getVariationConfidence() {
		return variantScore;
	}

	@Override
	public boolean isGenotype() {
		return true;
	}

	@Override
	public boolean isBiallelic() {
		return nNonref  < 2;
	}
}

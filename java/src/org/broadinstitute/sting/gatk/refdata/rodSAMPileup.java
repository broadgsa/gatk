package org.broadinstitute.sting.gatk.refdata;


import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.Pileup;

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
public class rodSAMPileup extends ReferenceOrderedDatum implements Genotype, Pileup {
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

    protected GenomeLoc loc;       // genomic location of this genotyped site
    // Reference sequence chromosome or scaffold
    // Start and stop positions in chrom


    protected char refBaseChar; // what we have set for the reference base (is set to a '*' for indel!)
    protected String refBases;        // the reference base sequence according to NCBI
    protected String observedString; // store the actual string representation of observed alleles

    protected String pileupQuals;     // the read base qualities
    protected String pileupBases;     // the read bases themselves

    protected List<String> observedAlleles = null;    // The sequences of the observed alleles from rs-fasta files
    protected int varType = NO_VARIANT;
    protected int ploidy = 2; // how many allelic variants we observe?
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

            parseBasesAndQuals(parts[8], parts[9]);

            observedAlleles = new ArrayList<String>(2);

            if ( refBaseChar == '*' ) {
                parseIndels(parts[3]) ;
                if ( varType == DELETION_VARIANT ) loc = new GenomeLoc(contig, start, start+eventLength-1);
                else loc = new GenomeLoc(contig, start, start-1); // if it's not a deletion and we are biallelic, this got to be an insertion; otherwise the state is inconsistent!!!!
            }
            else {
                // if the variant is a SNP or a reference base (i.e. no variant at all)
                assert parts[3].length() == 1 : "non-indel genotype is expected to be represented by a single letter";
                refBases = parts[2].toUpperCase();
                eventLength = 1;
                //loc = new GenomeLoc(contig, start, start+1);
                loc = new GenomeLoc(contig, start, start);

                char ch = parts[3].charAt(0);

                switch ( ch ) {
                    case 'A': observedAlleles.add(baseA); observedAlleles.add(baseA); break;
                    case 'C': observedAlleles.add(baseC); observedAlleles.add(baseC); break;
                    case 'G': observedAlleles.add(baseG); observedAlleles.add(baseG); break;
                    case 'T': observedAlleles.add(baseT); observedAlleles.add(baseT); break;
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
            System.out.printf("  Exception caught during parsing BasicPileup line:  %s%n", Utils.join(" <=> ", parts));
            throw e;
        }
  //      if ( nNonref > 1 ) System.out.println("SAM pileup: WARNING: multi-allelic variant :  ("+refBaseChar+") -->"+toMediumString());
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

            String varBases = obs[i].toUpperCase();

            switch ( obs[i].charAt(0) )  {
                case '+':
                    if ( varType != NO_VARIANT && varType != INSERTION_VARIANT )  varType = INDEL_VARIANT;
                    else varType = INSERTION_VARIANT;
                    refBases = emptyStr;
                    break;
                case '-' :
                    if ( varType != NO_VARIANT && varType != DELETION_VARIANT )  varType = INDEL_VARIANT;
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
            // or both deletions as compared to +ACC/-ACC which would still have the same bases 
            if ( observedAlleles.get(0).equals(observedAlleles.get(1))  ) nNonref = 1;
            else nNonref = 2;
        }
        // DONE with indels

    }

    private void parseBasesAndQuals(final String bases, final String quals)
    {
        //System.out.printf("%s%n%s%n", bases, quals);

        // needs to convert the base string with it's . and , to the ref base
        StringBuilder baseBuilder = new StringBuilder();
        StringBuilder qualBuilder = new StringBuilder();
        boolean done = false;
        for ( int i = 0, j = 0; i < bases.length() && ! done; i++ ) {
            //System.out.printf("%d %d%n", i, j);
            char c = (char)bases.charAt(i);

            switch ( c ) {
                case '.':   // matches reference
                case ',':   // matches reference
                    baseBuilder.append(refBaseChar);
                    qualBuilder.append((char)quals.charAt(j++));
                    break;
                case '$':   // end of read
                    break;
                case '*':   // end of indel?
                    j++;
                    break;
                case '^':   // mapping quality
                    i++;
                    break;
                case '+':   // start of indel
                case '-':   // start of indel
                    final Pattern regex = Pattern.compile("([0-9]+).*");             // matches case 1
                    final String rest = bases.substring(i+1);
                    //System.out.printf("sub is %s%n", rest);
                    Matcher match = regex.matcher(rest);
                    if ( ! match.matches() ) {
                        if ( refBaseChar != '*' )
                            throw new RuntimeException("Bad pileup format: " + bases + " at position " + i);
                        done = true;
                    }
                    else {
                        String g = match.group(1);
                        //System.out.printf("group is %d, match is %s%n", match.groupCount(), g);
                        int l = Integer.parseInt(g);
                        i += l + g.length();    // length of number + that many bases + +/- at the start (included in the next i++)
                        //System.out.printf("remaining is %d => %s%n", l, bases.substring(i+1));
                    }
                    break;
                default:   // non reference base
                    baseBuilder.append(c);
                    qualBuilder.append((char)quals.charAt(j++));
            }
        }
        
        pileupBases = baseBuilder.toString();
        pileupQuals = qualBuilder.toString();
    }

    public GenomeLoc getLocation()  { return loc; }
    public String getQuals()        { return pileupQuals; }
    public char getRef()            { return refBaseChar; }
    public int size()               { return pileupQuals.length(); }
    public String getBases()        { return pileupBases; }

    public String getPileupString()
    {
        return String.format("%s: %s %s %s", getLocation(), getRef(), getBases(), getQuals());
    }


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
  /*  public char getRefSnpFWD() throws IllegalStateException {
        if ( isIndel() ) throw new IllegalStateException("Variant is not a SNP");
        return refBaseChar;
    }
*/
    
    @Override
    public List<String> getFWDAlleles()  {
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

    public boolean isHom() {
    	// implementation-dependent: here we use the fact that for ref and snps we actually use fixed static strings to remember the genotype
    	if ( ! isIndel() ) return ( observedAlleles.get(0) == observedAlleles.get(1) ); 
    	return ( isInsertion() || isDeletion() ) && observedAlleles.get(0).equals(observedAlleles.get(1) );
    }
    
    public boolean isHet() {
    	// implementation-dependent: here we use the fact that for ref and snps we actually use fixed static strings to remember the genotype
    	if ( ! isIndel() ) return ( observedAlleles.get(0) != observedAlleles.get(1) ); 
    	return isIndel() || ( ! observedAlleles.get(0).equals(observedAlleles.get(1) ) );
    }

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

        String s = null;
        if ( refBaseChar == '*' ) s = String.format("%s:%s:%s", getLocation().toString(), name,observedString);
        else s = String.format("%s:%s:%s/%s", getLocation().toString(), name, observedAlleles.get(0),observedAlleles.get(1));
        if ( isSNP() ) s += ": SNP";
        else {
            if ( isInsertion() ) s += ": Insertion";
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

/*
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
*/
    
   /*
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
*/

 /*   
    @Override
    public int getPloidy() throws IllegalStateException {
        return 2; // ???
    }
*/
    
    @Override
    public double getVariantConfidence() {
        return variantScore;
    }


    @Override
    public boolean isBiallelic() {
        return nNonref  < 2;
    }

	@Override
	public double getConsensusConfidence() {
		return consensusScore;
	}

	@Override
	public String getFWDRefBases() {
		return refBases;
	}
}

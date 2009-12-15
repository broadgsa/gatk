package org.broadinstitute.sting.gatk.refdata;


import java.io.FileNotFoundException;
import java.util.*;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import org.broadinstitute.sting.utils.*;

import net.sf.picard.reference.ReferenceSequenceFileWalker;

/**
 *  This class wraps Maq/samtools allele calls from pileup format and presents them as a ROD.<br>
 *
 * Example format:<br>
 *   for SNP:<br>
 *      [chr]       [pos] [ref]  [consensus allele(s)]    [consensus confidence]  [snp confidence]   [max mapping qual]  [num reads in the pile] <br>
 *     chrX       466     T                     Y                                   170                                 170                          88                              32 ... (piles of read bases  and quals follow) <br>
 *     <br>
 *   for indel: <br>
 *    [chr]         [pos]     [always *]  [consensus alleles]    [consensus conf.?]        [indel conf.?]     [max mapping qual]     [num reads in the pile]  [indel]   [always *?] [reads haveindel]  [reads may have indel] [other reads?]<br>     
 *     chrX       141444        *           +CA/+CA                               32                                       468            255                                         25                                        +CA     *                         5                                   2                                   12             <br>
 * User: asivache
 * Date: Apr 03, 2009
 * Time: 2:58:33 PM
 * To change this template use File | Settings | File Templates.
 */


public class SAMPileupRecord implements Genotype, GenotypeList {
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
    protected String refBases;        // the reference base sequence according to NCBI; single base for point mutations, deleted bases for  deletions, empty string for insertions
    protected String observedString; // stores the actual string representation of observed alleles (for point mutations stores as A/C, not in extended alphabet!!)

    protected String pileupQuals;     // the read base qualities
    protected String pileupBases;     // the read bases themselves

    protected List<String> observedAlleles = null;    // The sequences of the observed alleles (e.g. {"A","C"} for point mutation or {"","+CC"} for het. insertion
    protected int varType = NO_VARIANT;
    protected int ploidy = 2; // how many allelic variants we have?
    protected int nNonref = 0; // number of non-reference alleles observed
    protected int eventLength = 0; // number of inserted or deleted bases

    protected double consensusScore = 0;
    protected double variantScore = 0;

    private char fldDelim = '\t'; 
    private String name = emptyStr;
    // ----------------------------------------------------------------------
    //
    // Constructors
    //
    // ----------------------------------------------------------------------
   public SAMPileupRecord(String name) {
    	this.name = name; // keeps track name
   }

   public boolean parseLine(final String line) {
//       0          1             2         3                  4         5            6         7
//*     chrX     466           T         Y                170      170       88      32 ... (piles of read bases  and quals follow)
//*    chrX    141444     *     +CA/+CA       32       468     255     25     +CA     *       5       2       12      6
	   try {
	
		   int[] pos = Utils.indexOfAll(line, fldDelim );
		   String contig = line.substring(0,pos[0] ); // field 0
		   long start = Long.parseLong(line.substring(pos[0]+1, pos[1])) ; // field 1 

		   refBaseChar = Character.toUpperCase(line.charAt(pos[1]+1)); // field 2 ( single char)

		   observedString = line.substring(pos[2]+1, pos[3]).toUpperCase(); // field 3
		   observedAlleles = new ArrayList<String>(2);

		   consensusScore = Double.parseDouble(line.substring(pos[3]+1,pos[4]));
		   variantScore = Double.parseDouble(line.substring(pos[4]+1, pos[5]));

		   if ( refBaseChar == '*' ) {
 	
			   parseIndels(observedString) ;
			   if ( varType == DELETION_VARIANT ) loc = GenomeLocParser.createGenomeLoc(contig, start, start+eventLength-1);
			   else loc = GenomeLocParser.createGenomeLoc(contig, start, start-1); // if it's not a deletion and we are biallelic, this got to be an insertion; otherwise the state is inconsistent!!!!
		   } else {
			   parseBasesAndQuals(line,pos[7]+1,pos[8], pos[8]+1, ( pos.length > 9 ? pos[9] : line.length()) );
//			   parseBasesAndQuals(line.substring(pos[7]+1,pos[8]), line.substring(pos[8]+1, ( pos.length > 9 ? pos[9] : line.length()) ) );
			   // if the variant is a SNP or a reference base (i.e. no variant at all)
			   if ( observedString.length() != 1 ) throw new RuntimeException( "point mutation genotype is expected to be represented by a single letter");
     
			   refBases = line.substring(pos[1]+1, pos[2]).toUpperCase();
			   eventLength = 1;
			   //loc = new GenomeLoc(contig, start, start+1);
			   loc = GenomeLocParser.createGenomeLoc(contig, start, start);

			   char ch = observedString.charAt(0);

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
				   // 	we know that at least one allele is non-ref;
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
		   System.out.printf("  Exception caught during parsing BasicPileup line: %s%n", line);
		   throw e;
	   }
	   return true;
	   //      	if ( nNonref > 1 ) System.out.println("SAM pileup: WARNING: multi-allelic variant :  ("+refBaseChar+") -->"+toMediumString());
   }

   
    /** Parses a line from SAM pileup file into this SAMPileupRecord object.. 
     * @param parts line from SAM pileup file, split on delimiter character (tab)
     */
    public boolean parseLine(final String[] parts) {
//    	          0          1             2         3                  4         5            6         7
//    	*     chrX     466           T         Y                170      170       88      32 ... (piles of read bases  and quals follow)
//    	 *    chrX    141444     *     +CA/+CA       32       468     255     25     +CA     *       5       2       12      6
        try {
        	
            String contig = parts[0];
            long start = Long.parseLong(parts[1]) ;

            refBaseChar = Character.toUpperCase(parts[2].charAt(0));

            parts[3] = parts[3].toUpperCase();
            observedString = parts[3];
            observedAlleles = new ArrayList<String>(2);

        	consensusScore = Double.parseDouble(parts[4]);
            variantScore = Double.parseDouble(parts[5]);

            if ( refBaseChar == '*' ) {
            	
                parseIndels(parts[3]) ;
                if ( varType == DELETION_VARIANT ) loc = GenomeLocParser.parseGenomeLoc(contig, start, start+eventLength-1);
                else loc = GenomeLocParser.parseGenomeLoc(contig, start, start-1); // if it's not a deletion and we are biallelic, this got to be an insertion; otherwise the state is inconsistent!!!!
            }
            else {
                parseBasesAndQuals(parts[8], parts[9]);
                // if the variant is a SNP or a reference base (i.e. no variant at all)
                if ( parts[3].length() != 1 ) throw new RuntimeException( "point mutation genotype is expected to be represented by a single letter");
                
                refBases = parts[2].toUpperCase();
                eventLength = 1;
                //loc = GenomeLoc.parseGenomeLoc(contig, start, start+1);
                loc = GenomeLocParser.parseGenomeLoc(contig, start, start);

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
            System.out.printf("  Exception caught during parsing BasicPileup line: %s%n", Utils.join(" <=> ", parts));
            throw e;
        }
        return true;
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
                observedAlleles.add(emptyStr);
                continue; 
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
                nNonref = 0; // no observations of non-reference allele at all
                refBases = emptyStr;
            } else {
                nNonref = 1; // hasRefAllele = true, so one allele was definitely ref, hence there is only one left
            }
        } else {
            // we observe two non-ref alleles; they better be the same variant, otherwise the site is not bi-allelic and at the moment we
            // fail to set data in a consistent way.
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

    private void parseBasesAndQuals(final String line, int basesStart, int basesStop, int qualsStart, int qualsStop)
    {
        //System.out.printf("%s%n%s%n", bases, quals);

        // needs to convert the base string with it's . and , to the ref base
        StringBuilder baseBuilder = new StringBuilder();
        StringBuilder qualBuilder = new StringBuilder();
        boolean done = false;
        for ( int i = basesStart, j = qualsStart; i < basesStop && ! done; i++ ) {
            //System.out.printf("%d %d%n", i, j);
            char c = (char)line.charAt(i);

            switch ( c ) {
                case '.':   // matches reference
                case ',':   // matches reference
                    baseBuilder.append(refBaseChar);
                    qualBuilder.append((char)line.charAt(j++));
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
                    final String rest = line.substring(i+1,basesStop);
                    //System.out.printf("sub is %s%n", rest);
                    Matcher match = regex.matcher(rest);
                    if ( ! match.matches() ) {
                        if ( refBaseChar != '*' )
                            throw new RuntimeException("Bad pileup format: " + line.substring(basesStart, basesStop) + " at position " + (i-basesStart));
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
                    qualBuilder.append((char)line.charAt(j++));
            }
        }
        
        pileupBases = baseBuilder.toString();
        pileupQuals = qualBuilder.toString();
    }

    
    public GenomeLoc getLocation()  { return loc; }
    public String getQualsAsString()        { return pileupQuals; }
    
    /** Returns reference base for point genotypes or '*' for indel genotypes, as a char.
     * 
     */
    public char getRef()            { return refBaseChar; }
    public int size()               { return pileupQuals.length(); }
    
    /** Returns pile of observed bases over the current genomic location.
     *  
     */
    public String getBasesAsString()        { return pileupBases; }

    /** Returns formatted pileup string for the current genomic location as
     * "location: reference_base observed_base_pile observed_qual_pile"
     */
    public String getPileupString()
    {
        return String.format("%s: %s %s %s", getLocation(), getRef(), getBasesAsString(), getQualsAsString());
    }

    public ArrayList<Byte> getBasesAsArrayList()        { throw new StingException("Not implemented"); }
    public ArrayList<Byte> getQualsAsArrayList()        { throw new StingException("Not implemented"); }
    public byte[] getBases()        { throw new StingException("Not implemented"); }
    public byte[] getQuals()        { throw new StingException("Not implemented"); }

    /** Returns bases in the reference allele as a String. For point genotypes, the string consists of a single
     * character (reference base). For indel genotypes, the string is empty for insertions into
     * the reference, or consists of deleted bases for deletions.
     *
     * @return reference allele, forward strand
     */
    @Override
    public String getFWDRefBases() {
        return refBases;
    }

        
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
                getLocation().getContig(), getLocation().getStart(), getLocation().getStop(), name, refBases, observedString);
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
    public int length() {
        return eventLength;
    }

	@Override
	public boolean isIndelGenotype() {
		return refBaseChar == '*';
	}


	@Override
	public boolean isPointGenotype() {
		return ! isIndelGenotype();
	}

	/** Implements method required by GenotypeList interface. If this object represents
	 * an indel genotype, then it returns itself through this method. If this object is a
	 * point genotype, this method returns null.
	 * @return
	 */
	@Override
	public Genotype getIndelGenotype() {
		if ( isIndelGenotype() ) return this;
		else return null;
	}

	/** Implements method required by GenotypeList interface. If this object represents
	 * a point genotype, then it returns itself through this method. If this object is an
	 * indel genotype, this method returns null.
	 * @return
	 */
	@Override
	public Genotype getPointGenotype() {
		if ( isPointGenotype() ) return this;
		else return null;
	}

	/** Returns true if this object \em is an indel genotype (and thus
	 * indel genotype is what it only has).
	 * @return
	 */
	@Override
	public boolean hasIndelGenotype() {
		return isIndelGenotype();
	}

	/** Returns true if this object \em is a point genotype (and thus
	 * point genotype is what it only has.
	 * @return
	 */
	@Override
	public boolean hasPointGenotype() {
		return isPointGenotype();
	}
	

	@Override
	public int compareTo(ReferenceOrderedDatum o) {
		return getLocation().compareTo(o.getLocation());
	}

	protected static class pileupFileIterator implements Iterator<SAMPileupRecord> {
//		private String fieldDelimiter = new String("\t");
		private String rodName = null;
		private xReadLines parser = null;
		private String lastProcessedLine = null;
		
//		private long z = 0;
//		private long t = 0;
		
		pileupFileIterator(String name, java.io.File f) {
			try {
				parser = new xReadLines(f);
			} catch ( FileNotFoundException e ) {
				Utils.scareUser("Couldn't open file: " + f);
			}
			rodName = name;
		}
		
		@Override
		public boolean hasNext() {
			return parser.hasNext();
		}

        /**
         * @return the next element in the iteration.
         * @throws NoSuchElementException - iterator has no more elements.
         */
		@Override
		public SAMPileupRecord next() {
             if (!this.hasNext()) throw new NoSuchElementException("SAMPileupRecord next called on iterator with no more elements");
//			 if ( z == 0 ) t = System.currentTimeMillis();
	          lastProcessedLine = parser.next();
	          SAMPileupRecord n = new SAMPileupRecord(rodName);
	           //System.out.printf("Line is %s%n", line);
//	          String parts[] = lastProcessedLine.split(fieldDelimiter);
//	          n.parseLine(parts);

	          n.parseLine(lastProcessedLine);
//	          if ( ++z == 1000000 ) { System.out.println(rodName + ": 1,000,000 records read in "+((double)(System.currentTimeMillis()-t)/1000.0) +" s");  z = 0; }
	          return n ;
		}


		@Override
		public void remove() {
			throw new UnsupportedOperationException("'remove' operation is not supported for file-backed SAM pileups");
		}
		
	}
	
	
	public static Iterator<SAMPileupRecord> createIterator(String name, java.io.File f) {
		return new pileupFileIterator(name,f);
	}
	

    
	public static void main(String argv[]) {
//		String testFile = "/humgen/gsa-scr1/asivache/TCGA/Ovarian/C2K/0805/normal.pileup";
		String testFile = "/humgen/gsa-scr1/asivache/trios/CEU/NA12891.12892.12878/mother.chr1.pileup.indel";
		
		Iterator<SAMPileupRecord> it = createIterator("test-normal", new java.io.File(testFile));
		
		ReferenceSequenceFileWalker reference = new ReferenceSequenceFileWalker(
              //new java.io.File(  "/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta")
              new java.io.File(  "/humgen/gsa-scr1/asivache/trios/CEU/NA12891.12892.12878/human_b36_both.fasta")
                     );

		if ( reference.getSequenceDictionary() == null ) {
			System.out.println("No reference sequence dictionary found. Abort.");
			System.exit(1);
		}

		GenomeLocParser.setupRefContigOrdering(reference.getSequenceDictionary());
		
		int counter = 0;
		
		while ( it.hasNext() && counter < 430 ) {
			SAMPileupRecord p = it.next();
			
			System.out.print(p.getLocation().toString());
			System.out.print('\t');
			
			if ( p.isIndel() && p.isSNP() ) { System.out.print("Indel+SNP"); }
			else {
				if ( p.isSNP() ) { System.out.print("SNP"); }
				else { 
					if ( p.isIndel() ) { System.out.print("Indel"); }
					else { System.out.print("REF"); }
				}
			}
			
			System.out.print('\t');
			System.out.print(p.getFWDAlleles().get(0)+"/"+p.getFWDAlleles().get(1));
			System.out.print('\t');
			System.out.println(p.getConsensusConfidence()+"\t"+p.getVariantConfidence());
			counter++;
		}
	}




}

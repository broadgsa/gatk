package org.broadinstitute.sting.utils;

import edu.mit.broad.picard.directed.IntervalList;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.util.Interval;
import net.sf.functionalj.Function1;
import net.sf.functionalj.FunctionN;
import net.sf.functionalj.Functions;
import net.sf.functionalj.reflect.JdkStdReflect;
import net.sf.functionalj.reflect.StdReflect;
import net.sf.functionalj.util.Operators;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import org.apache.log4j.Logger;

import java.io.File;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by IntelliJ IDEA.
 * User: mdepristo
 * Date: Mar 2, 2009
 * Time: 8:50:11 AM
 *
 * Genome location representation.  It is *** 1 *** based
 *
 *
 */
public class GenomeLoc implements Comparable<GenomeLoc>, Cloneable {
    private static Logger logger = Logger.getLogger(GenomeLoc.class);

    private Integer contigIndex;
    private long start;
    private long stop;

    // --------------------------------------------------------------------------------------------------------------
    //
    // Ugly global variable defining the optional ordering of contig elements
    //
    // --------------------------------------------------------------------------------------------------------------
    //public static Map<String, Integer> refContigOrdering = null;
    private static SAMSequenceDictionary contigInfo = null;

    public static boolean hasKnownContigOrdering() {
        return contigInfo != null;
    }

    
    public static SAMSequenceRecord getContigInfo( final String contig ) {
        return contigInfo.getSequence(contig);
    }

    /**
     * Returns the contig index of a specified string version of the contig
     * @param contig the contig string
     * @return the contig index, -1 if not found
     */
    public static int getContigIndex( final String contig ) {
        if (contigInfo.getSequenceIndex(contig) == -1)
            Utils.scareUser(String.format("Contig %s given as location, but this contig isn't present in the Fasta sequence dictionary", contig));

        return contigInfo.getSequenceIndex(contig);
    }

    public static boolean setupRefContigOrdering(final ReferenceSequenceFile refFile) {
        return setupRefContigOrdering(refFile.getSequenceDictionary());
    }

    public static boolean setupRefContigOrdering(final SAMSequenceDictionary seqDict) {
        if (seqDict == null)  { // we couldn't load the reference dictionary
        	logger.info("Failed to load reference dictionary, falling back to lexicographic order for contigs");
            Utils.scareUser("Failed to load reference dictionary");
            return false;
        } else if ( contigInfo == null ){
            contigInfo = seqDict;
            logger.debug(String.format("Prepared reference sequence contig dictionary"));
            for (SAMSequenceRecord contig : seqDict.getSequences() ) {
                logger.debug(String.format(" %s (%d bp)", contig.getSequenceName(), contig.getSequenceLength()));
            }
        }
        
        return true;
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // constructors
    //
    // --------------------------------------------------------------------------------------------------------------
    public GenomeLoc( int contigIndex, final long start, final long stop ) {
        if(contigInfo == null) {  throw new StingException("Contig info has not been setup in the GenomeLoc context yet."); }

        if (contigIndex < 0 || contigIndex >= contigInfo.size()) {
            throw new StingException("Contig info has not been setup in the GenomeLoc context yet.");
        }
        if (start < 0) { throw new StingException("Bad start position " + start);}
        if (stop  < -1) { throw new StingException("Bad stop position " + stop); }    // a negative -1 indicates it's not a meaningful end position

        this.contigIndex = contigIndex;
        this.start = start;
        this.stop = stop == -1 ? start : stop;
    }

    public GenomeLoc(final SAMRecord read) {
        this(read.getReferenceIndex(), read.getAlignmentStart(), read.getAlignmentEnd());
    }

    public GenomeLoc( final String contig, final long start, final long stop ) {
        this(contigInfo.getSequenceIndex(contig), start, stop);
    }

    public GenomeLoc( final String contig, final long pos ) {
        this(contig, pos, pos);
    }

    public GenomeLoc( final int contig, final long pos ) {
        this(contig, pos, pos );
    }

    public GenomeLoc( final GenomeLoc toCopy ) {
        this( toCopy.contigIndex, toCopy.getStart(), toCopy.getStop() );
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // Parsing string representations
    //
    // --------------------------------------------------------------------------------------------------------------
    private static long parsePosition( final String pos ) {
        String x = pos.replaceAll(",", "");
        return Long.parseLong(x);
    }

    public static GenomeLoc parseGenomeLoc( final String str ) {
        // 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
        //System.out.printf("Parsing location '%s'%n", str);

        final Pattern regex1 = Pattern.compile("([\\w&&[^:]]+)$");                      // matches case 1
        final Pattern regex2 = Pattern.compile("([\\w&&[^:]]+):([\\d,]+)$");            // matches case 2
        final Pattern regex3 = Pattern.compile("([\\w&&[^:]]+):([\\d,]+)-([\\d,]+)$");  // matches case 3
        final Pattern regex4 = Pattern.compile("([\\w&&[^:]]+):([\\d,]+)\\+");          // matches case 4

        String contig = null;
        long start = 1;
        long stop = Integer.MAX_VALUE;
        boolean bad = false;

        Matcher match1 = regex1.matcher(str);
        Matcher match2 = regex2.matcher(str);
        Matcher match3 = regex3.matcher(str);
        Matcher match4 = regex4.matcher(str);

        try {
            if ( match1.matches() ) {
                contig = match1.group(1);
            }
            else if ( match2.matches() ) {
                contig = match2.group(1);
                start = parsePosition(match2.group(2));
                stop = start;                
            }
            else if ( match4.matches() ) {
                contig = match4.group(1);
                start = parsePosition(match4.group(2));
            }
            else if ( match3.matches() ) {
                contig = match3.group(1);
                start = parsePosition(match3.group(2));
                stop = parsePosition(match3.group(3));

                if ( start > stop )
                    bad = true;
            }
            else {
                bad = true;
            }
        } catch ( Exception e ) {
            bad = true;
        }

        if ( bad ) {
            throw new StingException("Invalid Genome Location string: " + str);
        }

        if ( stop == Integer.MAX_VALUE && hasKnownContigOrdering() ) {
            // lookup the actually stop position!
            stop = getContigInfo(contig).getSequenceLength();
        }

        GenomeLoc loc = new GenomeLoc(contig, start, stop);
 //       System.out.printf("  => Parsed location '%s' into %s%n", str, loc);

        return loc;
    }

    /**
     * Useful utility function that parses a location string into a coordinate-order sorted
     * array of GenomeLoc objects
     *
     * @param str String representation of genome locs.  Null string corresponds to no filter.
     * @return Array of GenomeLoc objects corresponding to the locations in the string, sorted by coordinate order
     */
    public static ArrayList<GenomeLoc> parseGenomeLocs(final String str) {
        // Null string means no filter.
        if( str == null ) return null;

        // Of the form: loc1;loc2;...
        // Where each locN can be:
        // 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
        StdReflect reflect = new JdkStdReflect();
        FunctionN<GenomeLoc> parseOne = reflect.staticFunction(GenomeLoc.class, "parseGenomeLoc", String.class);
        Function1<GenomeLoc, String> f1 = parseOne.f1();
        try {
            Collection<GenomeLoc> result = Functions.map(f1, Arrays.asList(str.split(";")));
            ArrayList<GenomeLoc> locs = new ArrayList(result);
            Collections.sort(locs);
            //logger.info(String.format("Going to process %d locations", locs.length));
            locs = mergeOverlappingLocations(locs);
            logger.info("Locations are:" + Utils.join(", ", Functions.map(Operators.toString, locs)));
            return locs;
        } catch (Exception e) {
            e.printStackTrace();
            Utils.scareUser(String.format("Invalid locations string: %s, format is loc1;loc2; where each locN can be 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'", str));
            return null;
        }
    }

    public static ArrayList<GenomeLoc> mergeOverlappingLocations(final ArrayList<GenomeLoc> raw) {
        logger.debug("  Raw locations are:\n" + Utils.join("\n", Functions.map(Operators.toString, raw)));        
        if ( raw.size() <= 1 )
            return raw;
        else {
            ArrayList<GenomeLoc> merged = new ArrayList<GenomeLoc>();
            Iterator<GenomeLoc> it = raw.iterator();
            GenomeLoc prev = it.next();
            while ( it.hasNext() ) {
                GenomeLoc curr = it.next();
                if ( prev.contiguousP(curr) ) {
                    prev = prev.merge(curr);
                } else {
                    merged.add(prev);
                    prev = curr;
                }
            }
            merged.add(prev);
            return merged;
        }
    }


    /**
     * Move this Genome loc to the next contig, with a start
     * and stop of 1.
     * @return true if we are not out of contigs, otherwise false if we're
     *          at the end of the genome (no more contigs to jump to).
     */
    public boolean toNextContig() {
        if ((contigIndex + 1) < GenomeLoc.contigInfo.size()) {
            this.contigIndex++;
            this.start = 1;
            this.stop = 1;
            return true;
        }
        return false;
    }


    /**
     * Returns true iff we have a specified series of locations to process AND we are past the last
     * location in the list.  It means that, in a serial processing of the genome, that we are done.
     *
     * @param curr Current genome Location
     * @return true if we are past the last location to process
     */
    public static boolean pastFinalLocation(GenomeLoc curr, List<GenomeLoc> locs) {
        return (locs.size() > 0 && curr.isPast(locs.get(locs.size() - 1)));
    }

    /**
     * A key function that returns true if the proposed GenomeLoc curr is within the list of
     * locations we are processing in this TraversalEngine
     *
     * @param curr
     * @return true if we should process GenomeLoc curr, otherwise false
     */
    public static boolean inLocations(GenomeLoc curr, ArrayList<GenomeLoc> locs) {
        if ( locs.size() == 0 ) {
            return true;
        } else {
            for ( GenomeLoc loc : locs ) {
                //System.out.printf("  Overlap %s vs. %s => %b%n", loc, curr, loc.overlapsP(curr));
                if (loc.overlapsP(curr))
                    return true;
            }
            return false;
        }
    }

    public static void removePastLocs(GenomeLoc curr, List<GenomeLoc> locs) {
        while ( !locs.isEmpty() && curr.isPast(locs.get(0)) ) {
            //System.out.println("At: " + curr + ", removing: " + locs.get(0));
            locs.remove(0);
        }
    }

    public static boolean overlapswithSortedLocsP(GenomeLoc curr, List<GenomeLoc> locs, boolean returnTrueIfEmpty) {
        if ( locs.isEmpty() )
            return returnTrueIfEmpty;

        // skip loci before intervals begin
        if ( hasKnownContigOrdering() && curr.contigIndex < locs.get(0).contigIndex )
            return false;

        for ( GenomeLoc loc : locs ) {
            //System.out.printf("  Overlap %s vs. %s => %b%n", loc, curr, loc.overlapsP(curr));
            if ( loc.overlapsP(curr) )
                return true;
            if ( curr.compareTo(loc) < 0 )
                return false;
        }
        return false;
    }

    //
    // Accessors and setters
    //
    public final String getContig() {
        //this.contigIndex != -1;
        if (!(contigInfo != null && contigInfo.getSequences() != null)) {
            throw new StingException("The contig information or it's sequences are null");
        }
        if ((this.contigIndex < 0) || (this.contigIndex >= contigInfo.getSequences().size())) {
            throw new StingException("The contig index is not bounded by the zero and seqeunce count, contig index: " + contigIndex);
        }
        if (contigInfo.getSequence(this.contigIndex) == null ||
            contigInfo.getSequence(this.contigIndex).getSequenceName() == null) {
            throw new StingException("The associated sequence index for contig " + contigIndex + " is null");
        }
        return contigInfo.getSequence(this.contigIndex).getSequenceName();
        //if (contigInfo != null && contigInfo.getSequence(this.contigIndex) != null) {
        //    return contigInfo.getSequence(this.contigIndex).getSequenceName();
        //}

        //return null;
    }

    public final int getContigIndex() { return this.contigIndex; }
    public final long getStart()    { return this.start; }
    public final long getStop()     { return this.stop; }
    public final String toString()  {
        if ( throughEndOfContigP() && atBeginningOfContigP() )
            return getContig();
        else if ( throughEndOfContigP() || getStart() == getStop() )
            return String.format("%s:%d", getContig(), getStart());
        else
            return String.format("%s:%d-%d", getContig(), getStart(), getStop());
    }

    public final boolean isUnmapped() { return this.contigIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX; }
    public final boolean throughEndOfContigP() { return this.stop == Integer.MAX_VALUE; }
    public final boolean atBeginningOfContigP() { return this.start == 1; }

    public void setContig(String contig) {
        this.contigIndex = contigInfo.getSequenceIndex(contig);
    }

    public void setStart(long start) {
        this.start = start;
    }
    public void setStop(long stop) {
        this.stop = stop;
    }

    public final boolean isSingleBP() { return stop == start; }

    public final boolean disjointP(GenomeLoc that) {
        if ( this.contigIndex != that.contigIndex ) return true;    // different chromosomes
        if ( this.start > that.stop ) return true;                  // this guy is past that
        if ( that.start > this.stop ) return true;                  // that guy is past our start
        return false;
    }

    public final boolean discontinuousP(GenomeLoc that) {
        if ( this.contigIndex != that.contigIndex ) return true;    // different chromosomes
        if ( (this.start - 1) > that.stop ) return true;            // this guy is past that
        if ( (that.start - 1) > this.stop ) return true;            // that guy is past our start
        return false;
    }

    public final boolean overlapsP(GenomeLoc that) {
        return ! disjointP( that );
    }

    public final boolean contiguousP(GenomeLoc that) {
        return ! discontinuousP( that );
    }

    public GenomeLoc merge( GenomeLoc that ) throws StingException {
        if (!(this.contiguousP(that))) {
            throw new StingException("The two genome loc's need to be contigous");
        }

        return new GenomeLoc(getContig(),
                             Math.min(getStart(), that.getStart()),
                             Math.max( getStop(), that.getStop()) );
    }

    public final boolean containsP(GenomeLoc that) {
        if ( ! onSameContig(that) ) return false;
        return getStart() <= that.getStart() && getStop() >= that.getStop();
    }

    public final boolean onSameContig(GenomeLoc that) {
        return (this.contigIndex == that.contigIndex);
    }

    public final int minus( final GenomeLoc that ) {
        if ( this.contigIndex == that.contigIndex )
            return (int) (this.getStart() - that.getStart());
        else
            return Integer.MAX_VALUE;
    }

    public final int distance( final GenomeLoc that ) {
        return Math.abs(minus(that));
    }    

    public final boolean isBetween( final GenomeLoc left, final GenomeLoc right ) {
        return this.compareTo(left) > -1 && this.compareTo(right) < 1;
    }

    public final boolean isBefore( GenomeLoc that ) {
        int comparison = this.compareContigs(that);
        return ( comparison == -1 || ( comparison == 0 && this.getStop() < that.getStart() ));        
    }

    public final boolean isPast( GenomeLoc that ) {
        int comparison = this.compareContigs(that);
        return ( comparison == 1 || ( comparison == 0 && this.getStart() > that.getStop() ));
    }

    public final void incPos() {
        incPos(1);
    }
    public final void incPos(long by) {
        this.start += by;
        this.stop += by;
    }

    public final GenomeLoc nextLoc() {
        GenomeLoc n = new GenomeLoc(this);
        n.incPos();
        return n;
    }

    /**
     * Check to see whether two genomeLocs are equal.
     * Note that this implementation ignores the contigInfo object.
     * @param other Other contig to compare.
     */
    @Override
    public boolean equals(Object other) {
        if(other == null)
            return false;
        if(other instanceof GenomeLoc) {
            GenomeLoc otherGenomeLoc = (GenomeLoc)other;
            return this.contigIndex.equals(otherGenomeLoc.contigIndex) &&
                   this.start == otherGenomeLoc.start &&
                   this.stop == otherGenomeLoc.stop;
        }
        return false;
    }

    /**
     * Return a new GenomeLoc at this same position.
     * @return A GenomeLoc with the same contents as the current loc.
     */
    @Override
    public GenomeLoc clone() {
        return new GenomeLoc(this);
    }

    //
    // Comparison operations
    //
    // TODO: get rid of this method because it's sloooooooooooooow
    @Deprecated
    public static int compareContigs( final String thisContig, final String thatContig )
    {
        if ( thisContig == thatContig )
        {
            // Optimization.  If the pointers are equal, then the contigs are equal.
            return 0;
        }

        if ( hasKnownContigOrdering() )
        {
            int thisIndex = getContigIndex(thisContig);
            int thatIndex = getContigIndex(thatContig);

            if ( thisIndex == -1 )
            {
                if ( thatIndex == -1 )
                {
                    // Use regular sorted order
                    return thisContig.compareTo(thatContig);
                }
                else
                {
                    // this is always bigger if that is in the key set
                    return 1;
                }
            }
            else if ( thatIndex == -1 )
            {
                return -1;
            }
            else
            {
                if ( thisIndex < thatIndex ) return -1;
                if ( thisIndex > thatIndex ) return 1;
                return 0;
            }
        }
        else
        {
            return thisContig.compareTo(thatContig);
        }
    }

    public final int compareContigs( GenomeLoc that ) {
        return (this.contigIndex == that.contigIndex)?0:((this.contigIndex < that.contigIndex)?-1:1);
    }

    public int compareTo( GenomeLoc that ) {
        if ( this == that ) return 0;

        final int cmpContig = compareContigs(that);

        if ( cmpContig != 0 ) return cmpContig;
        if ( this.getStart() < that.getStart() ) return -1;
        if ( this.getStart() > that.getStart() ) return 1;

        // TODO: and error is being thrown because we are treating reads with the same start positions
        // but different stop as out of order
        //if ( this.getStop() < that.getStop() ) return -1;
        //if ( this.getStop() > that.getStop() ) return 1;
        return 0;
    }


    /**
     * Read a file of genome locations to process.
     * regions specified by the location string.  The string is of the form:
     * Of the form: loc1;loc2;...
     * Where each locN can be:
     * 'chr2', 'chr2:1000000' or 'chr2:1,000,000-2,000,000'
     *
     * @param file_name
     */
    public static ArrayList<GenomeLoc> IntervalFileToList(final String file_name) {
// first try to read it as an interval file since that's well structured
        // we'll fail quickly if it's not a valid file.  Then try to parse it as
        // a location string file
        ArrayList<GenomeLoc> ret = null;
        try {
            IntervalList il = IntervalList.fromFile(new File(file_name));

            // iterate through the list of merged intervals and add then as GenomeLocs
            ret = new ArrayList<GenomeLoc>();
            for(Interval interval : il.getUniqueIntervals()) {
                ret.add(new GenomeLoc(interval.getSequence(), interval.getStart(), interval.getEnd()));
            }
            return ret;

        } catch (Exception e) {
            try {
                xReadLines reader = new xReadLines(new File(file_name));
                List<String> lines = reader.readLines();
                reader.close();
                String locStr = Utils.join(";", lines);
                logger.debug("locStr: " + locStr);
                ret = parseGenomeLocs(locStr);
                return ret;
            } catch (Exception e2) {
                e2.printStackTrace();
                throw new StingException("Unable to parse out interval file in either format", e);
            }
        }
    }
}

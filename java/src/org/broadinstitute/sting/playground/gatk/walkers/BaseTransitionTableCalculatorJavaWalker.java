package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.LocusWalker;
import org.broadinstitute.sting.gatk.walkers.By;
import org.broadinstitute.sting.gatk.walkers.DataSource;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.GenotypeMetaData;

import java.util.*;
import java.io.PrintStream;
import java.io.PrintWriter;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 12, 2009
 * Time: 2:43:06 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
public class BaseTransitionTableCalculatorJavaWalker extends LocusWalker<ReferenceContextWindow,Integer>{
@Argument(fullName="usePreviousBases", doc="Use previous bases of the reference as part of the calculation, uses the specified number, defaults to 0", required=false)
    int nPreviousBases = 0;
@Argument(fullName="useSecondaryBase",doc="Use the secondary base of a read as part of the calculation", required=false)
    boolean useSecondaryBase = false;
@Argument(fullName="confidentRefThreshold",doc="Set the lod score that defines confidence in ref, defaults to 4", required=false)
    int confidentRefThreshold = 5;
@Argument(fullName="maxNumMismatches",doc="Set the maximum number of mismatches at a locus before choosing not to use it in calculation. Defaults to 1.", required=false)
    int maxNumMismatches = 1;
@Argument(fullName="minMappingQuality", doc ="Set the alignment quality below which to ignore reads; defaults to 30", required = false)
    int minMappingQuality = 30;
@Argument(fullName="minQualityScore", doc = "Set the base quality score below which to ignore bases in the pileup, defaults to 20", required = false)
    int minQualityScore = 20;
@Argument(fullName="usePileupMismatches", doc = "Use the number of mismatches in the pileup as a condition for the table", required=false)
    boolean usePileupMismatches = false;
@Argument(fullName="usePreviousReadBases", doc="Use previous bases of the read as part of the calculation. Will ignore reads if there aren't this many previous bases. Uses the specified number. Defaults to 0", required=false)
    int nPreviousReadBases = 0;

    private UnifiedGenotyper ug;
    private ReferenceContextWindow refWindow;
    private Set<BaseTransitionTable> conditionalTables;

    public void initialize() {
        ug = new UnifiedGenotyper();
        ug.initialize();
        refWindow = new ReferenceContextWindow(nPreviousBases);
        conditionalTables = new HashSet<BaseTransitionTable>();
    }

    public Integer reduceInit() {
        return 0;
    }

    public ReferenceContextWindow map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(),context);
        refWindow.update(ref,pileup,baseIsUsable(tracker,ref,pileup,context));

        return refWindow;
    }

    public Integer reduce ( ReferenceContextWindow map, Integer prevReduce ) {
        if ( map.isValidWindow() ) {
            prevReduce++;
            List<SAMRecord> reads = map.getPileup().getReads();
            List<Integer> offsets = map.getPileup().getOffsets();
            ReferenceContext ref = map.getMiddleReferenceContext();
            // ReadBackedPileup pileup = splitPileupNonref(map.getPileup());
            // List<SAMRecord> reads = pileup.getReads();
            // List<Integer> offsets = pileup.getOffsets();


            // System.out.println("Base and read are usable:");
            // System.out.println("Num Mismatches: "+countMismatches(map.getPileup()));
            // System.out.println("Ref: "+map.getMiddleReferenceContext().getBase()+" Pileup ref: "+map.getPileup().getRef());
            // System.out.println("Pileup: "+map.getPileup().getBases());

            for ( int r = 0; r < reads.size(); r ++ ) {
                if ( Character.toUpperCase(reads.get(r).getReadBases()[offsets.get(r)]) != ref.getBase() ) {
                    // System.out.println("Examining read. Mapping quality is: " + reads.get(r).getMappingQuality());
                    // System.out.println("Base quality is: "+reads.get(r).getBaseQualities()[offsets.get(r)]);
                    // System.out.println("Read base is: "+ (char) reads.get(r).getReadBases()[offsets.get(r)]);
                    // System.out.println("Pileup is: " + map.getPileup().getBases());
                    // if ( ! useRead(reads.get(r),offsets.get(r),ref) ) {
                    //    System.out.println("Read will not be used.");
                    //}
                }
                if ( useRead( reads.get(r), offsets.get(r), ref ) ) {
                    updateTables( reads.get(r), offsets.get(r), map );
                    // prevReduce++;
                }
            }
        }

        return prevReduce;
    }

    public void onTraversalDone( Integer numValidObservedMismatchingReads ) {
        logger.info(numValidObservedMismatchingReads);
        out.print(createHeaderFromConditions());
 	    for ( BaseTransitionTable t : conditionalTables )
		t.print(out);
    }

    public void updateTables ( SAMRecord read, int offset, ReferenceContextWindow map ) {
        List<Comparable> readConditions = buildConditions(read,offset,map);

        boolean createNewTable = true;

        for ( BaseTransitionTable t : conditionalTables ) {
            if ( t.conditionsMatch(readConditions) ) {
                updateTable(t,read,offset,map);
                createNewTable = false;
                break;
            }
        }

        if ( createNewTable ) {
            BaseTransitionTable t = new BaseTransitionTable(readConditions);
            updateTable(t,read,offset,map);
            conditionalTables.add(t);
        }
    }

    public void updateTable(BaseTransitionTable t, SAMRecord r, int o, ReferenceContextWindow map) {
        if ( r.getReadNegativeStrandFlag() ) {
            t.update(BaseUtils.simpleComplement((char) r.getReadBases()[o]),BaseUtils.simpleComplement(map.getMiddleReferenceContext().getBase()));
        } else {
            t.update(r.getReadBases()[o], map.getMiddleReferenceContext().getBase());
        }
    }

    public boolean useRead( SAMRecord read, int offset, ReferenceContext ref ) {

        if ( Character.toUpperCase(read.getReadBases()[offset]) == Character.toUpperCase(ref.getBase()) ) {
            return false;
        } else if ( read.getMappingQuality() <= minMappingQuality ) {
            return false;
        } else if ( ! BaseUtils.isRegularBase( (char) read.getReadBases()[offset]) ) {
            return false;
        } else if ( read.getBaseQualities()[offset] <= minQualityScore ) {
            return false;
        } else {
            return true;
        }
    }

    public List<Comparable> buildConditions( SAMRecord read, int offset, ReferenceContextWindow map ) {
        ArrayList<Comparable> conditions = new ArrayList<Comparable>();

        if ( nPreviousBases > 0 ) {
            if ( ! read.getReadNegativeStrandFlag() )
                conditions.add(map.getForwardRefString());
            else
                conditions.add(map.getReverseRefString());
        }

        if ( useSecondaryBase ) {
            conditions.add(getSecondaryBase(read,offset));
        }

        if ( nPreviousReadBases > 0 ) {
            conditions.add(read.getReadString().substring(offset-nPreviousReadBases,offset));
        }

        if ( usePileupMismatches ) {
            conditions.add(countMismatches(map.getPileup()));
        }

        return conditions;
    }

    public String createHeaderFromConditions() {
        String header = "True_base\tObserved_base";

        if ( nPreviousBases > 0) {
            header = header+"\tPrevious_"+nPreviousBases+"_bases";
        }

        if ( useSecondaryBase ) {
            header = header + "\tSecondary_base";
        }

        if ( nPreviousReadBases > 0 ) {
            header = header + "\tPrevious_"+nPreviousReadBases+"_read_bases";
        }

        if ( usePileupMismatches ) {
            header = header + "\tNumber_of_pileup_mismatches";
        }

        return String.format("%s\t%s%n",header,"Counts");
    }

    public int countMismatches(ReadBackedPileup p) {
        int refM = 0;

        switch (Character.toUpperCase(p.getRef())) {
            case 'A': refM = BasicPileup.countBases(p.getBases()).a;
                break;
            case 'C': refM = BasicPileup.countBases(p.getBases()).c;
                break;
            case 'G': refM = BasicPileup.countBases(p.getBases()).g;
                break;
            case 'T': refM = BasicPileup.countBases(p.getBases()).t;
                break;
            default:
                throw new StingException("Exception in countMismatches - this has been called for a non-reference base. Base character from pileup was " + p.getRef() );
        }

        return p.size()-refM;
    }

    public char getSecondaryBase ( SAMRecord read, int offset ) {
        return BaseUtils.baseIndexToSimpleBase(QualityUtils.compressedQualityToBaseIndex( ( (byte[]) read.getAttribute("SQ") )[offset] ) );
    }

    public boolean baseIsUsable ( RefMetaDataTracker tracker, ReferenceContext ref, ReadBackedPileup pileup, AlignmentContext context ) {
        return  pileupContainsNoNs(pileup) && baseIsConfidentRef(tracker,ref,context) && pileupBelowMismatchThreshold(ref,pileup);
    }

    public boolean pileupBelowMismatchThreshold( ReferenceContext ref, ReadBackedPileup pileup ) {
        return countMismatches(pileup) <= maxNumMismatches;
    }

    public boolean pileupContainsNoNs(ReadBackedPileup pileup) {
	for ( char c : pileup.getBases().toCharArray() ) {
	    if ( c == 'N' ) {
		return false;
	    }
	}
	
	return true;
    }

    public boolean baseIsConfidentRef( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        Pair<List<GenotypeCall>, GenotypeMetaData> calls = ug.map(tracker,ref,context);
        if (calls == null)
            return false;
        return  (! calls.first.get(0).isVariant()) && calls.first.get(0).getNegLog10PError() >= confidentRefThreshold && BaseUtils.isRegularBase(ref.getBase());

    }

}


class BaseTransitionTable {

    private int[][] table;
    private List<Comparable> conditions;

    public BaseTransitionTable(List<Comparable> conditions) {
        table = new int[BaseUtils.BASES.length][BaseUtils.BASES.length];
        for ( int i = 0; i < BaseUtils.BASES.length; i ++ ) {
            for ( int j = 0; j < BaseUtils.BASES.length; j ++ ) {
                table[i][j]=0;
            }
        }

        this.conditions = conditions;
    }

    public boolean conditionsMatch(Object obj) {
        if ( obj == null ) {
            return false;
        } else if ( ! (obj instanceof List) ) {
            return false;
        } else if ( this.numConditions() != ((List)obj).size() ){
            return false;
        } else {
            boolean eq = true;
            ListIterator thisIter = this.getConditionIterator();
            ListIterator thatIter = ((List)obj).listIterator();
            
            while ( thisIter.hasNext() ) {
                eq = eq && thisIter.next().equals(thatIter.next());
            }

            return eq;
        }
    }

    public void print( PrintStream out ) {

        for ( char observedBase : BaseUtils.BASES ) {
            for ( char refBase : BaseUtils.BASES ) {
                String outString = observedBase+"\t"+refBase;
                for ( Comparable c : conditions ) {
                    outString = outString+"\t"+c.toString();
                }
                out.printf("%s\t%d%n",outString,table[BaseUtils.simpleBaseToBaseIndex(observedBase)][BaseUtils.simpleBaseToBaseIndex(refBase)]);
            }
        }
    }

    public void update(char observedBase, char refBase ) {
        //if ( observedBase == refBase ) {
        //    throw new StingException("BaseTransitionTable received equal observed and reference bases, which should not happen.");
        //}
        table[BaseUtils.simpleBaseToBaseIndex(observedBase)][BaseUtils.simpleBaseToBaseIndex(refBase)] ++;
    }

    public void update(byte observed, char ref) {
        update( (char) observed, ref);
    }

    public int numConditions() {
        return conditions.size();
    }

    private Comparable getCondition(int offset) {
        return conditions.get(offset);
    }

    private ListIterator getConditionIterator() {
        return conditions.listIterator();
    }

}



class ReferenceContextWindow {

    protected int windowSize;
    protected int nPrevBases;
    protected LinkedList<ReadBackedPileup> prevAlignments;
    protected LinkedList<ReferenceContext> prevRefs;
    protected LinkedList<Boolean> usePrevious;
    protected boolean initialized;

    public ReferenceContextWindow( int nPrevBases ) {
        windowSize = 2*nPrevBases + 1;
        this.nPrevBases = nPrevBases;
        prevAlignments = new LinkedList<ReadBackedPileup>();
        prevRefs = new LinkedList<ReferenceContext>();
        usePrevious = new LinkedList<Boolean>();
        initialized = false;
    }

    public void update( ReferenceContext ref, ReadBackedPileup pileup, boolean useLocus ) {
        if ( ! initialized ) {
            prevAlignments.add(pileup);
            prevRefs.add(ref);
            usePrevious.add(useLocus);
            if ( prevAlignments.size() == windowSize ) {
                initialized = true;
            }
        } else {
            prevAlignments.removeFirst();
            prevRefs.removeFirst();
            usePrevious.removeFirst();
            prevAlignments.add(pileup);
            prevRefs.add(ref);
            usePrevious.add(useLocus);
        }
    }

    public String getReferenceString() {
        String ref = "";
        for ( ReferenceContext c : prevRefs ) {
            ref = ref + c.getBase();
        }

        return ref;
    }

    public String getForwardRefString() {
        String ref = "";
        for ( ReferenceContext c : prevRefs.subList(0,nPrevBases) ) {
            ref = ref + c.getBase();
        }

        return ref;
    }

    public String getReverseRefString() { // todo -- make sure we want to flip this done (yes we do)
        String ref = "";
        for ( int base = prevRefs.size()-1; base > nPrevBases; base -- ) {
            ref = ref + prevRefs.get(base).getBase();
        }

        return BaseUtils.simpleComplement(ref);
    }

    public ReadBackedPileup getPileup() {
        // because lists are 0-indexed, this returns the alignments
        // to the middle base in the window.
        return prevAlignments.get(nPrevBases);
    }

    public ReferenceContext getMiddleReferenceContext() {
        return prevRefs.get(nPrevBases);
    }

    public boolean isValidWindow() {
        boolean valid;
        if ( ! initialized ) {
            valid = false;
        } else {
            valid = true;
            // check if everything is confident ref
            for ( Boolean b : usePrevious ) {
                if ( !b ) {
                    valid = false;
                    break;
                }
            }
            // if still valid, check distances
            if ( valid ) {
                ListIterator<ReferenceContext> iter = prevRefs.listIterator();
                ReferenceContext prev = iter.next();
                while ( iter.hasNext() ) {
                    ReferenceContext cur = iter.next();

                    if ( cur.getLocus().distance(prev.getLocus()) > 1 ) {
                        valid = false;
                        break;
                    }

                    prev = cur;
                }
            }
        }

        return valid;
    }

    public int getWindowSize() {
        return windowSize;
    }

}

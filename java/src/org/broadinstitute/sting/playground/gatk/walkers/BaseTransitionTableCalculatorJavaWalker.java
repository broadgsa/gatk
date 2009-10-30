package org.broadinstitute.sting.playground.gatk.walkers;

import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.walkers.genotyper.*;
import org.broadinstitute.sting.gatk.refdata.RefMetaDataTracker;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;
import org.broadinstitute.sting.utils.cmdLine.Argument;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.*;

import java.util.*;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.OutputStream;
import java.io.FileNotFoundException;

import net.sf.samtools.SAMRecord;

/**
 * Created by IntelliJ IDEA.
 * User: chartl
 * Date: Oct 12, 2009
 * Time: 2:43:06 PM
 * To change this template use File | Settings | File Templates.
 */
@By(DataSource.REFERENCE)
@Reference(window=@Window(start=-3,stop=3))
public class BaseTransitionTableCalculatorJavaWalker extends LocusWalker<Set<BaseTransitionTable>,Set<BaseTransitionTable>> implements TreeReducible<Set<BaseTransitionTable>> {
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
    @Argument(fullName="useReadGroup", doc="Use the group number of the read as a condition of the table.", required = false)
    boolean useReadGroup = false;
    @Argument(fullName="outputFile", shortName="of", doc="Output to this file rather than standard out. Must be used with -nt.", required = false)
    String outFilePath = null;
    @Argument(fullName="forcePreviousReadBasesToMatchRef", doc="Forces previous read bases to match the reference", required = false)
    boolean readBasesMustMatchRef = false;

    private UnifiedGenotyper ug;
    // private ReferenceContextWindow refWindow;
    // private Set<BaseTransitionTable> conditionalTables;
    private List<Boolean> usePreviousBases;
    private List<GenomeLoc> previousBaseLoci;

    public void initialize() {
        if ( nPreviousBases > 3 || ( nPreviousReadBases > 3 && readBasesMustMatchRef ) ) {
            throw new StingException("You have opted to use a number of previous bases in excess of 3. In order to do this you must change the reference window size in the walker itself.");
        }
        ug = new UnifiedGenotyper();
        UnifiedArgumentCollection uac = new UnifiedArgumentCollection();
        ug.initialize();
        uac.baseModel = BaseMismatchModel.THREE_STATE;
        uac.ALL_BASES = true;
        ug.setUnifiedArgumentCollection(uac);
        // refWindow = new ReferenceContextWindow(nPreviousBases);
        usePreviousBases = new ArrayList<Boolean>();
        previousBaseLoci = new ArrayList<GenomeLoc>();

    }

    public Set<BaseTransitionTable> reduceInit() {
        return new TreeSet<BaseTransitionTable>();
    }

    public Set<BaseTransitionTable> map( RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context ) {
        ReadBackedPileup pileup = new ReadBackedPileup(ref.getBase(),context);
        Set<BaseTransitionTable> newCounts = null;
        //System.out.println(pileup.getBases());
        if ( baseIsUsable(tracker, ref, pileup, context) ) {
            //System.out.println("Pileup will be used");
            if ( previousLociCanBeUsed(usePreviousBases,previousBaseLoci,context.getLocation()) ) {
                for ( int r = 0; r < pileup.getReads().size(); r ++ ) {
                    if ( useRead ( pileup.getReads().get(r), pileup.getOffsets().get(r), ref ) ) {
                        newCounts = updateTables( newCounts, pileup.getReads().get(r), pileup.getOffsets().get(r), ref, pileup );
                    }
                }
            } else {
                updatePreviousBases(usePreviousBases,true,previousBaseLoci,context.getLocation() );
            }
        } else {
            updatePreviousBases( usePreviousBases,false,previousBaseLoci,context.getLocation() );
        }

        return newCounts;
    }

    public Set<BaseTransitionTable> reduce ( Set<BaseTransitionTable> map, Set<BaseTransitionTable> reduce  ) {
        if ( map != null && ! map.isEmpty() ) {
            for ( BaseTransitionTable t : map ) {
                boolean add = true;
                for ( BaseTransitionTable r : reduce ) {
                    if ( r.conditionsMatch(t) ) {
                        r.incorporateTable(t);
                        add = false;
                        break;
                    }
                }
                if ( add ) {
                    reduce.add(t);
                }
            }
        }
        // System.out.println("Reduce: size of TransitionTable set is " + reduce.size() + " -- size of Map: " + (map != null ? map.size() : "null"));
        return reduce;
    }

    public Set<BaseTransitionTable> treeReduce( Set<BaseTransitionTable> reduce1, Set<BaseTransitionTable> reduce2 ) {
        // check to see if this is a truly tree-reducable calculation
        if ( nPreviousBases >= 1 ) {
            String errMsg = "Parallelization cannot be used with UsePreviousBases due to the fact that internal walker data specifies whether a previous reference base is usable or not.";
            String errMsg2 = " This can cause cause concurrency issues and unpredictable behavior when used with parallelization. Either do not specify -nt, or try a the conjunction of ";
            String errMsg3 = "--usePreviousReadBases and --forcePreviousReadBasesToMatchRef.";
            throw new StingException(errMsg+errMsg2+errMsg3);
        }
        return reduce(reduce1,reduce2);
    }

    public void onTraversalDone( Set<BaseTransitionTable> conditionalTables ) {
        PrintStream output;
        if ( outFilePath == null ) {
            output = out;
        } else {
            try {
                output = new PrintStream(outFilePath);
            } catch ( FileNotFoundException e ) {
                throw new StingException("File given as input to -of, "+outFilePath+" could not be opened.",e);
            }
        }
        output.print(createHeaderFromConditions());
        for ( BaseTransitionTable t : conditionalTables )
            t.print(output);
    }

    public void updatePreviousBases(List<Boolean> usage, boolean canUse, List<GenomeLoc> loci, GenomeLoc locus) {
        // early return
        if ( nPreviousBases < 1 ) {
            return;
        }

        if ( usage.size() <= nPreviousBases ) {
            usage.add(canUse);
            loci.add(locus);
        } else {
            usage.remove(0);
            usage.add(canUse);
            loci.remove(0);
            loci.add(locus);
        }
    }

    public boolean previousLociCanBeUsed( List<Boolean> canUse, List<GenomeLoc> loci, GenomeLoc locus ) {
        if ( nPreviousBases < 1 ) {
            return true;
        }
        
        boolean use = true;
        for ( boolean b : canUse ) {
            use = use && b;
        }

        if ( use ) {
            use = use && ( loci.get(0).distance(locus) == 1 ); // truly is PREVIOUS base
        }

        return use;
    }

    public Set<BaseTransitionTable> updateTables ( Set<BaseTransitionTable> tables, SAMRecord read, int offset, ReferenceContext ref, ReadBackedPileup pileup ) {
        List<Comparable> readConditions = buildConditions(read,offset,ref, pileup);
        // System.out.println("Updating table with pileup: "+pileup.getBases()+ ( read.getReadNegativeStrandFlag() ? "-" : "+" ) + "  Quality: "+read.getBaseQualities()[offset] + "  MapQ:  "+read.getMappingQuality());

        if ( tables == null ) {
            tables = new TreeSet<BaseTransitionTable>();
        }

        boolean createNewTable = true;

        for ( BaseTransitionTable t : tables ) {
            if ( t.conditionsMatch(readConditions) ) {
                updateTable(t,read,offset,ref);
                createNewTable = false;
                break;
            }
        }

        if ( createNewTable ) {
            BaseTransitionTable t = new BaseTransitionTable(readConditions);
            updateTable(t,read,offset,ref);
            tables.add(t);
        }

        return tables;
    }

    public void updateTable(BaseTransitionTable t, SAMRecord r, int o, ReferenceContext ref) {
        // System.out.println("Update Table");
        if ( r.getReadNegativeStrandFlag() ) {
            t.update(BaseUtils.simpleComplement((char) r.getReadBases()[o]),BaseUtils.simpleComplement(ref.getBase()));
        } else {
            t.update(r.getReadBases()[o], ref.getBase());
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
        } else if ( useSecondaryBase && read.getAttribute("SQ") == null ) {
            return false;
        } else if ( nPreviousBases >= 1 && previousReadBasesMismatchRef(read, offset, ref) ) {
            return false;
        } else if ( nPreviousReadBases >= 1 && readLacksPreviousBases(read,offset,nPreviousReadBases) ) {
            return false;
        } else if ( nPreviousReadBases >= 1 && readBasesMustMatchRef && previousReadBasesMismatchRef(read, offset, ref) ) {
            return false;
        } else {
            return true;
        }
    }

    public boolean previousReadBasesMismatchRef( SAMRecord read, int offset, ReferenceContext ref ) {
        int c = read.getReadNegativeStrandFlag() ? 1 : -1;
        if ( offset + nPreviousBases*c < 0 ) {
            return true;
        } else if ( offset + nPreviousBases*c > read.getReadLength() ) {
            return true;
        }

        for ( int prevBase = 1; prevBase <= nPreviousBases; prevBase ++ ) {
            if ( Character.toUpperCase(read.getReadBases()[offset + prevBase*c]) != Character.toUpperCase(ref.getBases()[nPreviousBases+1+prevBase*c]) || ! BaseUtils.isRegularBase(ref.getBases()[nPreviousBases+1+prevBase*c])) {
                return true;
            }
        }

        return false;
    }

    public boolean readLacksPreviousBases( SAMRecord read, int offset, int prevBases ) {
        if ( ! read.getReadNegativeStrandFlag() ) {
            return offset - prevBases < 0;
        } else {
            return offset + prevBases + 1 >= read.getReadLength();
        }
    }

    public List<Comparable> buildConditions( SAMRecord read, int offset, ReferenceContext ref, ReadBackedPileup pileup ) {
        ArrayList<Comparable> conditions = new ArrayList<Comparable>();

        if ( nPreviousBases > 0 ) {
            conditions.add(buildRefString(ref,nPreviousBases, ! read.getReadNegativeStrandFlag()));

        }

        if ( useSecondaryBase ) {
            conditions.add(getSecondaryBase(read,offset));
        }

        if ( nPreviousReadBases > 0 ) {
            conditions.add(buildReadString(read, offset, nPreviousReadBases));
        }

        if ( usePileupMismatches ) {
            conditions.add(countMismatches(pileup));
        }

        if ( useReadGroup ) {
            conditions.add(read.getReadGroup().getReadGroupId());
        }

        return conditions;
    }

    public String buildRefString(ReferenceContext ref, int bases, boolean forwardRead) {
        if ( forwardRead ) {
            return ( new String(ref.getBases()) ).substring(0,nPreviousBases-1);
        } else {
            return BaseUtils.simpleReverseComplement( ( new String(ref.getBases()) ).substring(nPreviousBases+1) );
        }
    }

    public String buildReadString( SAMRecord read, int offset, int nPreviousReadBases ) {
        if ( ! read.getReadNegativeStrandFlag() ) {
            return read.getReadString().substring(offset-nPreviousReadBases,offset);
        } else {
            return BaseUtils.simpleReverseComplement( read.getReadString().substring(offset+1,offset+nPreviousReadBases+1) );
        }
    }

    public String createHeaderFromConditions() {
        String header = "Observed_base\tTrue_base";

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

        if ( useReadGroup ) {
            header = header + "\tRead_group";
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
        Pair<List<Genotype>, GenotypeMetaData> calls = ug.map(tracker,ref,context);
        if (calls == null || calls.first == null)
            return false;
        return  (! calls.first.get(0).isVariant(ref.getBase())) && calls.first.get(0).getNegLog10PError() > confidentRefThreshold && BaseUtils.isRegularBase(ref.getBase());

    }

}


class BaseTransitionTable implements Comparable {

    /*
     * no direct manipulation of these objects ever
     */
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
        } else if ( obj instanceof BaseTransitionTable ) {
            return ((BaseTransitionTable) obj).conditionsMatch(conditions);
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


    public int compareTo(Object obj) {
        if ( ! ( obj instanceof BaseTransitionTable ) ) {
            return -1;
        } else {
            BaseTransitionTable t = (BaseTransitionTable) obj;
            if ( this.conditionsMatch(t.conditions) ) {
                return 0;
            } else {
                if ( this.numConditions() == t.numConditions() ) {
                    ListIterator<Comparable> thisIter = this.conditions.listIterator();
                    ListIterator<Comparable> thatIter = t.conditions.listIterator();
                    int g = 0;
                    do {
                        g = thisIter.next().compareTo(thatIter.next());
                    } while ( g == 0 );

                    return g;
                    
                } else {
                    return (this.numConditions() > t.numConditions() ) ? 1 : -1;
                }
            }
        }

    }

    public void print( PrintStream out ) {
        StringBuilder s = new StringBuilder();
        for ( char observedBase : BaseUtils.BASES ) {
            for ( char refBase : BaseUtils.BASES ) {
                s.append(String.format("%s\t%s",observedBase,refBase));
                for ( Comparable c : conditions ) {
                    s.append(String.format("\t%s",c.toString()));
                }
                s.append(String.format("\t%d%n", table[BaseUtils.simpleBaseToBaseIndex(observedBase)][BaseUtils.simpleBaseToBaseIndex(refBase)]));
            }
        }

        out.print(s.toString());
    }

    public void update(char observedBase, char refBase ) {
        //if ( observedBase == refBase ) {
        //    throw new StingException("BaseTransitionTable received equal observed and reference bases, which should not happen.");
        //}
        // System.out.println("Table updating: Observed Base: "+observedBase+" Ref base: "+refBase);
        table[BaseUtils.simpleBaseToBaseIndex(observedBase)][BaseUtils.simpleBaseToBaseIndex(refBase)]++;
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

    public void incorporateTable(BaseTransitionTable t) {
        for ( int i = 0; i < BaseUtils.BASES.length; i ++ ) {
            for ( int j = 0; j < BaseUtils.BASES.length; j ++ ) {
                table[i][j] += t.observationsOf(i,j);
            }
        }
    }

    public int observationsOf( int observedBaseIndex, int referenceBaseIndex ) {
        return table[observedBaseIndex][referenceBaseIndex];
    }

}
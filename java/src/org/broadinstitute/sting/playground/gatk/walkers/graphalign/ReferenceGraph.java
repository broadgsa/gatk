package org.broadinstitute.sting.playground.gatk.walkers.graphalign;

import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleDirectedGraph;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.utils.genotype.Variation;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.IntervalTree;
import net.sf.samtools.util.StringUtil;

import java.util.*;
import java.io.Serializable;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.ObjectInputStream;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Oct 22, 2009
 * Time: 2:15:54 PM
 * To change this template use File | Settings | File Templates.
 */
class ReferenceGraph extends SimpleDirectedGraph<Fragment, DefaultEdge> implements Serializable {
    final private static boolean USE_IT = true;
    final private static boolean THROW_ERRORS_ON_BAD_INPUTS = false;

    private boolean DEBUG = false;
    int nSkippedIndels = 0;
    int nSkippedBecauseOfContinguousVariation = 0;
    int nBadPolymorphisms = 0;
    int nMultiStateAlleles = 0;

    GenomeLoc initialLoc = null;

    private transient IntervalTree<Fragment> loc2Fragment = new IntervalTree<Fragment>();

    public ReferenceGraph(boolean printDebuggingInfo) {
        super(DefaultEdge.class);
        DEBUG = printDebuggingInfo;
    }

    public ReferenceGraph() {
        this(false);
    }

    public void setDebugPrinting(boolean enable) {
        this.DEBUG = enable;
    }

    public void bindRefenceSequence(ReferenceSequence seq) {
        GenomeLoc refSeqLoc = GenomeLocParser.createGenomeLoc(seq.getContigIndex(), 1, seq.length());
        String refString = StringUtil.bytesToString(seq.getBases()).toUpperCase(); 
        Fragment frag = new Fragment(refSeqLoc, 1, StringUtil.stringToBytes(refString));
        addFragment(frag);
        initialLoc = refSeqLoc;
    }

    private void addFragmentToIntervalTree(Fragment frag) {
        loc2Fragment.put((int)frag.getLocation().getStart(), (int)frag.getLocation().getStop(), frag);
    }

    private void addFragment(Fragment frag) {
        addFragmentToIntervalTree(frag);
        addVertex(frag);
    }

    private void removeFragment(Fragment frag) {
        loc2Fragment.remove((int)frag.getLocation().getStart(), (int)frag.getLocation().getStop());
        removeVertex(frag);
    }

    public void validateGraph() {
        for ( Fragment v : this.vertexSet() ) {
            if ( this.inDegreeOf(v) == 0 && v.getLocation().getStart() != initialLoc.getStart() ) {
                throw new StingException(String.format("Fragment %s has no incoming edges but isn't at the start of the contig %s", v, initialLoc));
            }
            if ( this.outDegreeOf(v) == 0 && v.getLocation().getStop() != initialLoc.getStop() ) {
                throw new StingException(String.format("Fragment %s has no outgoing edges but isn't at the end of the contig %s", v, initialLoc));
            }
        }
        
        //System.out.printf("Passed validation: %s%n", this.toBriefString());
    }

    private void rebuildIntervalTree() {
        if ( DEBUG ) System.out.printf("rebuilding IntervalTree()%n");
        for ( Fragment v : this.vertexSet() ) {
            if ( DEBUG ) System.out.printf("  adding interval tree: %s%n", v);
            addFragmentToIntervalTree(v);
        }
    }

    private boolean allelesAreInExcisedFragment(Fragment cut, List<String> alleles) {
        boolean foundRef = false;
        for ( String allele : alleles ) {
            if ( allele.equals(cut.getBases()) ) 
                foundRef = true;
        }

        if ( ! foundRef && THROW_ERRORS_ON_BAD_INPUTS )
            throw new StingException(String.format("Polymorphic alleles %s do not contain the reference sequence %s", alleles, cut.getBases()));
        
        return foundRef;
    }

    public void addVariation(Variation variant, GenomeLoc loc, List<String> alleles) {
        if ( DEBUG ) System.out.printf("addVariation(%s, %s)%n", loc, alleles);
        //validateGraph();

        if ( variant.isSNP() ) {
            Fragment frag = getContainingFragment(loc);

            if ( frag == null ) {
                nMultiStateAlleles++;
                return;
            }

            if ( ! allelesAreInExcisedFragment(subFragment(frag, loc, 1), alleles)) {
                nBadPolymorphisms++;
                return;
            }

            List<Fragment> split = exciseFragment(frag, loc);
            if ( split != null ) {
                Fragment left = split.get(0);
                Fragment cut = split.get(1);
                Fragment right = split.get(2);

                if ( DEBUG ) System.out.printf("  cutFrag(%s, %s)%n", loc, cut);

                for ( String allele : alleles ) {
                    byte[] bases = StringUtil.stringToBytes(allele);
                    double freq = 1.0 / alleles.size();
                    Fragment alleleFrag = new Fragment(loc, freq, 0, bases.length, bases);
                    if ( DEBUG ) System.out.printf("  Creating allele fragment %s%n", alleleFrag);
                    addFragment(alleleFrag);
                    if ( left != null ) addEdge(left, alleleFrag);
                    if ( right != null ) addEdge(alleleFrag, right);
                }
            } else {
                nSkippedBecauseOfContinguousVariation++;
            }
        } else {
            nSkippedIndels++;
        }
    }


    private Fragment subFragment(Fragment frag, GenomeLoc loc, double freq ) {
        return new Fragment(loc, 1, frag.getUnderlyingBases());
    }

    public List<Fragment> exciseFragment(Fragment frag, GenomeLoc loc) {
        if ( DEBUG ) System.out.printf("  exciseFragment(%s, %s)%n", frag, loc);
        GenomeLoc fragLoc = frag.getLocation();

        Fragment cut = subFragment(frag, loc, 1);

        Set<DefaultEdge> inToFrag = incomingEdgesOf(frag);
        Set<DefaultEdge> outOfFrag = outgoingEdgesOf(frag);

        Fragment left = null;
        if ( fragLoc.getStart() == loc.getStart() ) {
            if ( ! inToFrag.isEmpty() ) {
                if ( THROW_ERRORS_ON_BAD_INPUTS )
                    throw new StingException(String.format("Attempting to create a variation at the start of a fragment %s %s", frag, loc));
                return null;
            }
        } else {
            GenomeLoc leftLoc = GenomeLocParser.createGenomeLoc(fragLoc.getContigIndex(), fragLoc.getStart(), loc.getStart()-1);
            left = new Fragment(leftLoc, 1, frag.getUnderlyingBases());
            addFragment(left);

            for ( DefaultEdge e : inToFrag ) {
                addEdge(getEdgeSource(e), left);
            }

            removeAllEdges(inToFrag);
        }

        Fragment right = null;
        if ( fragLoc.getStop() == loc.getStop() ) {
            if ( ! outOfFrag.isEmpty() ) {
                throw new StingException(String.format("Attempting to create a variation at the end of a fragment %s %s", frag, loc));
            }
        } else {
            GenomeLoc rightLoc = GenomeLocParser.createGenomeLoc(fragLoc.getContigIndex(), loc.getStop()+1, fragLoc.getStop());
            right = new Fragment(rightLoc, 1, frag.getUnderlyingBases());
            addFragment(right);

            for ( DefaultEdge e : outOfFrag ) {
                addEdge(right, getEdgeTarget(e));
            }

            removeAllEdges(outOfFrag);
        }

        if ( DEBUG ) System.out.printf("    removing %s%n", frag);
        removeFragment(frag);
        if ( DEBUG ) System.out.printf("    returning left=%s right=%s%n", left, right);
        return Arrays.asList(left, cut, right);
    }

    public Fragment getContainingFragment(GenomeLoc loc) {
        Fragment frag = USE_IT ? getContainingFragmentIT(loc) : getContainingFragmentG(loc);


        if ( frag == null ) {
            if ( THROW_ERRORS_ON_BAD_INPUTS )
                throw new StingException("No spanning fragment was found for " + loc);
            else
                return null;
        }
        else if ( frag.getLocation().getStart() > loc.getStart() || frag.getLocation().getStop() < loc.getStop() )
            throw new StingException("BUG: bad spanning fragment found for " + loc + " was " + frag.getLocation() );
        else
            return frag;
    }

    public Fragment getContainingFragmentG(GenomeLoc loc) {
        for ( Fragment v : this.vertexSet() ) {
            if ( v.getLocation().containsP(loc) ) {
                return v;
            }
        }

        return null;
    }

    public Fragment getContainingFragmentIT(GenomeLoc loc) {
        IntervalTree.Node<Fragment> node = loc2Fragment.minOverlapper((int)loc.getStart(), (int)loc.getStop());
        if ( node == null )
            return null;
        else
            return node.getValue();
    }

    public Collection<Fragment> getStartingFragment(GenomeLoc loc) {
        Collection<Fragment> frags = USE_IT ? getStartingFragmentIT(loc) : getStartingFragmentG(loc);
        //Collection<Fragment> frags = getStartingFragmentTest(loc);

        if ( frags == null || frags.size() == 0 )
            throw new StingException("No fragment contains location start of " + loc);
        if ( frags.size() == 1 && MathUtils.compareDoubles(frags.iterator().next().getFrequency(), 1.0) != 0 ) {
            Fragment bad = frags.iterator().next();
            throw new StingException(String.format("Only one fragment was found but it's frequency < 1 %s with %e", bad, bad.getFrequency()));
        }
        else
            return frags;
    }

    public Collection<Fragment> getStartingFragmentTest(GenomeLoc loc) {
        Collection<Fragment> fragsFromIT = getStartingFragmentIT(loc);
        Collection<Fragment> fragsFromG = getStartingFragmentG(loc);

        if ( fragsFromIT.size() != fragsFromG.size() ) {
            throw new StingException(String.format("Fragment sizes differ %d from IntervalTree, %d from graph", fragsFromIT.size(), fragsFromG.size()));
        }

        return USE_IT && false ? fragsFromIT : fragsFromG;
    }

    public Collection<Fragment> getStartingFragmentIT(GenomeLoc loc) {
        Collection<Fragment> frags = new HashSet<Fragment>();

        Iterator<IntervalTree.Node<Fragment>> it = loc2Fragment.overlappers((int)loc.getStart(), (int)loc.getStart());
        IntervalTree.Node<Fragment> node = null;
        while ( it.hasNext() ) {
            node = it.next();
            frags.add(node.getValue());
        }

        // todo -- painful bug work around -- should be removed
        if ( frags.size() == 1 && MathUtils.compareDoubles(node.getValue().getFrequency(), 1.0) != 0 ) {
            System.out.printf(">>> Using IT workaround at %s <<<%n", loc);
            return getStartingFragmentG(loc);
        }
        
        return frags;
//        IntervalTree.Node<Fragment> node = loc2Fragment.minOverlapper((int)loc.getStart(), (int)loc.getStart());
//        if ( node == null )
//            return null;
//        else
//            return node.getValue();
    }

    public Collection<Fragment> getStartingFragmentG(GenomeLoc loc) {
        Collection<Fragment> frags = new HashSet<Fragment>();
        for ( Fragment v : this.vertexSet() ) {
            //if ( v.getLocation().getStart() < loc.getStop() )
            //    System.out.printf("Checking %s vs. %s%n", loc, v.getLocation());
            if ( v.getLocation().containsStartPosition(loc.getStart()) ) {
            //    System.out.printf(" Adding %s%n", v.getLocation());
                frags.add(v);
            }
        }

        return frags;
    }

    public Set<Fragment> outgoingFragments( Fragment frag ) {
        Set<Fragment> outgoingFrags = new HashSet<Fragment>();
        
        for ( DefaultEdge e : outgoingEdgesOf(frag) ) {
            outgoingFrags.add(getEdgeTarget(e));
        }

        if ( outgoingFrags.size() == 0 && frag.getLocation().getStop() != initialLoc.getStop() ) {

        }

        return outgoingFrags;
    }

    public String toString() {
        StringBuilder s = new StringBuilder();

        for ( Fragment v : this.vertexSet() ) {
            s.append(String.format("Fragment: %s%n", v.toString()));
            for ( DefaultEdge e : this.incomingEdgesOf(v) ) {
                s.append(String.format("  [IN FROM] %s%n", this.getEdgeSource(e)));
            }
            for ( DefaultEdge e : this.outgoingEdgesOf(v) ) {
                s.append(String.format("  [OUT TO ] %s%n", this.getEdgeTarget(e)));
            }
        }

        return s.toString();
    }

    public String toBriefString() {
        return String.format("GraphRef: %d fragments, %d edges, skipped %d contingous variants, %d indels, %d polymorphisms w/o ref allele, %d multi-state",
                this.vertexSet().size(), this.edgeSet().size(), nSkippedBecauseOfContinguousVariation, nSkippedIndels, nBadPolymorphisms, nMultiStateAlleles);
    }

    // --------------------------------------------------------------------------------------------------------------
    //
    // serialization
    //
    // --------------------------------------------------------------------------------------------------------------
    private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
        //always perform the default de-serialization first
        stream.defaultReadObject();
        loc2Fragment = new IntervalTree<Fragment>();
        rebuildIntervalTree();
    }
}

class Fragment implements Serializable {
    GenomeLoc loc = null;   // Index position of this fragment into the reference
    int offset, stop;
    double freq = -1;
    byte[] bases = null;

    public Fragment( GenomeLoc loc, double freq, int offset, int stop, byte[] bases ) {
        this.loc = loc;
        this.freq = freq;
        this.bases = bases;
        this.offset = offset;
        this.stop = stop;
    }

    public Fragment( GenomeLoc loc, double freq, byte[] bases ) {
        this(loc, freq, (int)loc.getStart()-1, (int)loc.getStop(), bases);
    }

    public String toString() {
        //return String.format("%s:%.2f:%s", loc.toString(), getFrequency(), getBases());
        return String.format("%s:%.2f", loc.toString(), getFrequency());
    }

    public GenomeLoc getLocation() {
        return loc;
    }

    public double getFrequency() {
        return freq;
    }

    public int getUnderlyingOffset() { return offset; }
    public int getStop() { return stop; }
    public int getLength() { return getStop() - getUnderlyingOffset(); }

    public byte[] getUnderlyingBases() {
        return bases;
    }

    /**
     * how many bases over in the fragment are we over in this fragment?
     *
     * @param loc
     * @return
     */
    public int getFragOffsetFrom(GenomeLoc loc) {
        // todo -- ignores contigs -- can we fix this?
        if ( getLocation().getStart() > loc.getStart() )
            throw new StingException("BUG: Request for offset from " + loc + " in frag at " + getLocation() + " but this is beyond the location of the fragment");
        return (int)(loc.getStart() - getLocation().getStart());
    }

    public int getBaseLengthFrom( int fragOffset, int maxLength ) {
        int fragRemaining = getLength() - fragOffset;

        if ( fragRemaining < 0 )
            throw new StingException("BUG: Request for length from offset " + fragOffset + " but this is longer than the fragment itself");
        
        return Math.min(fragRemaining, maxLength);
    }

    public byte getBase(int fragOffset) {
        return bases[getUnderlyingOffset() + fragOffset];
    }

    public String getBases() {
        return StringUtil.bytesToString(getUnderlyingBases(), getUnderlyingOffset(), getLength());
    }
}

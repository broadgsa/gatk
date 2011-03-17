package org.broadinstitute.sting.playground.gatk.walkers.assembly;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import org.jgrapht.graph.*;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 14, 2011
 */
public class SimpleDeBruijnAssembler extends LocalAssemblyEngine {

    private static final boolean DEBUG = true;

    // k-mer length
    private static final int KMER_LENGTH = 19;

    // minimum base quality required in a contiguous stretch of a given read to be used in the assembly
    private static final int MIN_BASE_QUAL_TO_USE = 20;

    // minimum clipped sequence length to consider using
    private static final int MIN_SEQUENCE_LENGTH = 30;

    // the deBruijn graph object
    private DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = null;

    // simple node class for storing kmer sequences
    protected class DeBruijnVertex {
        protected byte[] sequence;

        public DeBruijnVertex(byte[] sequence) {
            this.sequence = sequence;
        }

        public boolean equals(DeBruijnVertex v) {
            return Arrays.equals(sequence, v.sequence);
        }
    }

    // simple edge class for connecting nodes in the graph
    protected class DeBruijnEdge {
        private int multiplicity;

        public DeBruijnEdge() {
            multiplicity = 1;
        }

        public int getMultiplicity() {
            return multiplicity;
        }

        public void setMultiplicity(int value) {
            multiplicity = value;
        }
    }


    public SimpleDeBruijnAssembler(PrintStream out, IndexedFastaSequenceFile referenceReader) {
        super(out, referenceReader);
    }

    public void runLocalAssembly(List<SAMRecord> reads) {

        // reset the graph
        graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        // clip the reads to get just the base sequences we want
        List<byte[]> sequences = clipReads(reads);

        // create the graph
        createDeBruijnGraph(sequences);

        // assign reads to the graph
        assignReadsToGraph(sequences);
    }

    // This method takes the base sequences from the SAM records and clips both ends so that
    // soft-clipped bases - plus all bases until the first Q20 (on either end) - are removed.
    // Clipped sequences that are overly clipped are not used.
    private List<byte[]> clipReads(List<SAMRecord> reads) {
        List<byte[]> sequences = new ArrayList<byte[]>(reads.size());

        for ( SAMRecord read : reads ) {
            byte[] sequencedReadBases = read.getReadBases();
            byte[] sequencedBaseQuals = read.getBaseQualities();

            int fromIndex = 0;
            boolean sawQ20 = false;
            for ( CigarElement ce : read.getCigar().getCigarElements() ) {
                if ( sawQ20 )
                    break;

                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        // skip soft-clipped bases
                        fromIndex += elementLength;
                        break;
                    case M:
                    case I:
                        for (int i = 0; i < elementLength; i++) {
                            if ( sequencedBaseQuals[fromIndex] >= MIN_BASE_QUAL_TO_USE ) {
                                sawQ20 = true;
                                break;
                            }
                            fromIndex++;
                        }
                    default:
                        break;
                }
            }

            int toIndex = sequencedReadBases.length - 1;
            sawQ20 = false;
            for (int ceIdx = read.getCigar().numCigarElements() - 1; ceIdx >= 0; ceIdx--) {
                if ( sawQ20 )
                    break;

                CigarElement ce = read.getCigar().getCigarElement(ceIdx);
                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        // skip soft-clipped bases
                        toIndex -= elementLength;
                        break;
                    case M:
                    case I:
                        for (int i = 0; i < elementLength; i++) {
                            if ( sequencedBaseQuals[toIndex] >= MIN_BASE_QUAL_TO_USE ) {
                                sawQ20 = true;
                                break;
                            }
                            toIndex--;
                        }
                    default:
                        break;
                }
            }

            int sequenceLength = toIndex - fromIndex + 1;
            if ( sequenceLength < MIN_SEQUENCE_LENGTH ) {
                if ( DEBUG ) System.out.println("Read " + read.getReadName() + " is overly clipped");
                continue;
            }

            //if ( DEBUG ) {
            //    System.out.println("Read " + read.getReadName() + ": " + read.getCigar());
            //    System.out.println("Clipping from " + fromIndex + " to " + toIndex);
            //}

            byte[] sequence = new byte[sequenceLength];
            System.arraycopy(sequencedReadBases, fromIndex, sequence, 0, sequenceLength);
            sequences.add(sequence);
        }

        return sequences;
    }

    private void createDeBruijnGraph(List<byte[]> reads) {

        // create the graph
        createGraphFromSequences(reads);

        // cleanup graph by merging nodes
        concatenateNodes();

        // cleanup the node sequences so that they print well
        cleanupNodeSequences();

        if ( DEBUG )
            printGraph();
    }

    private void createGraphFromSequences(List<byte[]> reads) {
        for ( byte[] sequence : reads ) {

            final int kmersInSequence = sequence.length - KMER_LENGTH + 1;
            for (int i = 0; i < kmersInSequence - 1; i++) {
                // get the kmers
                byte[] kmer1 = new byte[KMER_LENGTH];
                System.arraycopy(sequence, i, kmer1, 0, KMER_LENGTH);
                byte[] kmer2 = new byte[KMER_LENGTH];
                System.arraycopy(sequence, i+1, kmer2, 0, KMER_LENGTH);

                addEdgeToGraph(kmer1, kmer2);

                // TODO -- eventually, we'll need to deal with reverse complementing the sequences
            }
        }
    }

    private void addEdgeToGraph(byte[] kmer1, byte[] kmer2) {

        DeBruijnVertex v1 = addToGraphIfNew(kmer1);
        DeBruijnVertex v2 = addToGraphIfNew(kmer2);

        Set<DeBruijnEdge> edges = graph.outgoingEdgesOf(v1);
        DeBruijnEdge targetEdge = null;
        for ( DeBruijnEdge edge : edges ) {
            if ( graph.getEdgeTarget(edge).equals(v2) ) {
                targetEdge = edge;
                break;
            }
        }

        if ( targetEdge == null )
            graph.addEdge(v1, v2, new DeBruijnEdge());
        else
            targetEdge.setMultiplicity(targetEdge.getMultiplicity() + 1);
    }

    private DeBruijnVertex addToGraphIfNew(byte[] kmer) {
        // the graph.containsVertex() method is busted, so here's a hack around it

        DeBruijnVertex newV = new DeBruijnVertex(kmer);
        for ( DeBruijnVertex v : graph.vertexSet() ) {
            if ( v.equals(newV) )
                return v;
        }

        graph.addVertex(newV);
        return newV;
    }

    private void concatenateNodes() {

        while ( true ) {
            boolean graphWasModified = false;

            Set<DeBruijnVertex> vertexSet = graph.vertexSet();
            // convert to array because results of the iteration on a set are undefined when the graph is modified
            ArrayList<DeBruijnVertex> vertices = new ArrayList<DeBruijnVertex>(vertexSet);

            for (int i = 0; i < vertices.size(); i++) {

                DeBruijnVertex v1 = vertices.get(i);

                // try to merge v1 -> v2
                if ( graph.outDegreeOf(v1) == 1 ) {
                    DeBruijnEdge edge = graph.outgoingEdgesOf(v1).iterator().next();
                    DeBruijnVertex v2 = graph.getEdgeTarget(edge);

                    if ( graph.inDegreeOf(v2) == 1 ) {
                        mergeVertices(v1, v2);
                        graphWasModified = true;
                        break;
                    }
                }

                // try to merge v2 -> v1
                if ( graph.inDegreeOf(v1) == 1 ) {
                    DeBruijnEdge edge = graph.incomingEdgesOf(v1).iterator().next();
                    DeBruijnVertex v2 = graph.getEdgeSource(edge);

                    if ( graph.outDegreeOf(v2) == 1 ) {
                        mergeVertices(v2, v1);
                        graphWasModified = true;
                        break;
                    }
                }
            }

            if ( !graphWasModified )
                break;
        }
    }

    private void mergeVertices(DeBruijnVertex V1, DeBruijnVertex V2) {
        // (Vx -> V1 -> V2 -> Vy)
        //     should now be
        //   (Vx -> V12 -> Vy)

        // create V12
        int additionalSequenceFromV2 = V2.sequence.length - KMER_LENGTH + 1;
        byte[] newKmer = new byte[V1.sequence.length + additionalSequenceFromV2];
        System.arraycopy(V1.sequence, 0, newKmer, 0, V1.sequence.length);
        System.arraycopy(V2.sequence, KMER_LENGTH - 1, newKmer, V1.sequence.length, additionalSequenceFromV2);
        DeBruijnVertex V12 = new DeBruijnVertex(newKmer);
        graph.addVertex(V12);

        // copy edges coming from Vx to V12
        Set<DeBruijnEdge> Ex = graph.incomingEdgesOf(V1);
        for ( DeBruijnEdge edge : Ex ) {
            DeBruijnVertex Vx = graph.getEdgeSource(edge);
            DeBruijnEdge newEdge = new DeBruijnEdge();
            newEdge.setMultiplicity(edge.getMultiplicity());
            graph.addEdge(Vx, V12, newEdge);
        }

        // copy edges going to Vy from V12
        Set<DeBruijnEdge> Ey = graph.outgoingEdgesOf(V2);
        for ( DeBruijnEdge edge : Ey ) {
            DeBruijnVertex Vy = graph.getEdgeTarget(edge);
            DeBruijnEdge newEdge = new DeBruijnEdge();
            newEdge.setMultiplicity(edge.getMultiplicity());
            graph.addEdge(V12, Vy, newEdge);
        }

        // remove V1 and V2 and their associated edges
        graph.removeVertex(V1);
        graph.removeVertex(V2);
    }

    private void cleanupNodeSequences() {

        for ( DeBruijnVertex v :  graph.vertexSet() ) {

            // remove the first k-1 bases of the kmers
            if ( graph.inDegreeOf(v) > 0 )
                removeKmerPrefix(v);

            // move common suffixes from incoming nodes to this one
            if ( graph.inDegreeOf(v) > 1 )  {
                Set<DeBruijnVertex> connectedVs = new HashSet<DeBruijnVertex>();
                for ( DeBruijnEdge edge : graph.incomingEdgesOf(v) )
                    connectedVs.add(graph.getEdgeSource(edge));
                propagateCommonSuffix(v, connectedVs);
            }
        }

        removeEmptyNodes();
    }

    private void removeEmptyNodes() {

        // remember that results of an iteration on a set are undefined when the graph is modified
        while ( true ) {

            boolean graphWasModified = false;
            for ( DeBruijnVertex v :  graph.vertexSet() ) {
                if ( v.sequence.length == 0 ) {

                    Set<DeBruijnEdge> incoming = graph.incomingEdgesOf(v);
                    Set<DeBruijnEdge> outgoing = graph.outgoingEdgesOf(v);

                    // make edges from all incoming nodes to all outgoing nodes
                    for ( DeBruijnEdge Ex : incoming ) {
                        DeBruijnVertex Vx = graph.getEdgeSource(Ex);
                        for ( DeBruijnEdge Ey : outgoing ) {
                            DeBruijnVertex Vy = graph.getEdgeTarget(Ey);

                            DeBruijnEdge newEdge = new DeBruijnEdge();
                            newEdge.setMultiplicity(Ex.getMultiplicity());
                            graph.addEdge(Vx, Vy, newEdge);
                        }
                    }

                    // remove v and its associated edges
                    graph.removeVertex(v);

                    graphWasModified = true;
                    break;
                }
            }

            if ( !graphWasModified )
                break;
        }
    }

    private void removeKmerPrefix(DeBruijnVertex v) {
        int newLength = v.sequence.length - KMER_LENGTH + 1;
        byte[] newSequence = new byte[newLength];
        System.arraycopy(v.sequence, KMER_LENGTH - 1, newSequence, 0, newLength);
        v.sequence = newSequence;
    }

    private void propagateCommonSuffix(DeBruijnVertex Vx, Set<DeBruijnVertex> incoming) {

        // find the common matching suffix
        byte[] match = null;
        for ( DeBruijnVertex v : incoming ) {
            if ( match == null ) {
                match = v.sequence;
            } else {
                int idx = 0;
                while ( idx < match.length && idx < v.sequence.length && match[match.length - idx - 1] == v.sequence[v.sequence.length - idx - 1] )
                    idx++;

                if ( idx < match.length ) {
                    match = new byte[idx];
                    System.arraycopy(v.sequence, v.sequence.length - idx, match, 0, idx);
                }
            }
        }

        // if there is a common suffix...
        if ( match != null && match.length > 0 ) {

            // remove it from the end of the incoming nodes
            for ( DeBruijnVertex v : incoming ) {
                int newLength = v.sequence.length - match.length;
                byte[] newSequence = new byte[newLength];
                System.arraycopy(v.sequence, 0, newSequence, 0, newLength);
                v.sequence = newSequence;
            }

            // and put it at the front of this node
            byte[] newSequence = new byte[Vx.sequence.length + match.length];
            System.arraycopy(match, 0, newSequence, 0, match.length);
            System.arraycopy(Vx.sequence, 0, newSequence, match.length, Vx.sequence.length);
            Vx.sequence = newSequence;
        }
    }

    private void printGraph() {

        for ( DeBruijnVertex source : graph.vertexSet() ) {
            getOutputStream().print(new String(source.sequence) + " -> ");
            for ( DeBruijnEdge edge : graph.outgoingEdgesOf(source) ) {
                getOutputStream().print(new String(graph.getEdgeTarget(edge).sequence) + " (" + edge.getMultiplicity() + "), ");
            }
            getOutputStream().println();
        }
        getOutputStream().println("------------\n");
    }

    private void assignReadsToGraph(List<byte[]> reads) {

        // TODO -- implement me

    }
}

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
    private static final int KMER_LENGTH = 9;

    // minimum base quality required in a contiguous stretch of a given read to be used in the assembly
    private static final int MIN_BASE_QUAL_TO_USE = 20;

    // minimum clipped sequence length to consider using
    private static final int MIN_SEQUENCE_LENGTH = 25;

    // minimum multiplicity to consider using
    // private static final int MIN_MULTIPLICITY_TO_USE = 2;

    // FOR DEBUGGING
    private int numReadsToUse = -1;

    // the deBruijn graph object
    private DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph = null;


    public SimpleDeBruijnAssembler(PrintStream out, IndexedFastaSequenceFile referenceReader, int numReadsToUse) {
        super(out, referenceReader);
        this.numReadsToUse = numReadsToUse;
    }

    public void runLocalAssembly(List<SAMRecord> reads) {

        // reset the graph
        graph = new DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge>(DeBruijnEdge.class);

        // clip the reads to get just the base sequences we want
        List<byte[]> sequences = clipReads(reads);

        // create the graph
        createDeBruijnGraph(sequences);

        // find the 2 best paths in the graph
        findBestPaths();

        // assign reads to the graph
        assignReadsToGraph(sequences);
    }

    // This method takes the base sequences from the SAM records and pulls
    // out runs of bases that are not soft-clipped and are all at least Q20s.
    // Clipped sequences that are overly clipped are not used.
    private List<byte[]> clipReads(List<SAMRecord> reads) {
        List<byte[]> sequences = new ArrayList<byte[]>(reads.size());

        int counter = 0;

        for ( SAMRecord read : reads ) {

            // for debugging
            if ( numReadsToUse >= 0 && ++counter > numReadsToUse ) {
                System.out.println("Stopping before read: " + read.getReadName() + " at " + read.getAlignmentStart());
                break;
            }

            byte[] sequencedReadBases = read.getReadBases();
            byte[] sequencedBaseQuals = read.getBaseQualities();

            int curIndex = 0, firstQ20Index = -1;

            for ( CigarElement ce : read.getCigar().getCigarElements() ) {

                int elementLength = ce.getLength();
                switch ( ce.getOperator() ) {
                    case S:
                        // skip soft-clipped bases
                        curIndex += elementLength;
                        break;
                    case M:
                    case I:
                        for (int i = 0; i < elementLength; i++) {
                            if ( sequencedBaseQuals[curIndex] >= MIN_BASE_QUAL_TO_USE ) {
                                if ( firstQ20Index == -1 )
                                    firstQ20Index = curIndex;
                            } else if ( firstQ20Index != -1 ) {
                                int sequenceLength = curIndex - firstQ20Index;
                                if ( sequenceLength > MIN_SEQUENCE_LENGTH ) {
                                    byte[] sequence = new byte[sequenceLength];
                                    System.arraycopy(sequencedReadBases, firstQ20Index, sequence, 0, sequenceLength);
                                    sequences.add(sequence);
                                }
                                firstQ20Index = -1;
                            }
                            curIndex++;
                        }
                    case N:
                        // TODO -- implement me (cut the sequence)
                    default:
                        break;
                }
            }

            if ( firstQ20Index != -1 ) {
                int sequenceLength = curIndex - firstQ20Index;
                if ( sequenceLength > MIN_SEQUENCE_LENGTH ) {
                    byte[] sequence = new byte[sequenceLength];
                    System.arraycopy(sequencedReadBases, firstQ20Index, sequence, 0, sequenceLength);
                    sequences.add(sequence);
                }
            }
        }

        return sequences;
    }

    private void createDeBruijnGraph(List<byte[]> reads) {

        // create the graph
        createGraphFromSequences(reads);

        // remove nodes with incoming multiplicity of N
        // if ( MIN_MULTIPLICITY_TO_USE > 0 )
        //     removeNodesWithLowMultiplicity();

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
        int additionalSequenceFromV2 = V2.actualSequence.length - KMER_LENGTH + 1;
        byte[] newKmer = new byte[V1.actualSequence.length + additionalSequenceFromV2];
        System.arraycopy(V1.actualSequence, 0, newKmer, 0, V1.actualSequence.length);
        System.arraycopy(V2.actualSequence, KMER_LENGTH - 1, newKmer, V1.actualSequence.length, additionalSequenceFromV2);
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

        // remove the first k-1 bases of the kmers
        for ( DeBruijnVertex v :  graph.vertexSet() ) {
            if ( graph.inDegreeOf(v) > 0 )
                v.removePrefix(KMER_LENGTH - 1, true);
        }

        // move common suffixes from incoming nodes to this one

        while ( true ) {

            boolean graphWasModified = false;
            for ( DeBruijnVertex v :  graph.vertexSet() ) {

                if ( graph.inDegreeOf(v) > 1 )  {
                    Set<DeBruijnVertex> connectedVs = new HashSet<DeBruijnVertex>();
                    for ( DeBruijnEdge edge : graph.incomingEdgesOf(v) )
                        connectedVs.add(graph.getEdgeSource(edge));

                    if ( propagateCommonSuffix(v, connectedVs) ) {
                        removeEmptyNodes();
                        graphWasModified = true;
                        break;
                    }
                }
            }

            if ( !graphWasModified )
                break;
        }
    }

    private void removeEmptyNodes() {

        // remember that results of an iteration on a set are undefined when the graph is modified
        while ( true ) {

            boolean graphWasModified = false;
            for ( DeBruijnVertex v :  graph.vertexSet() ) {
                if ( v.printableSequence.length == 0 ) {
                    removeNode(v);
                    graphWasModified = true;
                    break;
                }
            }

            if ( !graphWasModified )
                break;
        }
    }

    private void removeNode(DeBruijnVertex v) {
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
    }

    private boolean propagateCommonSuffix(DeBruijnVertex Vx, Set<DeBruijnVertex> incoming) {

        // find the common matching suffix
        byte[] match = null;
        for ( DeBruijnVertex v : incoming ) {
            if ( match == null ) {
                match = v.printableSequence;
            } else {
                int idx = 0;
                while ( idx < match.length && idx < v.printableSequence.length && match[match.length - idx - 1] == v.printableSequence[v.printableSequence.length - idx - 1] )
                    idx++;

                if ( idx < match.length ) {
                    match = new byte[idx];
                    System.arraycopy(v.printableSequence, v.printableSequence.length - idx, match, 0, idx);
                }
            }
        }

        // if there is a common suffix...
        if ( match != null && match.length > 0 ) {

            // remove the suffix from the end of the incoming nodes...
            for ( DeBruijnVertex v : incoming )
                v.removeSuffix(match.length, false);

            // ...and put it at the front of this node
            Vx.addPrefix(match, false);
            return true;
        }

        return false;
    }

    private void printGraph() {

        for ( DeBruijnVertex source : graph.vertexSet() ) {
            if ( graph.inDegreeOf(source) == 0 )
                getOutputStream().print("* ");
            getOutputStream().print(source + " -> ");
            for ( DeBruijnEdge edge : graph.outgoingEdgesOf(source) ) {
                getOutputStream().print(graph.getEdgeTarget(edge) + " (" + edge.getMultiplicity() + "), ");
            }
            getOutputStream().println();
        }
        getOutputStream().println("------------\n");
    }

    private void findBestPaths() {

        // find them
        List<KBestPaths.Path> bestPaths = KBestPaths.getKBestPaths(graph, 2);

        // print them out
        for ( KBestPaths.Path path : bestPaths ) {

            List<DeBruijnEdge> edges = path.getEdges();
            for (int i = 0; i < edges.size(); i++) {

                DeBruijnEdge edge = edges.get(i);

                if ( i == 0 )
                    getOutputStream().print(graph.getEdgeSource(edge));

                getOutputStream().print(graph.getEdgeTarget(edge));
            }

            if ( edges.size() == 0 )
                getOutputStream().print(path.getLastVertexInPath());

            getOutputStream().println(" (score=" + path.getScore() + ", lowestEdge=" + path.getLowestEdge() +")");
        }
    }

    private void assignReadsToGraph(List<byte[]> reads) {

        // TODO -- implement me

    }

    /****
    private void removeNodesWithLowMultiplicity() {

        Set<DeBruijnVertex> vertexSet = graph.vertexSet();
        // convert to array because results of the iteration on a set are undefined when the graph is modified
        ArrayList<DeBruijnVertex> vertices = new ArrayList<DeBruijnVertex>(vertexSet);

        for (int i = 0; i < vertices.size(); i++) {

            DeBruijnVertex v = vertices.get(i);
            if ( graph.inDegreeOf(v) == 1 &&
                    graph.incomingEdgesOf(v).iterator().next().getMultiplicity() < MIN_MULTIPLICITY_TO_USE )
                removeNode(v);
        }
    }
    ****/
}

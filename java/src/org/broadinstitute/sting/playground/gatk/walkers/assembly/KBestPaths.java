package org.broadinstitute.sting.playground.gatk.walkers.assembly;

import org.jgrapht.graph.DefaultDirectedGraph;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */
// Class for finding the K best paths (as determined by the sum of multiplicities of the edges) in a graph.
// This is different from most graph traversals because we want to test paths from any source node to any sink node.
public class KBestPaths {

    // static access only
    protected KBestPaths() { }

    // class to keep track of paths
    protected static class Path {

        // the last vertex seen in the path
        private DeBruijnVertex lastVertex;

        // the list of edges comprising the path
        private List<DeBruijnEdge> edges;

        // the scores for the path
        private int totalScore = 0, lowestEdge = -1;

        public Path(DeBruijnVertex initialVertex) {
            lastVertex = initialVertex;
            edges = new ArrayList<DeBruijnEdge>(0);
        }

        public Path(Path p, DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, DeBruijnEdge edge) {
            lastVertex = graph.getEdgeTarget(edge);
            edges = new ArrayList<DeBruijnEdge>(p.edges);
            edges.add(edge);
            totalScore = p.totalScore + edge.getMultiplicity();
            lowestEdge = ( p.lowestEdge == -1 ) ? edge.getMultiplicity() : Math.min(p.lowestEdge, edge.getMultiplicity());
        }

        public boolean containsEdge(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, DeBruijnEdge edge) {
            for ( DeBruijnEdge e : edges ) {
                if ( e.equals(graph, edge))
                    return true;
            }

            return false;
        }

        public List<DeBruijnEdge> getEdges() { return edges; }

        public int getScore() { return totalScore; }

        public int getLowestEdge() { return lowestEdge; }

        public DeBruijnVertex getLastVertexInPath() { return lastVertex; }
    }

    protected static class PathComparator implements Comparator<Path> {
        public int compare(final Path path1, final Path path2) {
            return path1.totalScore - path2.totalScore;
        }
    }

    public static List<Path> getKBestPaths(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, int k) {
        PriorityQueue<Path> bestPaths = new PriorityQueue<Path>(k, new PathComparator());

        // run a DFS for best paths
        for ( DeBruijnVertex v : graph.vertexSet() ) {
            if ( graph.inDegreeOf(v) == 0 )
                findBestPaths(graph, new Path(v), k, bestPaths);
        }

        return new ArrayList<Path>(bestPaths);
    }

    private static void findBestPaths(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, Path path, int k, PriorityQueue<Path> bestPaths) {

        // did we hit the end of a path?
        if ( graph.outDegreeOf(path.lastVertex) == 0 ) {
            if ( bestPaths.size() < k ) {
                bestPaths.add(path);
            } else if ( bestPaths.peek().totalScore < path.totalScore ) {
                bestPaths.remove();
                bestPaths.add(path);
            }

            return;
        }

        // recursively run DFS
        for ( DeBruijnEdge edge : graph.outgoingEdgesOf(path.lastVertex) ) {

            // make sure the edge is not already in the path
            if ( path.containsEdge(graph, edge) )
                continue;

            Path newPath = new Path(path, graph, edge);
            findBestPaths(graph, newPath, k, bestPaths);
        }
    }

}

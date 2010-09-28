package org.broadinstitute.sting.playground.gatk.walkers.phasing;

import java.util.*;

/*
 * Copyright (c) 2010, The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// Represents an undirected graph with no self-edges:
public class Graph implements Iterable<GraphEdge> {
    private Neighbors[] adj;

    public Graph(int numVertices) {
        adj = new Neighbors[numVertices];
        for (int i = 0; i < numVertices; i++)
            adj[i] = new Neighbors();
    }

    public void addEdge(GraphEdge e) {
        if (e.v1 == e.v2) // do not permit self-edges
            return;

        adj[e.v1].addNeighbor(e);
        adj[e.v2].addNeighbor(e);
    }

    public void addEdges(Collection<GraphEdge> edges) {
        for (GraphEdge e : edges)
            addEdge(e);
    }

    public void removeEdge(GraphEdge e) {
        adj[e.v1].removeNeighbor(e);
        adj[e.v2].removeNeighbor(e);
    }

    public Collection<GraphEdge> removeAllIncidentEdges(int vertexIndex) {
        Collection<GraphEdge> incidentEdges = new TreeSet<GraphEdge>(adj[vertexIndex].neighbors); // implemented GraphEdge.compareTo()

        for (GraphEdge neighbEdge : incidentEdges) {
            if (vertexIndex != neighbEdge.v1) // vertexIndex == neighbEdge.v2
                adj[neighbEdge.v1].removeNeighbor(neighbEdge);
            else if (vertexIndex != neighbEdge.v2) // vertexIndex == neighbEdge.v1
                adj[neighbEdge.v2].removeNeighbor(neighbEdge);
        }
        adj[vertexIndex].clearAllNeighbors();

        return incidentEdges;
    }

    public DisjointSet getConnectedComponents() {
        DisjointSet cc = new DisjointSet(adj.length);

        for (GraphEdge e : this)
            cc.setUnion(e.v1, e.v2);

        return cc;
    }

    public Iterator<GraphEdge> iterator() {
        return new AllEdgesIterator();
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < adj.length; i++) {
            sb.append(i + ":");
            for (GraphEdge e : adj[i]) {
                sb.append(" " + (e.v1 == i ? e.v2 : e.v1));
            }
            sb.append("\n");
        }

        return sb.toString();
    }

    private class AllEdgesIterator implements Iterator<GraphEdge> {
        private int curInd;
        private Iterator<GraphEdge> innerIt;
        private GraphEdge nextEdge;

        public AllEdgesIterator() {
            curInd = 0;
            innerIt = null;
            nextEdge = null;
        }

        public boolean hasNext() {
            if (nextEdge != null)
                return true;

            for (; curInd < adj.length; curInd++) {
                if (innerIt == null)
                    innerIt = adj[curInd].iterator();

                while (innerIt.hasNext()) {
                    GraphEdge e = innerIt.next();
                    if (e.v1 == curInd) { // only want to see each edge once
                        nextEdge = e;
                        return true;
                    }
                }

                innerIt = null;
            }

            return false;
        }

        public GraphEdge next() {
            if (!hasNext())
                throw new NoSuchElementException();

            GraphEdge tmpEdge = nextEdge;
            nextEdge = null;
            return tmpEdge;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    private class Neighbors implements Iterable<GraphEdge> {
        private Set<GraphEdge> neighbors;

        public Neighbors() {
            this.neighbors = new TreeSet<GraphEdge>(); // implemented GraphEdge.compareTo()
        }

        public void addNeighbor(GraphEdge e) {
            neighbors.add(e);
        }

        public void removeNeighbor(GraphEdge e) {
            neighbors.remove(e);
        }

        public Iterator<GraphEdge> iterator() {
            return neighbors.iterator();
        }

        public void clearAllNeighbors() {
            neighbors.clear();
        }
    }
}

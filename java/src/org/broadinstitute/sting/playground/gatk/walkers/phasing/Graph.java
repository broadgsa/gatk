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

public class Graph implements Iterable<GraphEdge> {
    private Neighbors[] adj;

    public Graph(int numVertices) {
        adj = new Neighbors[numVertices];
        for (int i = 0; i < numVertices; i++)
            adj[i] = new Neighbors();
    }

    public void addEdge(GraphEdge e) {
        adj[e.v1].addNeighbor(e);
        adj[e.v2].addNeighbor(e);
    }

    public void removeEdge(GraphEdge e) {
        adj[e.v1].removeNeighbor(e);
        adj[e.v2].removeNeighbor(e);
    }

    public DisjointSet getConnectedComponents() {
        DisjointSet cc = new DisjointSet(adj.length);

        for (int i = 0; i < adj.length; i++)
            for (GraphEdge e : adj[i])
                cc.setUnion(e.v1, e.v2);

        return cc;
    }

    // Note that this will give each edge TWICE [since e=(v1,v2) is stored as a neighbor for both v1 and v2]
    public Iterator<GraphEdge> iterator() {
        return new AllEdgesIterator();
    }

    public List<GraphEdge> getAllEdges() {
        Set<GraphEdge> allEdges = new TreeSet<GraphEdge>(); // implemented GraphEdge.compareTo()
        for (GraphEdge e : this)
            allEdges.add(e);

        return new LinkedList<GraphEdge>(allEdges);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < adj.length; i++) {
            sb.append(i + ":");
            for (GraphEdge e : adj[i])
                sb.append(" " + e);
            sb.append("\n");
        }

        return sb.toString();
    }

    private class AllEdgesIterator implements Iterator<GraphEdge> {
        private int curInd;
        private Iterator<GraphEdge> innerIt;

        public AllEdgesIterator() {
            curInd = 0;
            innerIt = null;
        }

        public boolean hasNext() {
            for (; curInd < adj.length; curInd++) {
                if (innerIt == null)
                    innerIt = adj[curInd].iterator();
                if (innerIt.hasNext())
                    return true;

                innerIt = null;
            }

            return false;
        }

        public GraphEdge next() {
            if (!hasNext())
                throw new NoSuchElementException();

            return innerIt.next();
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
    }
}

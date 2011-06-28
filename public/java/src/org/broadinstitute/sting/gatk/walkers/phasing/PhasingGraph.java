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
package org.broadinstitute.sting.gatk.walkers.phasing;

import org.broadinstitute.sting.utils.DisjointSet;

import java.util.*;

// Represents an undirected graph with no self-edges:
public class PhasingGraph implements Iterable<PhasingGraphEdge> {
    private Neighbors[] adj;

    public PhasingGraph(int numVertices) {
        adj = new Neighbors[numVertices];
        for (int i = 0; i < numVertices; i++)
            adj[i] = new Neighbors();
    }

    public void addEdge(PhasingGraphEdge e) {
        if (e.v1 == e.v2) // do not permit self-edges
            return;

        adj[e.v1].addNeighbor(e);
        adj[e.v2].addNeighbor(e);
    }

    public void addEdges(Collection<PhasingGraphEdge> edges) {
        for (PhasingGraphEdge e : edges)
            addEdge(e);
    }

    public void removeEdge(PhasingGraphEdge e) {
        adj[e.v1].removeNeighbor(e);
        adj[e.v2].removeNeighbor(e);
    }

    public Collection<PhasingGraphEdge> removeAllIncidentEdges(int vertexIndex) {
        Collection<PhasingGraphEdge> incidentEdges = new TreeSet<PhasingGraphEdge>(adj[vertexIndex].neighbors); // implemented GraphEdge.compareTo()

        for (PhasingGraphEdge neighbEdge : incidentEdges) {
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

        for (PhasingGraphEdge e : this)
            cc.setUnion(e.v1, e.v2);

        return cc;
    }

    public Iterator<PhasingGraphEdge> iterator() {
        return new AllEdgesIterator();
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < adj.length; i++) {
            sb.append(i + ":");
            for (PhasingGraphEdge e : adj[i]) {
                sb.append(" " + (e.v1 == i ? e.v2 : e.v1));
            }
            sb.append("\n");
        }

        return sb.toString();
    }

    private class AllEdgesIterator implements Iterator<PhasingGraphEdge> {
        private int curInd;
        private Iterator<PhasingGraphEdge> innerIt;
        private PhasingGraphEdge nextEdge;

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
                    PhasingGraphEdge e = innerIt.next();
                    if (e.v1 == curInd) { // only want to see each edge once
                        nextEdge = e;
                        return true;
                    }
                }

                innerIt = null;
            }

            return false;
        }

        public PhasingGraphEdge next() {
            if (!hasNext())
                throw new NoSuchElementException();

            PhasingGraphEdge tmpEdge = nextEdge;
            nextEdge = null;
            return tmpEdge;
        }

        public void remove() {
            throw new UnsupportedOperationException();
        }
    }

    private class Neighbors implements Iterable<PhasingGraphEdge> {
        private Set<PhasingGraphEdge> neighbors;

        public Neighbors() {
            this.neighbors = new TreeSet<PhasingGraphEdge>(); // implemented GraphEdge.compareTo()
        }

        public void addNeighbor(PhasingGraphEdge e) {
            neighbors.add(e);
        }

        public void removeNeighbor(PhasingGraphEdge e) {
            neighbors.remove(e);
        }

        public Iterator<PhasingGraphEdge> iterator() {
            return neighbors.iterator();
        }

        public void clearAllNeighbors() {
            neighbors.clear();
        }
    }
}

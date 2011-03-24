package org.broadinstitute.sting.playground.gatk.walkers.assembly;

import org.jgrapht.graph.DefaultDirectedGraph;

/**
 * Created by IntelliJ IDEA.
 * User: ebanks
 * Date: Mar 23, 2011
 */
// simple edge class for connecting nodes in the graph
public class DeBruijnEdge {

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

    public boolean equals(DefaultDirectedGraph<DeBruijnVertex, DeBruijnEdge> graph, DeBruijnEdge edge) {
        return (graph.getEdgeSource(this) == graph.getEdgeSource(edge)) && (graph.getEdgeTarget(this) == graph.getEdgeTarget(edge));
    }
}

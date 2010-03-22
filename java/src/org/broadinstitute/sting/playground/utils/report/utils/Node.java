package org.broadinstitute.sting.playground.utils.report.utils;

import java.util.Collection;
import java.util.LinkedHashSet;

/**
* a node, extracted using the new report output system.
*/
public class Node {
    public String value;
    public String description; // if one is available 
    public Collection<Node> children;

    public Node(String value) {
        this.value = value;
    }

    public void addChild(Node child) {
        if (children == null) children = new LinkedHashSet<Node>();
        children.add(child);
    }

    public void addAllChildren(Collection<Node> children) {
        if (this.children == null) this.children = new LinkedHashSet<Node>();
        this.children.addAll(children);
    }

    public Boolean getComplex() { return (children != null && children.size() > 0); }

    /**
     * a convenience method for adding a new sub-node with the specified value
     * @param value the value of the sub-node
     */
    public void createSubNode(String value) {
        addChild(new Node(value));
    }

    public String getValue() {
        return value;
    }

    public String getDescription() {
        return description;
    }

    public Collection<Node> getChildren() {
        return children;
    }
}

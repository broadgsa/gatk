package org.broadinstitute.sting.playground.utils.report.utils;

import java.util.Collection;
import java.util.LinkedHashSet;

/**
 * a node, extracted using the new report output system.
 */
public class Node {
    public final String name;
    public final String value;
    public final String description;
    public final boolean display; // is this node an output node, or a node for tracking internal data? true if output node
    public Collection<Node> children;

    public Node(String name, String value, String description) {
        this.value = value;
        this.name = name;
        this.description = description;
        display = true;
    }

    public Node(String name, String value, String description, boolean display) {
        this.value = value;
        this.name = name;
        this.description = description;
        this.display = display;
    }

    public void addChild(Node child) {
        if (children == null) children = new LinkedHashSet<Node>();
        children.add(child);
    }

    public void addAllChildren(Collection<Node> children) {
        if (this.children == null) this.children = new LinkedHashSet<Node>();
        this.children.addAll(children);
    }

    public Boolean getComplex() {
        return (children != null && children.size() > 0);
    }

    /**
     * a convenience method for adding a new sub-node with the specified value
     *
     * @param value the value of the sub-node
     */
    public void createSubNode(String name, String value, String description) {
        addChild(new Node(name, value, description));
    }

    public String getValue() {
        return value;
    }

    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public Collection<Node> getChildren() {
        return children;
    }

    public boolean getDisplay() {
        return display;
    }
}

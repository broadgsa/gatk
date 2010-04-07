package org.broadinstitute.sting.playground.utils.report.utils;

import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * a node, extracted using the new report output system.
 */
public class Node {
    public String name;
    public String value;
    public String description;
    public boolean display; // is this node an output node, or a node for tracking internal data? true if output node
    public boolean table;   // this is a hack, but I needed a way to indicate that a node was a row root node
    public boolean tag;
    public Collection<Node> children;

    public Node(String name, String value, String description) {
        this.value = value;
        this.name = name;
        this.description = description;
        display = true;
        table = false;
        tag = false;
    }

    public Node(String name, String value, String description, boolean display) {
        this.value = value;
        this.name = name;
        this.description = description;
        this.display = display;
        table = false;
        tag = false;
    }

    public void setTable() {table = true;}

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
         return (children == null) ?  new ArrayList<Node>() : children;
    }

    public boolean getDisplay() {
        return display;
    }

    public boolean getTable() {
        return table;
    }

    public boolean getTag() {
        return tag;
    }

    public void setTag() {
        this.tag = true;
    }

    public void clone(Node n) {
        this.name = n.name;
        this.value = n.value;
        this.description = n.description;
        this.display = n.display;
        this.table = n.table;
        this.tag = n.tag;
        this.children = new LinkedHashSet<Node>();
        if (n.children != null) this.children.addAll(n.getChildren());
    }

    public List<List<Node>> getTableRows() {
        List<List<Node>> ret =  NodeUtils.flattenToRow(this,false);
        return ret;
    }

    public List<List<Node>> getTableRowsNoTables() {
        List<List<Node>> ret =  NodeUtils.flattenToRow(this,true);
        return ret;
    }

    public int getRowCount() {
        return NodeUtils.flattenToRowCount(this);
    }
}

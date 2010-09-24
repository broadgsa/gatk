package org.broadinstitute.sting.utils.report.utils;

import java.util.ArrayList;
import java.util.List;


/**
 * 
 * @author aaron 
 * 
 * Class NodeUtils
 *
 * utilities for working with nodes
 */
public class NodeUtils {

    static class NodeMarker {
        private Node node;
        
        public NodeMarker(Node n) {
            node = n;
        }

        public int rowCount() {
            int sum = (node.table) ? node.getChildren().size() : 1;
            for (Node n : node.getChildren()) {
                NodeMarker fn = new NodeMarker(n);
                sum = sum * fn.rowCount();
            }
            return sum;
        }

        private boolean validLeafNode() {
            return node.getChildren().size() == 0 && node.display && !node.tag;
        }

        private List<List<Node>> addToEachList(List<List<Node>> list) {
            for (List<Node> lt : list)
                    lt.add(node);
            return list;
        }

        public List<List<Node>> toRow(List<List<Node>> oldList, boolean excludeTables) {
            // if we're a leaf node that isn't a tag, add it to each list
            if (validLeafNode())
                addToEachList(oldList);

            // special case: if we've just got a single node, traverse into it
            else if (node.getChildren().size() > 0 && !node.table)
                for (Node n : node.children) {
                    oldList = new NodeMarker(n).toRow(oldList, excludeTables);
                }
            // when we encounter a table we want to branch into multiple rows
            else if (node.table && !excludeTables) {
                List<List<Node>> newList = new ArrayList<List<Node>>();
                for (Node child : node.children) {
                    if (child.display && !child.tag) {
                        List<List<Node>> tempList = new ArrayList<List<Node>>();
                        tempList.add(new ArrayList<Node>());
                        tempList.get(0).add(child);
                        NodeMarker marker = new NodeMarker(child);
                        List<List<Node>> carry = marker.toRow(tempList, excludeTables);
                        newList.addAll(carry);
                    }
                }
                List<List<Node>> ret = new ArrayList<List<Node>>();
                // permutations of each previous list and the new temp list
                for (List<Node> original : oldList)
                    for (List<Node> lst : newList) {
                        List<Node> temp = new ArrayList<Node>();
                        temp.addAll(original);
                        temp.addAll(lst);
                        ret.add(temp);
                    }
                return ret;
            }
            // be default return the old list
            return oldList;
        }

    }


    // given a node, get the number of rows it will generate
    public static int flattenToRowCount(Node n) {
        NodeMarker fn = new NodeMarker(n);
        return fn.rowCount();
    }

    // given a node, generate rows (flattening tables)
    public static List<List<Node>> flattenToRow(Node n, boolean excludeTables) {
        NodeMarker fn = new NodeMarker(n);
        List<List<Node>> nodesList = new ArrayList<List<Node>>();
        nodesList.add(new ArrayList<Node>());
        return fn.toRow(nodesList, excludeTables);
    }
}

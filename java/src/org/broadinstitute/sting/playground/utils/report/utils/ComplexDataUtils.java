package org.broadinstitute.sting.playground.utils.report.utils;

import java.util.*;


/**
 * 
 * @author aaron 
 * 
 * Class ComplexDataUtils
 *
 * This class contains methods and techniques for breaking down complex data in the output system
 */
public class ComplexDataUtils {

    /**
     * convert any string -> object pairing into a string keyed tree
     * @param obj the object
     * @return a mapping of the name to the associated value tree.  All non-leaf nodes will be Strings
     */
    public static Collection<Node> resolveObjects(Object obj) { // TODO: fix this, we need a way to get the name of the list from the data point
        Collection<Node> nodes = new ArrayList<Node>();

        // the simplest case
        if (obj == null) nodes.add(new Node("<null>")); // throw new IllegalStateException("object is null");
        else if (obj instanceof TableType)
            nodes.add(tableToNode((TableType)obj, ((TableType)obj).getName()));
        else if (obj instanceof Map) {
            for (Object key : ((Map)obj).keySet()) {
                Node keyNode = new Node(key.toString());
                nodes.add(keyNode);
                keyNode.addAllChildren(resolveObjects(((Map)obj).get(key)));
            }
        } else if (obj instanceof Collection)
            nodes.addAll(listToNode((Collection)obj, "collection"));
        else if (obj.getClass().isArray())
            nodes.addAll(listToNode(Arrays.asList(obj),"array"));
        else
            nodes.add(new Node(obj.toString()));
        return nodes;
    }


    /**
     * given a TableType object, make it into a tree using maps.
     * @param table the table type to convert into a map
     * @return
     */
    private static Node tableToNode(TableType table, String name) {
        Node root = new Node(name);
        Object[] rows = table.getRowKeys();
        Object[] cols = table.getColumnKeys();
        for (int x = 0; x < table.getRowKeys().length; x++) {
            Node row = new Node(rows[x].toString());
            root.addChild(row);
            for (int y = 0; y < table.getRowKeys().length; y++) {
                Node col = new Node(cols[y].toString());
                row.addChild(col);
                col.addChild(new Node(table.getCell(x,y)));
            }
        }
        return root;
    }

    /**
     * given a Collection object, make it into a tree using maps.
     * @param coll the collection to iterate, and turn into a list
     * @return a mapping of String to Object
     */
    private static Collection<Node> listToNode(Collection coll, String name) {
        Collection<Node> nodes = new ArrayList<Node>();
        Iterator<Object> iter = coll.iterator();
        for (int x = 0; x < coll.size(); x++) {
            Node value = new Node(String.valueOf(x));
            value.addChild(new Node(iter.next().toString()));
            nodes.add(value);
        }
        return nodes;
    }

}

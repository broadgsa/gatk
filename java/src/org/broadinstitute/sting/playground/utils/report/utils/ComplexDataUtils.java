package org.broadinstitute.sting.playground.utils.report.utils;

import org.broadinstitute.sting.utils.Utils;

import java.util.*;


/**
 * @author aaron
 *         <p/>
 *         Class ComplexDataUtils
 *         <p/>
 *         This class contains methods and techniques for breaking down complex data in the output system
 */
public class ComplexDataUtils {

    /**
     * convert any string -> object pairing into a string keyed tree
     *
     * @param obj the object
     * @return a mapping of the name to the associated value tree.  All non-leaf nodes will be Strings
     */
    public static Collection<Node> resolveObjects(Object obj) { // TODO: fix this, we need a way to get the name of the list from the data point
        Collection<Node> nodes = new ArrayList<Node>();

        // the simplest case
        if (obj == null)
            nodes.add(new Node("<null>", "<null>", "<null>")); // throw new IllegalStateException("object is null");
        else if (obj instanceof TableType)
            nodes.add(tableToNode((TableType) obj, ((TableType) obj).getName()));
        else if (obj instanceof Map) {
            for (Object key : ((Map) obj).keySet()) {
                Node keyNode = new Node("key", key.toString(), "map key");
                nodes.add(keyNode);
                keyNode.addAllChildren(resolveObjects(((Map) obj).get(key)));
            }
        } else if (obj instanceof Collection)
            nodes.addAll(listToNode((Collection) obj, "collection"));
        else if (obj.getClass().isArray())
            nodes.addAll(listToNode(Arrays.asList(obj), "array"));
        else
            nodes.add(new Node(obj.getClass().getSimpleName(), obj.toString(), "value"));
        return nodes;
    }


    /**
     * given a TableType object, make it into a tree using maps.
     *
     * @param table the table type to convert into a map
     * @return
     */
    private static Node tableToNode(TableType table, String name) {
        Node root = new Node("table", name, "Table");
        Object[] rows = table.getRowKeys();
        Object[] cols = table.getColumnKeys();

        // add the columns names
        Node col_names = new Node("col_names", "col_names", "the column names, ~ seperated", false);
        for (Object col : cols)
            col_names.addChild(new Node("col", col.toString(), "a column name", false));
        root.addChild(col_names);

        for (int x = 0; x < table.getRowKeys().length; x++) {
            Node row = new Node("row", rows[x].toString(), "a row");
            root.addChild(row);
            for (int y = 0; y < table.getColumnKeys().length; y++) {
                Node col = new Node("column", cols[y].toString(), "a column");
                row.addChild(col);
                col.addChild(new Node("cell(" + x + "," + y + ")", table.getCell(x, y), "a cell"));
            }
        }
        return root;
    }

    /**
     * given a Collection object, make it into a tree using maps.
     *
     * @param coll the collection to iterate, and turn into a list
     * @return a mapping of String to Object
     */
    private static Collection<Node> listToNode(Collection coll, String name) {
        Collection<Node> nodes = new ArrayList<Node>();
        Iterator<Object> iter = coll.iterator();
        for (int x = 0; x < coll.size(); x++) {
            Node value = new Node("column " + x, String.valueOf(x), "column");
            value.addChild(new Node("value " + x, iter.next().toString(), "value"));
            nodes.add(value);
        }
        return nodes;
    }
}

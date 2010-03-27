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

        // the simplest case, the object is null
        if (obj == null) nodes.add(new Node("<null>", "<null>", "<null>"));
        // capture objects of type TableTable
        else if (obj instanceof TableType)
            nodes.add(tableToNode((TableType) obj, ((TableType) obj).getName()));

        // try to handle maps
        else if (obj instanceof Map) {
            extractMap(obj, nodes);

        // handle collections
        } else if (obj instanceof Collection)
            nodes.addAll(listToNode((Collection) obj, "collection"));

        // arrays
        else if (obj.getClass().isArray())
            nodes.addAll(listToNode(Arrays.asList(obj), "array"));

        // else we have a simple object (at least try to handle it that way
        else
            nodes.add(extractPlainObjectOrPrimitive(obj.getClass().getSimpleName(),obj));

        // return the collection of nodes we've parsed out
        return nodes;
    }

    /**
     * extract a map object
     * @param obj the object (instance of Map)
     * @param nodes the node list to add our key->values to
     */
    private static void extractMap(Object obj, Collection<Node> nodes) {
        for (Object key : ((Map) obj).keySet()) {
            Node keyNode = new Node("key", key.toString(), "map key");
            nodes.add(keyNode);
            keyNode.addAllChildren(resolveObjects(((Map) obj).get(key)));
        }
        // special case: if the map is empty, add a null node
        if (nodes.isEmpty()) nodes.add(new Node("<null>", "<null>", "<null>"));
    }

    /**
     * extract a (hopefully) primitive value
     * @param obj the object
     */
    private static Node extractPlainObjectOrPrimitive(String name, Object obj) {
        String value = "<null>";
        if (obj instanceof Float || obj instanceof Double)
            value = String.format("%.4f",(Double)obj);
        else
            value = obj.toString();
        return new Node(name, value, "value");
    }

    /**
     * given a TableType object, make it into a tree using maps.
     *
     * @param table the table type to convert into a map
     * @return a node representing this table
     */
    private static Node tableToNode(TableType table, String name) {
        Node root = new Node("table", name, "Table");
        root.setTable();
        Object[] rows = table.getRowKeys();
        Object[] cols = table.getColumnKeys();

        // add the columns names
        for (int x = 0; x < table.getRowKeys().length; x++) {
            Node row = new Node("row", rows[x].toString(), "a row in a table");
            root.addChild(row);
            for (int y = 0; y < table.getColumnKeys().length; y++) {
                Node col = new Node("column", cols[y].toString(), "columns in a table");
                row.addChild(col);
                col.addChild(extractPlainObjectOrPrimitive("cell(" + x + "," + y + ")", table.getCell(x, y)));
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

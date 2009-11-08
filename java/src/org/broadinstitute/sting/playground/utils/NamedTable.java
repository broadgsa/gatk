package org.broadinstitute.sting.playground.utils;

import java.util.*;

/**
 * NamedTable is a utility class for maintaining a table and accessing rows and columns
 * with named identifiers, rather than indicies that must be remembered.  It also grows
 * dynamically; you needn't specify the rows and columns before accessing them, so you can
 * continuously expand the table in situations where you don't necessarily know how many
 * rows or columns you'll have in the end.
 */
public class NamedTable {
    // If in the future, this class gets templatized, the counter variable should really become a CountedObject<T>.
    private HashMap<String, HashMap<String, Double>> table;
    private HashSet<String> rowNames;
    private HashSet<String> colNames;

    public NamedTable() {
        table = new HashMap<String, HashMap<String, Double>>();
        rowNames = new HashSet<String>();
        colNames = new HashSet<String>();
    }

    /**
     * Copy another table into this new table
     * @param ct  the table to copy
     */
    public NamedTable(NamedTable ct) {
        table = new HashMap<String, HashMap<String, Double>>();
        rowNames = new HashSet<String>();
        colNames = new HashSet<String>();

        for (String rowName : ct.getRowNames()) {
            for (String colName : ct.getColumnNames()) {
                this.set(rowName, colName, ct.get(rowName, colName));
            }
        }
    }

    /**
     * If the entry we're trying to access doesn't exist, create it.
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     */
    private void verifyEntry(String rowName, String colName) {
        rowNames.add(rowName);
        colNames.add(colName);

        if (!table.containsKey(rowName)) {
            table.put(rowName, new HashMap<String, Double>());
        }

        if (!table.get(rowName).containsKey(colName)) {
            table.get(rowName).put(colName, 0.0);
        }
    }

    /**
     * Set an entry in the table
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     * @param value    the value to set for the (row,column)-th entry
     */
    public void set(String rowName, String colName, double value) {
        verifyEntry(rowName, colName);

        table.get(rowName).put(colName, value);
    }

    /**
     * Get an entry in the table
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     * @return the value of the (row,column)-th entry
     */
    public double get(String rowName, String colName) {
        verifyEntry(rowName, colName);

        return table.get(rowName).get(colName);
    }

    /**
     * For convenience, increment the (row,column)-th entry in the table so that the
     * user doesn't need to extract, increment, and then reassign the value.  One day,
     * this should probably be rewritten to use Andrey's CountedObject class.
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     */
    public void increment(String rowName, String colName) {
        double value = get(rowName, colName);
        
        table.get(rowName).put(colName, value + 1.0);
    }

    /**
     * For convenience, decrement the (row,column)-th entry in the table so that the
     * user doesn't need to extract, increment, and then reassign the value.  One day,
     * this should probably be rewritten to use Andrey's CountedObject class.
     *
     * @param rowName  the name of the row
     * @param colName  the name of the column
     */
    public void decrement(String rowName, String colName) {
        double value = get(rowName, colName);

        table.get(rowName).put(colName, value - 1.0);
    }

    /**
     * Get a sorted list of all the rows in the table
     *
     * @return a sorted list of all the rows in the table
     */
    public ArrayList<String> getRowNames() {
        ArrayList<String> rows = new ArrayList<String>();
        Iterator<String> rowit = rowNames.iterator();
        while (rowit.hasNext()) {
            rows.add(rowit.next());
        }
        Collections.sort(rows);

        return rows;
    }

    /**
     * Get a sorted list of all the columns in the table
     *
     * @return a sorted list of all the columns in the table
     */
    public ArrayList<String> getColumnNames() {
        ArrayList<String> cols = new ArrayList<String>();
        Iterator<String> colit = colNames.iterator();
        while (colit.hasNext()) {
            cols.add(colit.next());
        }
        Collections.sort(cols);

        return cols;
    }

    /**
     * Get a new table representing a subset of the current table
     * @param rowNames  a list of rows to extract
     * @param colNames  a list of columns to extract
     * @return  the subsetted table
     */
    public NamedTable getSubset(ArrayList<String> rowNames, ArrayList<String> colNames) {
        NamedTable ct = new NamedTable();

        for (String rowName : rowNames) {
            for (String colName : colNames) {
                ct.set(rowName, colName, get(rowName, colName));
            }
        }

        return ct;
    }

    /* This stuff doesn't belong in this class, but I don't want
       to delete the code until it's moved somewhere appropriate */
    /*
    public boolean twoTailedFisherExactTest(double pValueLimit) {
        return (twoTailedFisherExactTest() < pValueLimit);
    }

    public double oneTailedFisherExactTestRight() {
        NamedTable ct = new NamedTable(this);

        double pCutoff = pValue();
        double pValue = pCutoff;

        while (ct.rotateRight()) {
            double pValuePiece = ct.pValue();

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        return pValue;
    }

    public double oneTailedFisherExactTestLeft() {
        NamedTable ct = new NamedTable(this);

        double pCutoff = pValue();
        double pValue = pCutoff;

        while (ct.rotateLeft()) {
            double pValuePiece = ct.pValue();

            if (pValuePiece <= pCutoff) {
                pValue += pValuePiece;
            }
        }

        return pValue;
    }

    public double twoTailedFisherExactTest() {
        return oneTailedFisherExactTestLeft() + oneTailedFisherExactTestRight();
    }

    public double pValue() {
        double p = 0.0;
        
        if (rowNames.size() == 2 && colNames.size() == 2) {
            String[] rows = rowNames.toArray(new String[1]);
            String[] columns = colNames.toArray(new String[1]);

            double a = get(rows[0], columns[0]);
            double b = get(rows[0], columns[1]);
            double c = get(rows[1], columns[0]);
            double d = get(rows[1], columns[1]);
            double n = a + b + c + d;

            double p1 = Arithmetic.binomial(a + b, (long) a);
            double p2 = Arithmetic.binomial(c + d, (long) c);
            double pn = Arithmetic.binomial(n, (long) (a + c));

            p = p1*p2/pn;
        }

        return p;
    }

    public boolean rotateRight() {
        String[] rows = rowNames.toArray(new String[1]);
        String[] columns = colNames.toArray(new String[1]);

        decrement(rows[0], columns[0]);
        increment(rows[1], columns[0]);

        increment(rows[0], columns[1]);
        decrement(rows[1], columns[1]);

        return (get(rows[0], columns[0]) >= 0 && get(rows[1], columns[1]) >= 0);
    }

    public boolean rotateLeft() {
        String[] rows = rowNames.toArray(new String[1]);
        String[] columns = colNames.toArray(new String[1]);

        increment(rows[0], columns[0]);
        decrement(rows[1], columns[0]);

        decrement(rows[0], columns[1]);
        increment(rows[1], columns[1]);

        return (get(rows[0], columns[1]) >= 0 && get(rows[1], columns[0]) >= 0);
    }
    */

    /**
     * Get a nicely-formatted version of the contents of the table.
     *
     * @return  a String representing the contents of the table
     */
    public String toString() {
        String tableString = "";
        boolean headerPrinted = false;

        ArrayList<String> rows = getRowNames();
        ArrayList<String> cols = getColumnNames();

        for (String rowName : rows) {
            if (!headerPrinted) {
                tableString += "rowName  ";
                for (String colName : cols) {
                    tableString += "\t" + colName;
                }
                tableString += "\n";

                headerPrinted = true;
            }

            tableString += rowName;

            for (String colName : cols) {
                tableString += String.format("\t%7.7f", get(rowName, colName));
            }

            tableString += "\n";
        }

        return tableString;
    }
}

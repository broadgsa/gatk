package org.broadinstitute.sting.playground.utils.report.utils;


/**
 * 
 * @author aaron 
 * 
 * Class TableType
 *
 * an interface for turning arbritary objects into tables
 */
public interface TableType {
    public Object[] getRowKeys();
    public Object[] getColumnKeys();
    public String getCell(int x, int y);
    public String getName();
}

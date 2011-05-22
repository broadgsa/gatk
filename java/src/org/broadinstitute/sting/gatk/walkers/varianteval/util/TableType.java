package org.broadinstitute.sting.gatk.walkers.varianteval.util;


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
    public Object getCell(int x, int y);
    public String getName();
}

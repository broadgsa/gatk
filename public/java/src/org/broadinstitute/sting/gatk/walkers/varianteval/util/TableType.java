package org.broadinstitute.sting.gatk.walkers.varianteval.util;


/**
 * 
 * @author aaron 
 * 
 * Class TableType
 *
 * an interface for turning arbritary objects into tables
 */
public abstract class TableType {
    public abstract Object[] getRowKeys();
    public abstract Object[] getColumnKeys();
    public abstract Object getCell(int x, int y);
    public String getName() { return getClass().getSimpleName(); }
    public String getRowName() { return "row"; }
    
}

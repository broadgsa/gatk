package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.gatk.refdata.utils.RODRecordList;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 10, 2009
 * Time: 6:10:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class RODRecordListImpl extends AbstractList<ReferenceOrderedDatum> implements Comparable<RODRecordList>, Cloneable, RODRecordList {
    private List<ReferenceOrderedDatum> records;
    private GenomeLoc location = null;
    private String name = null;

    public RODRecordListImpl(String name) {
        records = new ArrayList<ReferenceOrderedDatum>();
        this.name = name;
    }

    /**
     * Fully qualified constructor: instantiates a new ReferenceOrderedDatumRecordList object with specified ReferenceOrderedDatum track name, location on the
     * reference, and list of associated ReferenceOrderedDatums. This is a knee-deep COPY constructor: passed name, loc, and data element
     * objects will be referenced from the created ReferenceOrderedDatumRecordList (so that changing them from outside will affect data
     * in this object), however, the data elements will be copied into a newly
     * allocated list, so that the 'data' collection argument can be modified afterwards without affecting the state
     * of this record list. WARNING: this constructor is (semi-)validating: passed name and location
     * are allowed to be nulls (although it maybe unsafe, use caution), but if they are not nulls, then passed non-null ReferenceOrderedDatum data
     * elements must have same track name, and their locations must overlap with the passed 'location' argument. Null
     * data elements or null 'data' collection argument are allowed as well.
     * @param name
     * @param data
     * @param loc
     */
    public RODRecordListImpl(String name, Collection<ReferenceOrderedDatum> data, GenomeLoc loc) {
        this.records = new ArrayList<ReferenceOrderedDatum>(data==null?0:data.size());
        this.name = name;
        this.location = loc;
        if ( data == null || data.size() == 0 ) return; // empty dataset, nothing to do
        for ( ReferenceOrderedDatum r : data ) {
            records.add(r);
            if ( r == null ) continue;
            if ( ! this.name.equals(r.getName() ) ) {
                throw new StingException("Attempt to add ReferenceOrderedDatum with non-matching name "+r.getName()+" to the track "+name);
            }
            if ( location != null && ! location.overlapsP(r.getLocation()) ) {
                    throw new StingException("Attempt to add ReferenceOrderedDatum that lies outside of specified interval "+location+"; offending ReferenceOrderedDatum:\n"+r.toString());
            }
        }
    }


    public GenomeLoc getLocation() { return location; }
    public String getName() { return name; }
    public List<ReferenceOrderedDatum> getRecords() { return records; }
    public Iterator<ReferenceOrderedDatum> iterator() { return records.iterator() ; }
    public void clear() { records.clear(); }
    public boolean isEmpty() { return records.isEmpty(); }

    public boolean add(ReferenceOrderedDatum record) { add(record, false); return true;}

    @Override
    public ReferenceOrderedDatum get(int i) {
        return records.get(i);
    }

    public void add(ReferenceOrderedDatum record, boolean allowNameMismatch) {
        if ( record != null ) {
            if ( ! allowNameMismatch && ! name.equals(record.getName() ) )
                throw new StingException("Attempt to add ReferenceOrderedDatum with non-matching name "+record.getName()+" to the track "+name);
        }
        records.add(record);
    }

    public void add(RODRecordList records ) { add( records, false ); }

    public void add(RODRecordList records, boolean allowNameMismatch) {
        for ( ReferenceOrderedDatum record : records )
            add(record, allowNameMismatch);
    }    

    public int size() { return records.size() ; }

    /**
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object.
     *
     * @param that the object to be compared.
     * @return a negative integer, zero, or a positive integer as this object
     *         is less than, equal to, or greater than the specified object.
     * @throws ClassCastException if the specified object's type prevents it
     *                            from being compared to this object.
     */
    public int compareTo(RODRecordList that) {
        return getLocation().compareTo(that.getLocation());  //To change body of implemented methods use File | Settings | File Templates.
    }
}

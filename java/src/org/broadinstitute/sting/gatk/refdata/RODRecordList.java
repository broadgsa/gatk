package org.broadinstitute.sting.gatk.refdata;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.GenomeLocParser;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 10, 2009
 * Time: 6:10:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class RODRecordList<ROD extends ReferenceOrderedDatum> implements Iterable<ROD>, Comparable<RODRecordList<ROD>>, Cloneable {
    private List<ROD> records;
    private GenomeLoc location = null;
    private String name = null;

    private RODRecordList() {} // dummy constructor for internal use; does not initialize/allocate anything

    public RODRecordList(String name) {
        records = new ArrayList<ROD>();
        this.name = name;
    }

    /**
     * Fully qualified constructor: instantiates a new RODRecordList object with specified ROD track name, location on the
     * reference, and list of associated RODs. This is a knee-deep COPY constructor: passed name, loc, and data element
     * objects will be referenced from the created RODRecordList (so that changing them from outside will affect data
     * in this object), however, the data elements will be copied into a newly
     * allocated list, so that the 'data' collection argument can be modified afterwards without affecting the state
     * of this record list. WARNING: this constructor is (semi-)validating: passed name and location
     * are allowed to be nulls (although it maybe unsafe, use caution), but if they are not nulls, then passed non-null ROD data
     * elements must have same track name, and their locations must overlap with the passed 'location' argument. Null
     * data elements or null 'data' collection argument are allowed as well.
     * @param name
     * @param data
     * @param loc
     */
    public RODRecordList(String name, Collection<ROD> data, GenomeLoc loc) {
        this.records = new ArrayList<ROD>(data==null?0:data.size());
        this.name = name;
        this.location = loc;
        if ( data == null || data.size() == 0 ) return; // empty dataset, nothing to do
        for ( ROD r : data ) {
            records.add(r);
            if ( r == null ) continue;
            if ( ! this.name.equals(r.getName() ) ) {
                throw new StingException("Attempt to add ROD with non-matching name "+r.getName()+" to the track "+name);
            }
            if ( location != null && ! location.overlapsP(r.getLocation()) ) {
                    throw new StingException("Attempt to add ROD that lies outside of specified interval "+location+"; offending ROD:\n"+r.toString());
            }
        }
    }


    public GenomeLoc getLocation() { return location; }
    public void setLocation(GenomeLoc location) { this.location = location; }
    public String getName() { return name; }
    public List<ROD> getRecords() { return records; }
    public Iterator<ROD> iterator() { return records.iterator() ; }
    public void clear() { records.clear(); }
    public boolean isEmpty() { return records.isEmpty(); }

    public void add(ROD record) { add(record, false); }

    public void add(ROD record, boolean allowNameMismatch) {
        if ( record != null ) {
            if ( ! allowNameMismatch && ! name.equals(record.getName() ) )
                throw new StingException("Attempt to add ROD with non-matching name "+record.getName()+" to the track "+name);
        }
        records.add(record);
    }

    public void add(RODRecordList<ROD> records ) { add( records, false ); }

    public void add(RODRecordList<ROD> records, boolean allowNameMismatch) {
        for ( ROD record : records )
            add(record, allowNameMismatch);
    }    

    public int size() { return records.size() ; }

    /**
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object.
     * <p/>
     * <p>The implementor must ensure <tt>sgn(x.compareTo(y)) ==
     * -sgn(y.compareTo(x))</tt> for all <tt>x</tt> and <tt>y</tt>.  (This
     * implies that <tt>x.compareTo(y)</tt> must throw an exception iff
     * <tt>y.compareTo(x)</tt> throws an exception.)
     * <p/>
     * <p>The implementor must also ensure that the relation is transitive:
     * <tt>(x.compareTo(y)&gt;0 &amp;&amp; y.compareTo(z)&gt;0)</tt> implies
     * <tt>x.compareTo(z)&gt;0</tt>.
     * <p/>
     * <p>Finally, the implementor must ensure that <tt>x.compareTo(y)==0</tt>
     * implies that <tt>sgn(x.compareTo(z)) == sgn(y.compareTo(z))</tt>, for
     * all <tt>z</tt>.
     * <p/>
     * <p>It is strongly recommended, but <i>not</i> strictly required that
     * <tt>(x.compareTo(y)==0) == (x.equals(y))</tt>.  Generally speaking, any
     * class that implements the <tt>Comparable</tt> interface and violates
     * this condition should clearly indicate this fact.  The recommended
     * language is "Note: this class has a natural ordering that is
     * inconsistent with equals."
     * <p/>
     * <p>In the foregoing description, the notation
     * <tt>sgn(</tt><i>expression</i><tt>)</tt> designates the mathematical
     * <i>signum</i> function, which is defined to return one of <tt>-1</tt>,
     * <tt>0</tt>, or <tt>1</tt> according to whether the value of
     * <i>expression</i> is negative, zero or positive.
     *
     * @param that the object to be compared.
     * @return a negative integer, zero, or a positive integer as this object
     *         is less than, equal to, or greater than the specified object.
     * @throws ClassCastException if the specified object's type prevents it
     *                            from being compared to this object.
     */
    public int compareTo(RODRecordList<ROD> that) {
//        if ( this.getLocation() == null ) {
 //           if ( that.getLocation() == null )
  //      }
        return getLocation().compareTo(that.getLocation());  //To change body of implemented methods use File | Settings | File Templates.
    }
}

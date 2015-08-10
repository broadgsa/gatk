/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.refdata;

import org.broadinstitute.gatk.utils.refdata.utils.GATKFeature;
import org.broadinstitute.gatk.utils.refdata.utils.RODRecordList;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.HasGenomeLocation;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 10, 2009
 * Time: 6:10:48 PM
 * To change this template use File | Settings | File Templates.
 */
public class RODRecordListImpl extends AbstractList<GATKFeature> implements Comparable<RODRecordList>, Cloneable, RODRecordList, HasGenomeLocation {
    private List<GATKFeature> records;
    private GenomeLoc location = null;
    private String name = null;

    public RODRecordListImpl(String name) {
        records = new ArrayList<GATKFeature>();
        this.name = name;
    }

    /**
     * Fully qualified constructor: instantiates a new GATKFeatureRecordList object with specified GATKFeature track name, location on the
     * reference, and list of associated GATKFeatures. This is a knee-deep COPY constructor: passed name, loc, and data element
     * objects will be referenced from the created GATKFeatureRecordList (so that changing them from outside will affect data
     * in this object), however, the data elements will be copied into a newly
     * allocated list, so that the 'data' collection argument can be modified afterwards without affecting the state
     * of this record list. WARNING: this constructor is (semi-)validating: passed name and location
     * are allowed to be nulls (although it maybe unsafe, use caution), but if they are not nulls, then passed non-null GATKFeature data
     * elements must have same track name, and their locations must overlap with the passed 'location' argument. Null
     * data elements or null 'data' collection argument are allowed as well.
     * @param name the name of the track
     * @param data the collection of features at this location
     * @param loc the location
     */
    public RODRecordListImpl(String name, Collection<GATKFeature> data, GenomeLoc loc) {
        this.records = new ArrayList<GATKFeature>(data==null?0:data.size());
        this.name = name;
        this.location = loc;
        if ( data == null || data.size() == 0 ) return; // empty dataset, nothing to do
        for ( GATKFeature r : data ) {
            records.add(r);
            if ( r == null ) continue;
            if ( ! this.name.equals(r.getName() ) ) {
                throw new ReviewedGATKException("Attempt to add GATKFeature with non-matching name "+r.getName()+" to the track "+name);
            }
            if ( location != null && ! location.overlapsP(r.getLocation()) ) {
                    throw new ReviewedGATKException("Attempt to add GATKFeature that lies outside of specified interval "+location+"; offending GATKFeature:\n"+r.toString());
            }
        }
    }


    public GenomeLoc getLocation() { return location; }
    public String getName() { return name; }
    public Iterator<GATKFeature> iterator() { return records.iterator() ; }
    public void clear() { records.clear(); }
    public boolean isEmpty() { return records.isEmpty(); }

    public boolean add(GATKFeature record) { add(record, false); return true;}

    @Override
    public GATKFeature get(int i) {
        return records.get(i);
    }

    public void add(GATKFeature record, boolean allowNameMismatch) {
        if ( record != null ) {
            if ( ! allowNameMismatch && ! name.equals(record.getName() ) )
                throw new ReviewedGATKException("Attempt to add GATKFeature with non-matching name "+record.getName()+" to the track "+name);
        }
        records.add(record);
    }

    public void add(RODRecordList records ) { add( records, false ); }

    public void add(RODRecordList records, boolean allowNameMismatch) {
        for ( GATKFeature record : records )
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

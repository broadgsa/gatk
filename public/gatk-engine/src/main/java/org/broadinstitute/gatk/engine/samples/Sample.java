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

package org.broadinstitute.gatk.engine.samples;


import org.broadinstitute.gatk.utils.exceptions.UserException;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 *
 */
public class Sample implements Comparable<Sample> { // implements java.io.Serializable {
    final private String familyID, paternalID, maternalID;
    final private Gender gender;
    final private String otherPhenotype;
    final private Affection affection;
    final private String ID;
    final private SampleDB infoDB;
    final private Map<String, Object> properties = new HashMap<String, Object>();

    public final static String UNSET_QT = null;

    public Sample(final String ID, final SampleDB infoDB,
                  final String familyID, final String paternalID, final String maternalID,
                  final Gender gender, final Affection affection, final String otherPhenotype) {
        this.familyID = familyID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.gender = gender;
        this.otherPhenotype = otherPhenotype;
        this.affection = affection;
        this.ID = ID;
        this.infoDB = infoDB;
    }

    protected Sample(final String ID,
                     final String familyID, final String paternalID, final String maternalID,
                     final Gender gender, final Affection affection, final String otherPhenotype) {
        this(ID, null, familyID, paternalID, maternalID, gender, affection, otherPhenotype);
    }

    protected Sample(final String ID,
                     final String familyID, final String paternalID, final String maternalID,
                     final Gender gender, final Affection affection) {
        this(ID, null, familyID, paternalID, maternalID, gender, affection, UNSET_QT);
    }


    public Sample(final String ID, final SampleDB infoDB,
                  final String familyID, final String paternalID, final String maternalID, final Gender gender) {
        this(ID, infoDB, familyID, paternalID, maternalID, gender, Affection.UNKNOWN, UNSET_QT);
    }

    public Sample(final String ID, final SampleDB infoDB, final Affection affection, final String otherPhenotype) {
        this(ID, infoDB, null, null, null, Gender.UNKNOWN, affection, otherPhenotype);
    }

    public Sample(String id, SampleDB infoDB) {
        this(id, infoDB, null, null, null,
                Gender.UNKNOWN, Affection.UNKNOWN, UNSET_QT);
    }

    // -------------------------------------------------------------------------------------
    //
    // standard property getters
    //
    // -------------------------------------------------------------------------------------

    public String getID() {
        return ID;
    }

    public String getFamilyID() {
        return familyID;
    }

    public String getPaternalID() {
        return paternalID;
    }

    public String getMaternalID() {
        return maternalID;
    }

    public Affection getAffection() {
        return affection;
    }

    public boolean hasOtherPhenotype() {
        return affection == Affection.OTHER;
    }

    public String getOtherPhenotype() {
        return otherPhenotype;
    }

    /**
     * Get the sample's mother
     * @return sample object with relationship mother, if exists, or null
     */
    public Sample getMother() {
        return infoDB.getSample(maternalID);
    }

    /**
     * Get the sample's father
     * @return sample object with relationship father, if exists, or null
     */
    public Sample getFather() {
        return infoDB.getSample(paternalID);
    }

    public ArrayList<Sample> getParents(){
        ArrayList<Sample> parents = new ArrayList<Sample>(2);
        Sample parent = getMother();
        if(parent != null)
            parents.add(parent);
        parent = getFather();
        if(parent != null)
            parents.add(parent);
        return parents;
    }

    /**
     * Get gender of the sample
     * @return property of key "gender" - must be of type Gender
     */
    public Gender getGender() {
        return gender;
    }

    @Override
    public int compareTo(final Sample sample) {
        return ID.compareTo(sample.getID());
    }

    @Override
    public String toString() {
        return String.format("Sample %s fam=%s dad=%s mom=%s gender=%s affection=%s qt=%s props=%s",
                getID(), getFamilyID(), getPaternalID(), getMaternalID(), getGender(), getAffection(),
                getOtherPhenotype(), properties);
    }

//    // -------------------------------------------------------------------------------------
//    //
//    // code for working with additional -- none standard -- properites
//    //
//    // -------------------------------------------------------------------------------------
//
//    public Map<String, Object> getExtraProperties() {
//        return Collections.unmodifiableMap(properties);
//    }
//
//    /**
//     * Get one property
//     * @param key key of property
//     * @return value of property as generic object
//     */
//    public Object getExtraPropertyValue(final String key) {
//        return properties.get(key);
//    }
//
//    /**
//     *
//     * @param key property key
//     * @return true if sample has this property (even if its value is null)
//     */
//    public boolean hasExtraProperty(String key) {
//        return properties.containsKey(key);
//    }

    @Override
    public int hashCode() {
        return ID.hashCode();
    }

    @Override
    public boolean equals(final Object o) {
        if(o == null)
            return false;
        if(o instanceof Sample) {
            Sample otherSample = (Sample)o;
            return ID.equals(otherSample.ID) &&
                    equalOrNull(familyID, otherSample.familyID) &&
                    equalOrNull(paternalID, otherSample.paternalID) &&
                    equalOrNull(maternalID, otherSample.maternalID) &&
                    equalOrNull(gender, otherSample.gender) &&
                    equalOrNull(otherPhenotype, otherSample.otherPhenotype) &&
                    equalOrNull(affection, otherSample.affection) &&
                    equalOrNull(properties, otherSample.properties);
        }
        return false;
    }

    private final static boolean equalOrNull(final Object o1, final Object o2) {
        if ( o1 == null )
            return o2 == null;
        else
            return o2 == null ? false : o1.equals(o2);
    }

    private final static <T> T mergeValues(final String name, final String field, final T o1, final T o2, final T emptyValue) {
        if ( o1 == null || o1.equals(emptyValue) ) {
            // take o2 if both are null, otherwise keep o2
            return o2 == null ? null : o2;
        } else {
            if ( o2 == null || o2.equals(emptyValue) )
                return o1; // keep o1, since it's a real value
            else {
                // both o1 and o2 have a value
                if ( o1 instanceof String && o1.equals(o2) )
                    return o1;
                else if ( o1 == o2 )
                    return o1;
                else
                    throw new UserException("Inconsistent values detected for " + name + " for field " + field + " value1 " + o1 + " value2 " + o2);
            }
        }
    }

    public final static Sample mergeSamples(final Sample prev, final Sample next) {
        if ( prev.equals(next) )
            return next;
        else {
            return new Sample(prev.getID(), prev.infoDB,
                    mergeValues(prev.getID(), "Family_ID", prev.getFamilyID(), next.getFamilyID(), null),
                    mergeValues(prev.getID(), "Paternal_ID", prev.getPaternalID(), next.getPaternalID(), null),
                    mergeValues(prev.getID(), "Material_ID", prev.getMaternalID(), next.getMaternalID(), null),
                    mergeValues(prev.getID(), "Gender", prev.getGender(), next.getGender(), Gender.UNKNOWN),
                    mergeValues(prev.getID(), "Affection", prev.getAffection(), next.getAffection(), Affection.UNKNOWN),
                    mergeValues(prev.getID(), "OtherPhenotype", prev.getOtherPhenotype(), next.getOtherPhenotype(), UNSET_QT));
                    //mergeValues(prev.getID(), "ExtraProperties", prev.getExtraProperties(), next.getExtraProperties(), Collections.emptyMap()));
        }
    }
}

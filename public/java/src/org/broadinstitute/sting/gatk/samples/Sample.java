package org.broadinstitute.sting.gatk.samples;


import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 *
 */
public class Sample implements java.io.Serializable {
    final private String familyID, paternalID, maternalID;
    final private Gender gender;
    final private double quantitativePhenotype;
    final private Affection affection;
    final private String ID;
    final private SampleDataSource dataSource;
    final private Map<String, Object> properties = new HashMap<String, Object>();

    public final static double UNSET_QT = Double.NaN;

    public Sample(final String ID, final SampleDataSource dataSource,
                  final String familyID, final String paternalID, final String maternalID,
                  final Gender gender, final Affection affection, final double quantitativePhenotype) {
        this.familyID = familyID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.gender = gender;
        this.quantitativePhenotype = quantitativePhenotype;
        this.affection = affection;
        this.ID = ID;
        this.dataSource = dataSource;
    }

    protected Sample(final String ID,
                     final String familyID, final String paternalID, final String maternalID,
                     final Gender gender, final Affection affection, final double quantitativePhenotype) {
        this(ID, null, familyID, paternalID, maternalID, gender, affection, quantitativePhenotype);
    }

    protected Sample(final String ID,
                     final String familyID, final String paternalID, final String maternalID,
                     final Gender gender, final Affection affection) {
        this(ID, null, familyID, paternalID, maternalID, gender, affection, UNSET_QT);
    }


    public Sample(final String ID, final SampleDataSource dataSource,
                  final String familyID, final String paternalID, final String maternalID, final Gender gender) {
        this(ID, dataSource, familyID, paternalID, maternalID, gender, Affection.UNKNOWN, UNSET_QT);
    }

    public Sample(final String ID, final SampleDataSource dataSource, final Affection affection, final double quantitativePhenotype) {
        this(ID, dataSource, null, null, null, Gender.UNKNOWN, affection, quantitativePhenotype);
    }

    public Sample(String id, SampleDataSource dataSource) {
        this(id, dataSource, null, null, null,
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

    public boolean hasQuantitativeTrait() {
        return affection == Affection.QUANTITATIVE;
    }

    public double getQuantitativePhenotype() {
        return quantitativePhenotype;
    }

    /**
     * Get the sample's mother
     * @return sample object with relationship mother, if exists, or null
     */
    public Sample getMother() {
        return dataSource.getSample(maternalID);
    }

    /**
     * Get the sample's father
     * @return sample object with relationship father, if exists, or null
     */
    public Sample getFather() {
        return dataSource.getSample(paternalID);
    }

    /**
     * Get gender of the sample
     * @return property of key "gender" - must be of type Gender
     */
    public Gender getGender() {
        return gender;
    }

    @Override
    public String toString() {
        return String.format("Sample %s fam=%s dad=%s mom=%s gender=%s affection=%s qt=%s props=%s",
                getID(), getFamilyID(), getPaternalID(), getMaternalID(), getGender(), getAffection(),
                getQuantitativePhenotype(), getExtraProperties());
    }

    // -------------------------------------------------------------------------------------
    //
    // code for working with additional -- none standard -- properites
    //
    // -------------------------------------------------------------------------------------

    public Map<String, Object> getExtraProperties() {
        return Collections.unmodifiableMap(properties);
    }

    /**
     * Get one property
     * @param key key of property
     * @return value of property as generic object
     */
    public Object getExtraPropertyValue(final String key) {
        return properties.get(key);
    }

    /**
     *
     * @param key property key
     * @return true if sample has this property (even if its value is null)
     */
    public boolean hasExtraProperty(String key) {
        return properties.containsKey(key);
    }

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
                    equalOrNull(quantitativePhenotype, otherSample.quantitativePhenotype) &&
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
}

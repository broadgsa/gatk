package org.broadinstitute.sting.gatk.samples;


import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Jul 26, 2010
 * Time: 3:31:38 PM
 */
public class Sample implements java.io.Serializable {
    final private String familyID, paternalID, maternalID;
    final private Sample.Gender gender;
    final private double quantitativePhenotype;
    final private Sample.Affection affection;
    final private String ID;
    final private SampleDataSource dataSource;

    // todo -- conditionally add the property map -- should be empty by default
    private final Map<String, Object> properties = new HashMap<String, Object>();

    public enum Gender {
        MALE,
        FEMALE,
        UNKNOWN
    }

    public enum Affection {
        /** Status is unknown */
        UNKNOWN,
        /** Suffers from the disease */
        AFFECTED,
        /** Unaffected by the disease */
        UNAFFECTED,
        /** A quantitative trait: value of the trait is stored elsewhere */
        QUANTITATIVE
    }

    public final static double UNSET_QUANTITIATIVE_TRAIT_VALUE = Double.NaN;

    public Sample(final String ID, final SampleDataSource dataSource,
                  final String familyID, final String paternalID, final String maternalID,
                  final Gender gender, final double quantitativePhenotype, final Affection affection) {
        this.familyID = familyID;
        this.paternalID = paternalID;
        this.maternalID = maternalID;
        this.gender = gender;
        this.quantitativePhenotype = quantitativePhenotype;
        this.affection = affection;
        this.ID = ID;
        this.dataSource = dataSource;
    }

    public Sample(final String ID, final SampleDataSource dataSource,
                  final String familyID, final String paternalID, final String maternalID, final Gender gender) {
        this(ID, dataSource, familyID, paternalID, maternalID, gender,
                UNSET_QUANTITIATIVE_TRAIT_VALUE, Affection.UNKNOWN);
    }

    public Sample(final String ID, final SampleDataSource dataSource, final double quantitativePhenotype, final Affection affection) {
        this(ID, dataSource, null, null, null, Gender.UNKNOWN, quantitativePhenotype, affection);
    }

    public Sample(String id, SampleDataSource dataSource) {
        this(id, dataSource,
                null, null, null,
                Gender.UNKNOWN, UNSET_QUANTITIATIVE_TRAIT_VALUE, Affection.UNKNOWN);
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

    public String getFamilyId() {
        return familyID;
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

//    @Override
//    public boolean equals(Object o) {
//        if (this == o) return true;
//        if (o == null || getClass() != o.getClass()) return false;
//
//        Sample sample = (Sample) o;
//        if (ID != null ? !ID.equals(sample.ID) : sample.ID != null) return false;
//        if (properties != null ? !properties.equals(sample.properties) : sample.properties != null) return false;
//
//        return true;
//    }
//
//    @Override
//    public int hashCode() {
//        return ID != null ? ID.hashCode() : "".hashCode();
//    }
}

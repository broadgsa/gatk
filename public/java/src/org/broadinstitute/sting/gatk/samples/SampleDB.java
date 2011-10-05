package org.broadinstitute.sting.gatk.samples;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.util.*;

/**
 *
 */
public class SampleDB {
    /**
     * This is where Sample objects are stored. Samples are usually accessed by their ID, which is unique, so
     * this is stored as a HashMap.
     */
    private final HashMap<String, Sample> samples = new HashMap<String, Sample>();

    /**
     * Constructor takes both a SAM header and sample files because the two must be integrated.
     */
    public SampleDB() {

    }

    /**
     * Protected function to add a single sample to the database
     *
     * @param sample to be added
     */
    protected SampleDB addSample(Sample sample) {
        Sample prev = samples.get(sample.getID());
        if ( prev != null )
            sample = Sample.mergeSamples(prev, sample);
        samples.put(sample.getID(), sample);
        return this;
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for getting a sample from the DB
    //
    // --------------------------------------------------------------------------------

    /**
     * Get a sample by its ID
     * If an alias is passed in, return the main sample object 
     * @param id
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(String id) {
        return samples.get(id);
    }

    /**
     *
     * @param read
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(final SAMRecord read) {
        return getSample(read.getReadGroup());
    }

    /**
     *
     * @param rg
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(final SAMReadGroupRecord rg) {
        return getSample(rg.getSample());
    }

    /**
     * @param g Genotype
     * @return sample Object with this ID, or null if this does not exist
     */
    public Sample getSample(final Genotype g) {
        return getSample(g.getSampleName());
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for accessing samples in the DB
    //
    // --------------------------------------------------------------------------------

    /**
     * Get number of sample objects
     * @return size of samples map
     */
    public int sampleCount() {
        return samples.size();
    }

    public Set<Sample> getSamples() {
        return new HashSet<Sample>(samples.values());
    }

    public Collection<String> getSampleNames() {
        return Collections.unmodifiableCollection(samples.keySet());
    }


    /**
     * Takes a collection of sample names and returns their corresponding sample objects
     * Note that, since a set is returned, if you pass in a list with duplicates names there will not be any duplicates in the returned set
     * @param sampleNameList Set of sample names
     * @return Corresponding set of samples
     */
    public Set<Sample> getSamples(Collection<String> sampleNameList) {
        HashSet<Sample> samples = new HashSet<Sample>();
        for (String name : sampleNameList) {
            try {
                samples.add(getSample(name));
            }
            catch (Exception e) {
                throw new StingException("Could not get sample with the following ID: " + name, e);
            }
        }
        return samples;
    }

    // --------------------------------------------------------------------------------
    //
    // Higher level pedigree functions
    //
    // --------------------------------------------------------------------------------

    /**
     * Returns a sorted set of the family IDs in all samples (excluding null ids)
     * @return
     */
    public final Set<String> getFamilyIDs() {
        return getFamilies().keySet();
    }

    /**
     * Returns a map from family ID -> set of family members for all samples with
     * non-null family ids
     *
     * @return
     */
    public final Map<String, Set<Sample>> getFamilies() {
        final Map<String, Set<Sample>> families = new TreeMap<String, Set<Sample>>();

        for ( final Sample sample : samples.values() ) {
            final String famID = sample.getFamilyID();
            if ( famID != null ) {
                if ( ! families.containsKey(famID) )
                    families.put(famID, new TreeSet<Sample>());
                families.get(famID).add(sample);
            }
        }

        return families;
    }

    /**
     * Return all samples with a given family ID
     * @param familyId
     * @return
     */
    public Set<Sample> getFamily(String familyId) {
        return getFamilies().get(familyId);
    }

    /**
     * Returns all children of a given sample
     * See note on the efficiency of getFamily() - since this depends on getFamily() it's also not efficient
     * @param sample
     * @return
     */
    public Set<Sample> getChildren(Sample sample) {
        final HashSet<Sample> children = new HashSet<Sample>();
        for ( final Sample familyMember : getFamily(sample.getFamilyID())) {
            if ( familyMember.getMother() == sample || familyMember.getFather() == sample ) {
                children.add(familyMember);
            }
        }
        return children;
    }
}

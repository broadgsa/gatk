package org.broadinstitute.sting.gatk.samples;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.utils.SampleUtils;
import org.broadinstitute.sting.utils.exceptions.StingException;
import org.broadinstitute.sting.utils.variantcontext.Genotype;

import java.io.File;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: brett
 * Date: Jul 26, 2010
 * Time: 3:30:09 PM
 *
 * This class stores and manages sample metadata. This data is encoded in a sample file, which can be included
 * in the GATK by the "--samples" argument. This class reads and parses those files.
 *
 * Although there are a set of public methods for accessing sample data, they aren't used by walkers - they are really
 * only used by GenomeAnalysisEngine. An instance of GenomeAnalysisEngine has one SampleDataSource. When a walker
 * wants to access sample data, it asks GenomeAnalysis to fetch this data from its SampleDataSource.
 *
 */
public class SampleDataSource {
    /**
     * This is where Sample objects are stored. Samples are usually accessed by their ID, which is unique, so
     * this is stored as a HashMap.
     */
    private final HashMap<String, Sample> samples = new HashMap<String, Sample>();

    /**
     * Constructor takes both a SAM header and sample files because the two must be integrated.
     */
    public SampleDataSource() {
        samples.put(null, new Sample(null, this));
    }

    public SampleDataSource(final SAMFileHeader header, final List<File> sampleFiles) {
        this();
        addSamples(header);
        addSamples(sampleFiles);
    }

    // --------------------------------------------------------------------------------
    //
    // Functions for adding samples to the DB
    //
    // TODO: these should be protected, really
    //
    // --------------------------------------------------------------------------------

    /**
     * Hallucinates sample objects for all the samples in the SAM file and stores them
     */
    public SampleDataSource addSamples(SAMFileHeader header) {
        for (String sampleName : SampleUtils.getSAMFileSamples(header)) {
            if (getSample(sampleName) == null) {
                Sample newSample = new Sample(sampleName, this);
                samples.put(sampleName, newSample);
            }
        }
        return this;
    }

    public SampleDataSource addSamples(final List<File> sampleFiles) {
        // add files consecutively
        for (File file : sampleFiles) {
            addSamples(file);
        }
        return this;
    }

    /**
     * Parse one sample file and integrate it with samples that are already there
     * Fail quickly if we find any errors in the file
     */
    public SampleDataSource addSamples(File sampleFile) {
        return this;
    }

    /**
     * Add a sample to the collection
     * @param sample to be added
     */
    private SampleDataSource addSample(Sample sample) {
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

    /**
     * Return all samples with a given family ID
     * Note that this isn't terribly efficient (linear) - it may be worth adding a new family ID data structure for this
     * @param familyId
     * @return
     */
    public Set<Sample> getFamily(String familyId) {
        HashSet<Sample> familyMembers = new HashSet<Sample>();

        for (Sample sample : samples.values()) {
            if (sample.getFamilyId() != null) {
                if (sample.getFamilyId().equals(familyId))
                    familyMembers.add(sample);
            }
        }
        return familyMembers;
    }

    /**
     * Returns all children of a given sample
     * See note on the efficiency of getFamily() - since this depends on getFamily() it's also not efficient
     * @param sample
     * @return
     */
    public Set<Sample> getChildren(Sample sample) {
        HashSet<Sample> children = new HashSet<Sample>();
        for (Sample familyMember : getFamily(sample.getFamilyId())) {
            if (familyMember.getMother() == sample || familyMember.getFather() == sample) {
                children.add(familyMember);
            }
        }
        return children;
    }

    public Collection<Sample> getSamples() {
        return Collections.unmodifiableCollection(samples.values());
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
}

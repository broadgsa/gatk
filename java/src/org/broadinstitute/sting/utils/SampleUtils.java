package org.broadinstitute.sting.utils;

import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMFileHeader;

import java.util.*;

import org.broadinstitute.sting.utils.genotype.Genotype;
import org.broadinstitute.sting.utils.genotype.SampleBacked;


/**
 * SampleUtils is a static class (no instantiation allowed!) with some utility methods for getting samples
 * quality scores.
 *
 * @author ebanks
 */
public class SampleUtils {
    /**
     * Private constructor.  No instantiating this class!
     */
    private SampleUtils() {}

    /**
     * Pull out the samples from a SAMFileHeader;
     * note that we use a TreeSet so that they are sorted
     *
     * @param header  the sam file header
     * @return list of strings representing the sample names
     */
    public static Set<String> getSAMFileSamples(SAMFileHeader header) {
        // get all of the unique sample names
        Set<String> samples = new TreeSet<String>();
        List<SAMReadGroupRecord> readGroups = header.getReadGroups();
        for ( SAMReadGroupRecord readGroup : readGroups )
            samples.add(readGroup.getSample());
        return samples;
    }

    /**
     * get the samples names from genotype objects if they are backed by samples
     *
     * @param genotypes the genotype list
     * @return list of strings representing the sample names
     */
    public static List<String> getGenotypeSamples(List<Genotype> genotypes) {
        List<String> samples = new ArrayList<String>();
        for ( Genotype genotype : genotypes ) {
            if ( genotype instanceof SampleBacked )
                samples.add(((SampleBacked)genotype).getSampleName());
        }
        return samples;
    }
}
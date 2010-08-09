package org.broadinstitute.sting.gatk.refdata.tracks.builders;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.broad.tribble.Feature;
import org.broad.tribble.index.Index;
import org.broad.tribble.index.linear.LinearIndex;
import org.broad.tribble.iterators.CloseableTribbleIterator;
import org.broad.tribble.source.BasicFeatureSource;
import org.broad.tribble.util.LittleEndianOutputStream;
import org.broad.tribble.vcf.VCFCodec;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableCodec;
import org.broadinstitute.sting.gatk.refdata.features.annotator.AnnotatorInputTableFeature;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.collections.Pair;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.io.*;
import java.util.*;

/**
 * performance tests for different index types
 */
public class IndexPerformanceTests extends BaseTest {
    // the RMD track builder
    private TribbleRMDTrackBuilder builder;

    // set the logger level
    Logger logger = Logger.getLogger(IndexPerformanceTests.class);

    // the input files to test
    Map<String, File> inputFiles = new LinkedHashMap<String,File>();

    // the input types
    Map<String, Class> inputTypes = new HashMap<String,Class>();

    // where the vcf files are located
    String fileLocation = validationDataLocation + "Index_Performance_Data/";

    // bin sizes to try
    int[] binSizes = {10, 100, 1000, 5000, 10000, 50000};

    PrintWriter writer;
    PrintWriter writer2;
    /** setup the files we're going to run with, including their names */
    @Before
    public void setupFilesAndIndexes() {
        logger.setLevel(Level.INFO);
        builder = new TribbleRMDTrackBuilder();
        IndexedFastaSequenceFile seq = new IndexedFastaSequenceFile(new File(hg18Reference));
        GenomeLocParser.setupRefContigOrdering(seq);

        // the input files
        inputFiles.put("\"10\"",new File(fileLocation + "tip10.vcf"));
        inputFiles.put("\"100\"",new File(fileLocation + "tip100.vcf"));
        inputFiles.put("\"1,000\"",new File(fileLocation + "tip1000.vcf"));
        inputFiles.put("\"10,000\"",new File(fileLocation + "tip10000.vcf"));
        inputFiles.put("\"100,000\"",new File(fileLocation + "tip100000.vcf"));
        inputFiles.put("\"1,000,000\"",new File(fileLocation + "tip1000000.vcf"));

        for (String name : inputFiles.keySet()) {
            inputTypes.put(name,VCFCodec.class);
        }
        inputFiles.put("Big Table",new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/slowAnnotator/big.table.txt"));
        inputTypes.put("Big Table", AnnotatorInputTableCodec.class);
    }

    @Test
    public void emptyTest() {
        // do nothing
    }

    //@Test
    public void performanceTest() {
        try {
            writer = new PrintWriter(new FileWriter("testOutput_linear.txt"));
            writer2 = new PrintWriter(new FileWriter("testOutput_tree.txt"));
        } catch (IOException e) {
            Assert.fail("Unable to open file testOutput.txt");
        }
        writer.println("name,index,binSize,createTime,seekTime,thousandPerThousand,record_count,index_size");
        writer2.println("name,index,binSize,createTime,seekTime,thousandPerThousand,record_count,index_size");
        for (String name : inputFiles.keySet()) {
            for (int size : binSizes) {
                System.err.println("running " + name + " with bin size " + size);
                printTestLine(name, true, size);
                printTestLine(name, false, size);
            }
        }
        writer.close();
        writer2.close();
    }

    private void printTestLine(String name, boolean useLinear, int size) {
        PrintWriter wr = (useLinear) ? writer : writer2;
        List<Long> values = performIndexTest(name,useLinear, size);
        wr.print(name + "," + ((useLinear) ? "linear" : "tree") + "," + size);
        for (Long l : values) {
            wr.print(",");
            wr.print(l);
        }
        wr.println();
    }

    /**
     * time various tasks using the specified index
     * @param   name the name to get
     * @return  a five-piece: the time to create the index, the time to seek to chromosome 1, and the time to process reading
     *          every other 1000 bases of chr1 (of the first 100M), the count of records seen in the last operation, and the index size
     */
    public List<Long> performIndexTest(String name, boolean useLinear, int size) {
        TribbleRMDTrackBuilder.useLinearIndex = useLinear;
        TribbleRMDTrackBuilder.binSize = size;

        deleteIndex(inputFiles.get(name));
        // time creating the index
        long createTime = System.currentTimeMillis();
        Pair<BasicFeatureSource, SAMSequenceDictionary> pairing = builder.createFeatureReader(inputTypes.get(name),inputFiles.get(name));
        createTime = System.currentTimeMillis() - createTime;
        //System.err.println("index creation took " + createTime);

        // seek to chr1
        long seekTo1 = seekToChr1(pairing);

        // seek every 1000 bases in Chr1
        long count = 0;
        long thousandEveryThousand = System.currentTimeMillis();
        try {
            for (int x = 1; x < 1000000; x = x + 1000) {
                //CloseableTribbleIterator<Feature> iter = pairing.first.query("chr1", x+(int)Math.floor(Math.random()*1000), x+1000); // query
                CloseableTribbleIterator<Feature> iter = pairing.first.query("chr1", x, x+1000); // query
                for (Feature feat : iter) {
                    count++;
                }
            }

        } catch (IOException e) {
            Assert.fail("Unable to load file for query!!");
        }
        thousandEveryThousand = System.currentTimeMillis() - thousandEveryThousand;
        //System.err.println("thousand every thousand (for first million) took " + thousandEveryThousand);
        return Arrays.asList(createTime,seekTo1,thousandEveryThousand,count,new File(inputFiles.get(name) + ".idx").length());
    }

    private long seekToChr1(Pair<BasicFeatureSource, SAMSequenceDictionary> pairing) {
        // time seeking to the first 1M bases of Chr1
        long seekTo1 = System.currentTimeMillis();
        try {
            CloseableTribbleIterator iter = pairing.first.query("chr1",1,10000000); // query
        } catch (IOException e) {
            Assert.fail("Unable to load file for query!!");
        }
        seekTo1 = System.currentTimeMillis() - seekTo1;
        //System.err.println("seeking to chr1 took " + seekTo1);
        return seekTo1;
    }

    //@Test
    public void testBigTable() {
        // store the mapping of location -> variant; this will get big
        Map<GenomeLoc, Integer> features = new TreeMap<GenomeLoc,Integer>();

        Map<Integer,Integer> bucketToCount = getMapOfFeatures(features,false);
        Pair<BasicFeatureSource, SAMSequenceDictionary> pairing;

        Map<GenomeLoc, Integer> features2 = new TreeMap<GenomeLoc,Integer>();
        Map<Integer,Integer> bucketToCount2 = getMapOfFeatures(features2,true);

        System.err.println("Summary: ");
        System.err.println("Summary: tree " + features.size());
        System.err.println("Summary: linear " + features2.size());

        // compare the two
        for (Map.Entry<GenomeLoc,Integer> entry: features.entrySet()) {
            if (!features2.containsKey(entry.getKey())) {
                System.err.println("key " + entry + " missing from linear, count " + entry.getValue());

            }
            else if (features2.get(entry.getKey()) != entry.getValue()) {
                /*System.err.println("counts are not equal at " +
                        entry.getKey() +
                        " features2.get(entry.getKey()) = " +
                        features2.get(entry.getKey()) +
                        " feature1 = " + entry.getValue());*/
            }
            if (features2.containsKey(entry.getKey())) features2.remove(entry.getKey());
        }
        System.err.println("Missing from the tree :");
        for (Map.Entry<GenomeLoc,Integer> entry2: features2.entrySet()) {
            System.err.println("Position " + entry2.getKey() + " count = " + entry2.getValue());
        }

        for (Integer bucket : bucketToCount.keySet()) {
            if (!bucketToCount2.get(bucket).equals(bucketToCount.get(bucket))) {
                System.err.println("Bucket " + bucket + " tree != linear, " + bucketToCount2.get(bucket) + " " + bucketToCount.get(bucket));
            }
        }
    }

    private Map<Integer,Integer> getMapOfFeatures(Map<GenomeLoc, Integer> features, boolean useLinear) {
        File bigTable = new File("/humgen/gsa-hpprojects/dev/depristo/oneOffProjects/slowAnnotator/big.table.txt");
        TribbleRMDTrackBuilder.useLinearIndex = useLinear;
        TribbleRMDTrackBuilder.binSize = 1000;

        deleteIndex(inputFiles.get("Big Table"));
        // time creating the index
        logger.warn("creating index");

        Map<Integer,Integer> bucketToCount = new TreeMap<Integer,Integer>();
        Pair<BasicFeatureSource, SAMSequenceDictionary> pairing = builder.createFeatureReader(inputTypes.get("Big Table"),inputFiles.get("Big Table"));
        try {
            for (int x = 5000; x < 6000; x = x + 1000) {
                int bucketCount = 0;
                CloseableTribbleIterator<Feature> iter = pairing.first.query("chr1", x, x+1000); // query
                for (Feature feat : iter) {
                    GenomeLoc loc = GenomeLocParser.createGenomeLoc(feat.getChr(),feat.getStart(),feat.getEnd());
                    if (loc.getStop() < 5000 || loc.getStart() > 6000) continue;
                    int count = 0;
                    if (features.containsKey(loc))
                        count = features.get(loc)+1;
                    features.put(loc,count);
                    bucketCount++;
                }
                bucketToCount.put(x,bucketCount);
            }
        } catch (IOException e) {
            Assert.fail("Unable to load file for query!!");
        }
        return bucketToCount;
    }

    //@Test
    public void testGetTreeIndexLocation() {
        File bigTable = new File("small.table.txt");
        TribbleRMDTrackBuilder.useLinearIndex = false;
        TribbleRMDTrackBuilder.binSize = 1000;

        deleteIndex(bigTable);
        // time creating the index
        logger.warn("creating index");

        Map<Integer,Integer> bucketToCount = new TreeMap<Integer,Integer>();
        Pair<BasicFeatureSource, SAMSequenceDictionary> pairing = builder.createFeatureReader(inputTypes.get("Big Table"),bigTable);
        try {
            int count= 0;
            CloseableTribbleIterator<Feature> iter = null;
            for (int x = 5000; x < 6000; x = x + 1000)
                iter = pairing.first.query("chr1", x, x+1000); // query
                for (Feature feat : iter) {
                    GenomeLoc loc = GenomeLocParser.createGenomeLoc(feat.getChr(),feat.getStart(),feat.getEnd());
                    if (loc.getStop() < 5000 || loc.getStart() > 6000) continue;
                    count++;
                }
        System.err.println(count);
        } catch (IOException e) {
            Assert.fail("Unable to load file for query!!");
        }
    }


    private void deleteIndex(File fl) {
        File indexFile = new File(fl + TribbleRMDTrackBuilder.indexExtension);
        boolean deleted = true;
        if (indexFile.exists())
            deleted = indexFile.delete();
        if (!deleted)
            Assert.fail("Unable to delete index file");
    }

}

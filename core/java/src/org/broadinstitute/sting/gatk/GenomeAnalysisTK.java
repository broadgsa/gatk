package org.broadinstitute.sting.gatk;

import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.RuntimeIOException;
import edu.mit.broad.picard.cmdline.CommandLineProgram;
import edu.mit.broad.picard.cmdline.Usage;
import edu.mit.broad.picard.cmdline.Option;
import edu.mit.broad.picard.reference.ReferenceSequenceFileFactory;
import edu.mit.broad.picard.reference.ReferenceSequenceFile;
import edu.mit.broad.picard.reference.ReferenceSequence;

import org.broadinstitute.sting.utils.*;
import org.broadinstitute.sting.gatk.walkers.*;
import org.broadinstitute.sting.gatk.refdata.ReferenceOrderedData;
import org.broadinstitute.sting.gatk.refdata.rodDbSNP;
import org.broadinstitute.sting.gatk.refdata.rodGFF;
import org.broadinstitute.sting.gatk.iterators.ReferenceIterator;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.FastaSequenceFile2;

import java.io.*;
import java.util.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Random;

public class GenomeAnalysisTK extends CommandLineProgram {
    // Usage and parameters
    @Usage(programVersion="0.1") public String USAGE = "Genome Analysis Toolkit\n";
    @Option(shortName="I", doc="SAM or BAM file for validation") public File INPUT_FILE;
    @Option(shortName="M", doc="Maximum number of reads to process before exiting", optional=true) public String MAX_READS_ARG = "-1";
    @Option(shortName="S", doc="How strict should we be with validation", optional=true) public String STRICTNESS_ARG = "strict";
    @Option(shortName="R", doc="Reference sequence file", optional=true) public File REF_FILE_ARG = null;
    @Option(shortName="B", doc="Debugging output", optional=true) public String DEBUGGING_STR = null;
    @Option(shortName="L", doc="Genome region to operation on: from chr:start-end", optional=true) public String REGION_STR = null;

    @Option(shortName="INT", doc="File containing list of genomic intervals to operate on. line := <contig> <start> <end>\n", optional=true) public String INTERVALS_FILE = null;

    @Option(shortName="T", doc="Type of analysis to run") public String Analysis_Name = null;
    @Option(shortName="DBSNP", doc="DBSNP file", optional=true) public String DBSNP_FILE = null;
    @Option(shortName="THREADED_IO", doc="If true, enables threaded I/O operations", optional=true) public String ENABLED_THREADED_IO = "false";
    @Option(shortName="U", doc="If true, enables unsafe operations, nothing will be checked at runtime.  You better know what you are doing if you set this flag.", optional=false) public String UNSAFE = "false";
    @Option(shortName="SORT_ON_FLY", doc="If true, enables on fly sorting of reads file.", optional=false) public String ENABLED_SORT_ON_FLY = "false";

    @Option(shortName="PLUGINS", doc="Directory where plugin class files live.", optional=true)
    public String pluginPathName = null;

    private WalkerManager walkerManager = null;
    
    private TraversalEngine engine = null;
    public boolean DEBUGGING = false;

    /** Required main method implementation. */
    public static void main(String[] argv) {
        System.exit(new GenomeAnalysisTK().instanceMain(argv));
    }

    protected int doWork() {
        walkerManager = new WalkerManager(pluginPathName);

        //testNewReferenceFeatures(REF_FILE_ARG);

        final boolean TEST_ROD = false;
        ReferenceOrderedData[] rods = null;

        if ( TEST_ROD ) {
            ReferenceOrderedData gff = new ReferenceOrderedData(new File("trunk/data/gFFTest.gff"), rodGFF.class );
            gff.testMe();

            //ReferenceOrderedData dbsnp = new ReferenceOrderedData(new File("trunk/data/dbSNP_head.txt"), rodDbSNP.class );
            ReferenceOrderedData dbsnp = new ReferenceOrderedData(new File("/Volumes/Users/mdepristo/broad/ATK/exampleSAMs/dbSNP_chr20.txt"), rodDbSNP.class );
            //dbsnp.testMe();
            rods = new ReferenceOrderedData[] { dbsnp }; // { gff, dbsnp };
        }
        else if ( DBSNP_FILE != null ) {
            ReferenceOrderedData dbsnp = new ReferenceOrderedData(new File(DBSNP_FILE), rodDbSNP.class );
            //dbsnp.testMe();
            rods = new ReferenceOrderedData[] { dbsnp }; // { gff, dbsnp };
        }
        else {
            rods = new ReferenceOrderedData[] {}; // { gff, dbsnp };
        }

        this.engine = new TraversalEngine(INPUT_FILE, REF_FILE_ARG, rods);

        // Prepare the sort ordering w.r.t. the sequence dictionary
        final ReferenceSequenceFile refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(REF_FILE_ARG);
        List<SAMSequenceRecord> refContigs = refFile.getSequenceDictionary().getSequences();
        HashMap<String, Integer> refContigOrdering = new HashMap<String, Integer>();
        int i = 0;
        for ( SAMSequenceRecord contig : refContigs ) {
            refContigOrdering.put(contig.getSequenceName(), i);
            i++;
        }
        GenomeLoc.setContigOrdering(refContigOrdering);

        ValidationStringency strictness;
    	if ( STRICTNESS_ARG == null ) {
            strictness = ValidationStringency.STRICT;
    	}
    	else if ( STRICTNESS_ARG.toLowerCase().equals("lenient") ) {
    		strictness = ValidationStringency.LENIENT;
    	}
    	else if ( STRICTNESS_ARG.toLowerCase().equals("silent") ) {
    		strictness = ValidationStringency.SILENT;
    	}
    	else {
            strictness = ValidationStringency.STRICT;
    	}
        System.err.println("Strictness is " + strictness);
        engine.setStrictness(strictness);

        engine.setDebugging(! ( DEBUGGING_STR == null || DEBUGGING_STR.toLowerCase().equals("true")));
        engine.setMaxReads(Integer.parseInt(MAX_READS_ARG));

        if ( REGION_STR != null ) {
            engine.setLocation(REGION_STR);
        }

        if (INTERVALS_FILE != null)
        {
            engine.setLocationFromFile(INTERVALS_FILE);
        }

        engine.setSafetyChecking(! UNSAFE.toLowerCase().equals("true"));
        engine.setSortOnFly(ENABLED_SORT_ON_FLY.toLowerCase().equals("true"));

        engine.initialize(ENABLED_THREADED_IO.toLowerCase().equals("true"));

        //engine.testReference();

        //LocusWalker<Integer,Integer> walker = new PileupWalker();

        // Try to get the walker specified
        Object my_walker;
        if (walkerManager.doesWalkerExist(Analysis_Name)) {
            my_walker = walkerManager.getWalkerByName(Analysis_Name);
        } else {
            System.out.println("Could not find walker "+Analysis_Name);
            return 0;
        }

        try {
            LocusWalker<?, ?> walker = (LocusWalker<?, ?>)my_walker;
            engine.traverseByLoci(walker);
        }
        catch ( java.lang.ClassCastException e ) {
            // I guess we're a read walker LOL
            ReadWalker<?, ?> walker = (ReadWalker<?, ?>)my_walker;
            engine.traverseByRead(walker);
        }

        return 0;
    }

    /**
     * An inappropriately placed validation and performance testing routine for jumping
     * around in the fasta sequence file.
     * @param refFileName
     */
    private static void testNewReferenceFeatures(final File refFileName) {
        final FastaSequenceFile2 refFile = new FastaSequenceFile2(refFileName);
        Utils.setupRefContigOrdering(refFile);

        List<SAMSequenceRecord> refContigs = refFile.getSequenceDictionary().getSequences();

        /*
        for ( SAMSequenceRecord refContig: refContigs ) {
            System.out.printf("  Traversing from chr1 to %s would require jumping %d bytes%n",
                    refContig.getSequenceName(), refFile.getDistanceBetweenContigs("chr1", refContig.getSequenceName()));
        }
        */

        String lastContig = null;
        List<Double> timings = new ArrayList<Double>();
        for ( SAMSequenceRecord startContig : refFile.getSequenceDictionary().getSequences() ) {
            final String startContigName = startContig.getSequenceName();
            for ( SAMSequenceRecord targetContig : refFile.getSequenceDictionary().getSequences() ) {
                refFile.seekToContig(startContigName, true);
                System.out.printf("Seeking: current=%s, target=%s%n", startContigName, targetContig.getSequenceName());
                long lastTime = System.currentTimeMillis();
                final boolean success = refFile.seekToContig(targetContig.getSequenceName(), true);
                long curTime = System.currentTimeMillis();
                final double elapsed = (curTime - lastTime) / 1000.0;
                timings.add(elapsed);
                System.out.printf("  -> Elapsed time %.2f, averaging %.2f sec / seek for %d seeks%n",
                        elapsed, Utils.averageDouble(timings), timings.size());

                if ( ! success ) {
                    System.out.printf("Failured to seek to %s from %s%n", targetContig.getSequenceName(), lastContig );
                }
                //System.exit(1);
            }
        }
        System.exit(1);

        // code for randomly sampling the seeks
//        Random rnd = new Random();
//        String lastContig = null;
//        List<Double> timings = new ArrayList<Double>();
//        final int N_SAMPLES = 1000;
//        //try { refFile.seekToContig("chr3"); } catch ( IOException e ) {}
//        for ( int i = 0; i < N_SAMPLES; i++ ) {
//            final int nextIndex = rnd.nextInt(refContigs.size());
//            String nextContig = refFile.getSequenceDictionary().getSequence(nextIndex).getSequenceName();
//            //nextContig = "chr2";
//            try {
//                System.out.printf("Seeking: current=%s, target=%s%n", refFile.getContigName(), nextContig);
//                long lastTime = System.currentTimeMillis();
//                final boolean success = refFile.seekToContig(nextContig, true);
//                long curTime = System.currentTimeMillis();
//                final double elapsed = (curTime - lastTime) / 1000.0;
//                timings.add(elapsed);
//                System.out.printf("  -> Elapsed time %.2f, averaging %.2f sec / seek for %d seeks%n",
//                        elapsed, Utils.averageDouble(timings), timings.size());
//
//                if ( ! success ) {
//                    System.out.printf("Failured to seek to %s from %s%n", nextContig, lastContig );
//                }
//                //System.exit(1);
//            } catch ( IOException e ) {
//                System.out.printf("Failured to seek to %s from %s%n", nextContig, lastContig );
//                e.printStackTrace();
//            }
//
//            lastContig = nextContig;
//        }
//        System.exit(1);

/*
        final String targetChr = "chr10";
        try {
            refFile.seekToContig(targetChr);
        } catch ( IOException e ){
            System.out.printf("Failured to seek to %s%n", targetChr);
            e.printStackTrace();            
        }
        System.exit(1);
*/

        //List<Double> timings = new ArrayList<Double>();
        final long startTime = System.currentTimeMillis();
        long lastTime = System.currentTimeMillis();

        int i = 0;
        String prevNextContigName = null;
        System.out.printf("Walking reference sequence:%n");
        for ( SAMSequenceRecord refContig: refContigs ) {
            long curTime = System.currentTimeMillis();
            ReferenceSequence contig = refFile.nextSequence();
            final double elapsed = (curTime - lastTime) / 1000.0;
            timings.add(elapsed);
            System.out.printf("%2d : expected %s contig, found %s with next of %s after %.2f seconds, average is %.2f%n", i,
                    refContig.getSequenceName(), contig.getName(), refFile.getNextContigName(), elapsed, Utils.averageDouble(timings));
            if ( prevNextContigName != null && contig.getName() != null && ! prevNextContigName.equals(contig.getName()) )
                throw new RuntimeIOException(String.format("Unexpected contig ordering %s was expected next, but I found %s?",
                        prevNextContigName, contig.getName()));

            prevNextContigName = refFile.getNextContigName();
            lastTime = curTime;
            i++;

            System.out.printf("  Traversing from chr1 to %s would require jumping %d bytes%n",
                    contig.getName(), refFile.getDistanceBetweenContigs("chr1", contig.getName()));
        }
    }
}

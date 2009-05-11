package org.broadinstitute.sting.gatk;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.GATKArgumentCollection;
import org.junit.After;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * User: aaron
 * Date: May 7, 2009
 * Time: 1:12:58 PM
 *
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT 
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 *
 */


/**
 * @author aaron
 * @version 1.0
 * @date May 7, 2009
 * <p/>
 * Class GATKArgumentCollection
 * <p/>
 * A descriptions should go here. Blame aaron if it's missing.
 */
public class GATKArgumentCollectionTest extends BaseTest {

    // our collection of arguments
    private GATKArgumentCollection collect;

    // where to write our xml file
    private String xmlFileLoc = "testfile.xml";

    /** setup our test */
    @Before
    public void setup() {
        collect = new GATKArgumentCollection();
    }

    /** destroy the temp file */
    @After
    public void takedown() {
        File f = new File(xmlFileLoc);
        if (f.exists()) {
            f.delete();
        }
    }

    private void setupCollection() {
        // parameters and their defaults
        Map<String, String> wArgs = new HashMap<String, String>();
        wArgs.put("wArgType1", "Arg1");
        wArgs.put("wArgType2", "Arg2");
        wArgs.put("wArgType3", "Arg3");
        collect.walkerArgs = wArgs;

        List<File> input = new ArrayList<File>();
        input.add(new File("test.file"));
        collect.samFiles = input;
        collect.maximumReads = "-1";
        collect.strictnessLevel = "strict";
        collect.referenceFile = new File("referenceFile".toLowerCase());
        collect.genomeRegion = "genomeRegion".toLowerCase();
        collect.analysisName = "analysisName".toLowerCase();
        collect.DBSNPFile = "DBSNPFile".toLowerCase();
        collect.HAPMAPFile = "HAPMAPFile".toLowerCase();
        collect.HAPMAPChipFile = "HAPMAPChipFile".toLowerCase();
        collect.enabledThreadedIO = true;
        collect.unsafe = false;
        collect.maximumReadSorts = "maximumReadSorts".toLowerCase();
        collect.downsampleFraction = "downsampleFraction".toLowerCase();
        collect.downsampleCoverage = "downsampleCoverage".toLowerCase();
        collect.intervalsFile = "intervalsFile".toLowerCase();
        collect.walkAllLoci = true;
        collect.disableThreading = false;
        collect.outFileName = "outFileName".toLowerCase();
        collect.errFileName = "errFileName".toLowerCase();
        collect.outErrFileName = "outErrFileName".toLowerCase();
        collect.numberOfThreads = 1;

        // make some rod bindings up
        ArrayList<String> fakeBindings = new ArrayList<String>();
        fakeBindings.add("Bind1");
        fakeBindings.add("Bind2");
        fakeBindings.add("Bind3");

        collect.RODBindings = fakeBindings;
    }


    /** test the output of an XML file in the arg collection */
    @Test
    public void testOutput() {
        setupCollection();

        GATKArgumentCollection.marshal(collect, xmlFileLoc);
        GATKArgumentCollection collection = GATKArgumentCollection.unmarshal(xmlFileLoc);
        if (!collect.equals(collection)) {
            fail("Collections not equal");
        }
    }


    /** test the output of an XML file in the arg collection */
    @Test
    public void testInput() {
        setupCollection();
        GATKArgumentCollection.marshal(collect, xmlFileLoc);
    }
}

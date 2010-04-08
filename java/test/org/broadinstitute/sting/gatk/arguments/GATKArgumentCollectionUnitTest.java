package org.broadinstitute.sting.gatk.arguments;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.arguments.GATKArgumentCollection;
import org.junit.After;
import static org.junit.Assert.fail;
import org.junit.Before;
import org.junit.Test;

import java.io.File;
import java.util.*;

import net.sf.samtools.SAMFileReader;

/*
 * Copyright (c) 2009 The Broad Institute
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
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/**
 * @author aaron
 * @version 1.0
 * @date May 7, 2009
 * <p/>
 * Class GATKArgumentCollection
 * <p/>
 * Test out the argument collection class
 */
public class GATKArgumentCollectionUnitTest extends BaseTest {

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
        collect.maximumEngineIterations = -1;
        collect.strictnessLevel = SAMFileReader.ValidationStringency.STRICT;
        collect.referenceFile = new File("referenceFile".toLowerCase());
        collect.DBSNPFile = "DBSNPFile".toLowerCase();
        collect.HAPMAPFile = "HAPMAPFile".toLowerCase();
        collect.HAPMAPChipFile = "HAPMAPChipFile".toLowerCase();
        collect.unsafe = ValidationExclusion.TYPE.ALL;
        collect.downsampleFraction = null;
        collect.downsampleCoverage = null;
        collect.intervals = new ArrayList<String>();
        collect.intervals.add("intervals".toLowerCase());
        collect.excludeIntervals = new ArrayList<String>();
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

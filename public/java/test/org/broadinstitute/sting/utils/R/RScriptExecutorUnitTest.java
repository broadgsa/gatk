/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.R;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.walkers.diffengine.DiffElement;
import org.broadinstitute.sting.gatk.walkers.diffengine.DiffEngine;
import org.broadinstitute.sting.gatk.walkers.diffengine.DiffNode;
import org.broadinstitute.sting.gatk.walkers.diffengine.Difference;
import org.broadinstitute.sting.utils.exceptions.UserException;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Basic unit test for RScriptExecutor in reduced reads
 */
public class RScriptExecutorUnitTest extends BaseTest {
    final static String testrscript = "print(\"hello, world\")\n";
    final static String publicRScript = "plot_Tranches.R";
    final static String privateRScript = "variantCallQC.R";

    // --------------------------------------------------------------------------------
    //
    // Difference testing routines
    //
    // --------------------------------------------------------------------------------

    private void testOne(String script, String pathToRscript, String anotherSearchPath, boolean exceptOnError) {
        RScriptExecutor.RScriptArgumentCollection collection =
                new RScriptExecutor.RScriptArgumentCollection();
        if ( pathToRscript != null )
            collection.PATH_TO_RSCRIPT = pathToRscript;
        if ( anotherSearchPath != null ) {
            List<String> x = new ArrayList<String>(collection.PATH_TO_RESOURCES);
            x.add(anotherSearchPath);
            collection.PATH_TO_RESOURCES = x;
        }
        RScriptExecutor executor = new RScriptExecutor(collection, exceptOnError);
        executor.callRScripts(script);
    }

    @Test
    public void testPublic() { testOne(publicRScript, null, null, true); }

    @Test
    public void testPrivate() { testOne(privateRScript, null, null, true); }

    // make sure we don't break finding something in private by adding another directory
    @Test
    public void testPrivateWithAdditionalPath1() { testOne(privateRScript, null, "dist", true); }

    // make sure we don't break finding something in private by adding another directory
    @Test
    public void testPrivateWithAdditionalPath2() { testOne(privateRScript, null, "doesNotExist", true); }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistantScriptException() { testOne("does_not_exist.R", null, null, true); }

    @Test()
    public void testNonExistantScriptNoException() { testOne("does_not_exist.R", null, null, false); }

    @Test(expectedExceptions = UserException.class)
    public void testNonExistantRScriptException() { testOne(publicRScript, "badRScriptValue", null, true); }

    @Test()
    public void testNonExistantRScriptNoException() { testOne(publicRScript, "badRScriptValue", null, false); }

    @Test()
    public void testScriptInNewPath() throws IOException {
        File t = createTempFile("myTestScript", ".R");
        FileUtils.writeStringToFile(t, testrscript);
        testOne(t.getName(), null, t.getParent(), true);
    }
}
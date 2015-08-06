/*
* Copyright 2012-2015 Broad Institute, Inc.
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
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.commandline;

import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collection;

public class RodBindingCollectionUnitTest extends BaseTest {

    private ParsingEngine parsingEngine;
    private Tags mytags;

    private static final String defaultTagString = "VCF";
    private static final String testVCFFileName = privateTestDir + "empty.vcf";
    private static final String testListFileName = createTempListFile("oneVCF", testVCFFileName).getAbsolutePath();

    @BeforeMethod
    public void setUp() {
        parsingEngine = new ParsingEngine(null);
        RodBinding.resetNameCounter();
        mytags = new Tags();
        mytags.addPositionalTag(defaultTagString);
    }

    private class RodBindingCollectionArgProvider {
        @Argument(fullName="input",doc="input",shortName="V")
        public RodBindingCollection<VariantContext> input;
    }

    @Test
    public void testStandardVCF() {
        final String[] commandLine = new String[] {"-V", testVCFFileName};

        parsingEngine.addArgumentSource( RodBindingCollectionArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        final RodBindingCollectionArgProvider argProvider = new RodBindingCollectionArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.input.getRodBindings().iterator().next().getSource(), testVCFFileName, "Argument is not correctly initialized");
    }

    @Test
    public void testList() {
        final String[] commandLine = new String[] {"-V", testListFileName};

        parsingEngine.addArgumentSource(RodBindingCollectionArgProvider.class);
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        final RodBindingCollectionArgProvider argProvider = new RodBindingCollectionArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.input.getRodBindings().iterator().next().getSource(), testVCFFileName, "Argument is not correctly initialized");
    }

    @Test
    public void testDefaultTagsInFile() throws IOException {

        final File testFile = createTempListFile("RodBindingCollectionUnitTest.defaultTags", testVCFFileName);

        ArgumentTypeDescriptor.getRodBindingsCollection(testFile, parsingEngine, VariantContext.class, "foo", mytags, "input");

        final Collection<RodBinding> bindings = parsingEngine.getRodBindings();
        Assert.assertNotNull(bindings);
        Assert.assertEquals(bindings.size(), 1);

        final RodBinding binding = bindings.iterator().next();
        Assert.assertEquals(parsingEngine.getTags(binding), mytags);
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testDuplicateEntriesInFile() throws IOException {

        final File testFile = createTempListFile("RodBindingCollectionUnitTest.variantListWithDuplicates", testVCFFileName, testVCFFileName);

        ArgumentTypeDescriptor.getRodBindingsCollection(testFile, parsingEngine, VariantContext.class, "foo", mytags, "input");
    }

    @Test(expectedExceptions = UserException.BadArgumentValue.class)
    public void testValidateEmptyFile() throws IOException {
        final File testFile = createTempListFile("RodBindingCollectionUnitTest.emptyVCFList");

        ArgumentTypeDescriptor.getRodBindingsCollection(testFile, parsingEngine, VariantContext.class, "foo", mytags, "input");
    }

    @Test
    public void testOverrideTagsInFile() throws IOException {
        final File testFile = createTempListFile("RodBindingCollectionUnitTest.overrideTags", "foo " + testVCFFileName);

        ArgumentTypeDescriptor.getRodBindingsCollection(testFile, parsingEngine, VariantContext.class, "foo", mytags, "input");

        final Collection<RodBinding> bindings = parsingEngine.getRodBindings();
        Assert.assertNotNull(bindings);
        Assert.assertEquals(bindings.size(), 1);

        final RodBinding binding = bindings.iterator().next();
        Assert.assertNotEquals(parsingEngine.getTags(binding), mytags);
    }
}

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

package org.broadinstitute.gatk.engine.io.stubs;

import htsjdk.variant.variantcontext.VariantContext;
import it.unimi.dsi.fastutil.objects.ObjectArrayList;
import htsjdk.samtools.SAMFileWriter;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.engine.GenomeAnalysisEngine;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Collection;


public class ArgumentTypeDescriptorUnitTest extends BaseTest {

    ////////////////////////////////////////////////////////////////////
    // This section tests the functionality of the @Output annotation //
    ////////////////////////////////////////////////////////////////////

    private class ATDTestCommandLineProgram extends CommandLineProgram {
        public int execute() { return 0; }

        @Override
        public Collection<ArgumentTypeDescriptor> getArgumentTypeDescriptors() {
            final GenomeAnalysisEngine engine = new GenomeAnalysisEngine();
            return Arrays.asList( new SAMFileWriterArgumentTypeDescriptor(engine, System.out),
                    new OutputStreamArgumentTypeDescriptor(engine, System.out),
                    new VCFWriterArgumentTypeDescriptor(engine, System.out, null));
        }

        protected abstract class ATDTestOutputArgumentSource {
            public abstract Object getOut();
        }

        protected class OutputRequiredSamArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = true)
            public SAMFileWriter out;
            public Object getOut() { return out; }
        }

        protected class OutputRequiredVcfArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = true)
            public VariantContextWriter out;
            public Object getOut() { return out; }
        }

        protected class OutputRequiredStreamArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = true)
            public PrintStream out;
            public Object getOut() { return out; }
        }

        protected class OutputNotRequiredNoDefaultSamArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = false, defaultToStdout = false)
            public SAMFileWriter out;
            public Object getOut() { return out; }
        }

        protected class OutputNotRequiredNoDefaultVcfArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = false, defaultToStdout = false)
            public VariantContextWriter out;
            public Object getOut() { return out; }
        }

        protected class OutputNotRequiredNoDefaultStreamArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = false, defaultToStdout = false)
            public PrintStream out;
            public Object getOut() { return out; }
        }

        protected class OutputNotRequiredSamArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = false)
            public SAMFileWriter out;
            public Object getOut() { return out; }
        }

        protected class OutputNotRequiredVcfArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = false)
            public VariantContextWriter out;
            public Object getOut() { return out; }
        }

        protected class OutputNotRequiredStreamArgumentSource extends ATDTestOutputArgumentSource {
            @Output(shortName="o", doc="output file", required = false)
            public PrintStream out;
            public Object getOut() { return out; }
        }
    }

    @DataProvider(name = "OutputProvider")
    public Object[][] OutputProvider() {

        ObjectArrayList<Object[]> tests = new ObjectArrayList<Object[]>();

        final ATDTestCommandLineProgram clp = new ATDTestCommandLineProgram();

        for ( final Object obj : Arrays.asList(clp.new OutputRequiredSamArgumentSource(), clp.new OutputRequiredVcfArgumentSource(), clp.new OutputRequiredStreamArgumentSource()) ) {
            for ( final boolean provided : Arrays.asList(true, false) ) {
                tests.add(new Object[]{obj, true, true, provided});
            }
        }

        for ( final Object obj : Arrays.asList(clp.new OutputNotRequiredSamArgumentSource(), clp.new OutputNotRequiredVcfArgumentSource(), clp.new OutputNotRequiredStreamArgumentSource()) ) {
            for ( final boolean provided : Arrays.asList(true, false) ) {
                tests.add(new Object[]{obj, false, true, provided});
            }
        }

        for ( final Object obj : Arrays.asList(clp.new OutputNotRequiredNoDefaultSamArgumentSource(), clp.new OutputNotRequiredNoDefaultVcfArgumentSource(), clp.new OutputNotRequiredNoDefaultStreamArgumentSource()) ) {
            for ( final boolean provided : Arrays.asList(true, false) ) {
                tests.add(new Object[]{obj, false, false, provided});
            }
        }

        return tests.toArray(new Object[][]{});
    }

    @Test(dataProvider = "OutputProvider")
    public void testOutput(final ATDTestCommandLineProgram.ATDTestOutputArgumentSource argumentSource, final boolean required, final boolean hasDefault, final boolean provided) {

        final ParsingEngine parser = new ParsingEngine(new ATDTestCommandLineProgram());
        parser.addArgumentSource(argumentSource.getClass());
        parser.parse(provided ? new String[] {"out", "foo"} : new String[] {});

        try {
            parser.loadArgumentsIntoObject(argumentSource);

            if ( !provided && (required || !hasDefault) )
                Assert.assertEquals(argumentSource.getOut(), null);
            else if ( !provided )
                Assert.assertNotEquals(argumentSource.getOut(), null);
            else if ( argumentSource.getOut() == null || !(argumentSource.getOut() instanceof SAMFileWriterStub) ) // can't test this one case
                Assert.assertEquals(!provided, outputIsStdout(argumentSource.getOut()));

        } catch (Exception e) {
            throw new ReviewedGATKException(e.getMessage());
        }
    }

    @Test
    public void testRodBindingsCollection() {

        final ParsingEngine parser = new ParsingEngine(new ATDTestCommandLineProgram());

        //A list file containing a single VCF
        final File listFile = createTempListFile("oneVCF", privateTestDir + "empty.vcf");

        try {
            Object result = ArgumentTypeDescriptor.getRodBindingsCollection(listFile,
                    parser,
                    VariantContext.class,
                    "variant",
                    new Tags(),
                    "variantTest");
            if (!(result instanceof RodBindingCollection))
                throw new ReviewedGATKException("getRodBindingsCollection did not return a RodBindingCollection");
            RodBindingCollection<?> rbc = (RodBindingCollection) result;

            Assert.assertEquals(rbc.getType(), VariantContext.class);
            Assert.assertEquals(rbc.getRodBindings().size(), 1);

        } catch (IOException e) {
            throw new ReviewedGATKException(e.getMessage(), e);
        }

        //The same file, now with an extra blank line
        final File listFileWithBlank = createTempListFile("oneVCFwithBlankLine", privateTestDir + "empty.vcf", "");
        try {
            Object result = ArgumentTypeDescriptor.getRodBindingsCollection(listFileWithBlank,
                    parser,
                    VariantContext.class,
                    "variant",
                    new Tags(),
                    "variantTest");
            if (!(result instanceof RodBindingCollection))
                throw new ReviewedGATKException("getRodBindingsCollection did not return a RodBindingCollection");
            RodBindingCollection<?> rbc = (RodBindingCollection) result;

            Assert.assertEquals(rbc.getType(), VariantContext.class);
            Assert.assertEquals(rbc.getRodBindings().size(), 1);

        } catch (IOException e) {
            throw new ReviewedGATKException(e.getMessage(), e);
        }
    }

    private static boolean outputIsStdout(final Object out) {
        if ( out == null ) {
            return false;
        } else if ( out instanceof SAMFileWriterStub ) {
            return ((SAMFileWriterStub)out).getOutputStream() != System.out;
        } else if ( out instanceof VariantContextWriterStub ) {
            return ((VariantContextWriterStub)out).getOutputStream() == System.out;
        } else if ( out instanceof OutputStreamStub ) {
            return ((OutputStreamStub)out).getOutputStream() == System.out;
        }
        return false;
    }

}
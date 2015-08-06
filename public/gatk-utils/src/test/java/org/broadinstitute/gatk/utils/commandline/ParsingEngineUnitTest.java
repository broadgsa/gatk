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

import org.apache.commons.io.FileUtils;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import htsjdk.variant.variantcontext.VariantContext;
import org.testng.Assert;
import org.broadinstitute.gatk.utils.BaseTest;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.EnumSet;
import java.util.Set;

/**
 * Test suite for the parsing engine.
 */
public class ParsingEngineUnitTest extends BaseTest {
    /** we absolutely cannot have this file existing, or we'll fail the UnitTest */
    private final static String NON_EXISTANT_FILENAME_VCF = "this_file_should_not_exist_on_disk_123456789.vcf";
    private ParsingEngine parsingEngine;

    @BeforeMethod
    public void setUp() {
        parsingEngine = new ParsingEngine(null);
        RodBinding.resetNameCounter();
    }

    private class InputFileArgProvider {
        @Argument(fullName="input_file",doc="input file",shortName="I")
        public String inputFile;
    }

    @Test
    public void shortNameArgumentTest() {
        final String[] commandLine = new String[] {"-I","na12878.bam"};

        parsingEngine.addArgumentSource( InputFileArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.inputFile,"na12878.bam","Argument is not correctly initialized");
    }

    @Test
    public void multiCharShortNameArgumentTest() {
        final String[] commandLine = new String[] {"-out","out.txt"};

        parsingEngine.addArgumentSource( MultiCharShortNameArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        MultiCharShortNameArgProvider argProvider = new MultiCharShortNameArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.outputFile,"out.txt","Argument is not correctly initialized");
    }


    private class MultiCharShortNameArgProvider {
        @Argument(shortName="out", doc="output file")
        public String outputFile;
    }

    @Test
    public void longNameArgumentTest() {
        final String[] commandLine = new String[] {"--input_file", "na12878.bam"};

        parsingEngine.addArgumentSource( InputFileArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.inputFile,"na12878.bam","Argument is not correctly initialized");
    }

    @Test
    public void extraWhitespaceTest() {
        final String[] commandLine = new String[] {"  --input_file ", "na12878.bam"};

        parsingEngine.addArgumentSource( InputFileArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.inputFile,"na12878.bam","Argument is not correctly initialized");
    }

    @Test
    public void primitiveArgumentTest() {
        final String[] commandLine = new String[] {"--foo", "5"};

        parsingEngine.addArgumentSource( PrimitiveArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        PrimitiveArgProvider argProvider = new PrimitiveArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.foo, 5, "Argument is not correctly initialized");
    }

    @Test(expectedExceptions=InvalidArgumentValueException.class)
    public void primitiveArgumentNoValueTest() {
        final String[] commandLine = new String[] {"--foo"};

        parsingEngine.addArgumentSource( PrimitiveArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        PrimitiveArgProvider argProvider = new PrimitiveArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.foo, 5, "Argument is not correctly initialized");
    }

    private class PrimitiveArgProvider {
        @Argument(doc="simple integer")
        int foo;
    }

    @Test
    public void flagTest() {
        final String[] commandLine = new String[] {"--all_loci"};

        parsingEngine.addArgumentSource( AllLociArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        AllLociArgProvider argProvider = new AllLociArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertTrue(argProvider.allLoci,"Argument is not correctly initialized");
    }

    private class AllLociArgProvider {
        @Argument(fullName="all_loci",shortName="A", doc="all loci")
        public boolean allLoci = false;
    }

    @Test
    public void arrayTest() {
        final String[] commandLine = new String[] {"-I", "foo.txt", "--input_file", "bar.txt"};

        parsingEngine.addArgumentSource( MultiValueArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        MultiValueArgProvider argProvider = new MultiValueArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.inputFile.length, 2, "Argument array is of incorrect length");
        Assert.assertEquals(argProvider.inputFile[0],"foo.txt","1st filename is incorrect");
        Assert.assertEquals(argProvider.inputFile[1],"bar.txt","2nd filename is incorrect");
    }

    private class MultiValueArgProvider {
        @Argument(fullName="input_file",shortName="I", doc="input file")
        public String[] inputFile;
    }

    @Test
    public void enumTest() {
        final String[] commandLine = new String[] {  "--test_enum", "TWO" };

        parsingEngine.addArgumentSource( EnumArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        EnumArgProvider argProvider = new EnumArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.testEnum, TestEnum.TWO, "Enum value is not correct");
    }

    @Test
    public void enumMixedCaseTest() {
        final String[] commandLine = new String[] {  "--test_enum", "oNe" };

        parsingEngine.addArgumentSource( EnumArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        EnumArgProvider argProvider = new EnumArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.testEnum, TestEnum.ONE, "Enum value is not correct");
    }

    @Test
    public void enumDefaultTest() {
        final String[] commandLine = new String[] {};

        parsingEngine.addArgumentSource( EnumArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        EnumArgProvider argProvider = new EnumArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.testEnum, TestEnum.THREE, "Enum value is not correct");
    }

    public enum TestEnum { ONE, TWO, THREE }

    private class EnumArgProvider {
        @Argument(fullName="test_enum",shortName="ti",doc="test enum",required=false)
        public TestEnum testEnum = TestEnum.THREE;
    }

    @Test
    public void typedCollectionTest() {
        final String[] commandLine = new String[] { "-N","2","-N","4","-N","6","-N","8","-N","10" };

        parsingEngine.addArgumentSource( IntegerListArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        IntegerListArgProvider argProvider = new IntegerListArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertNotNull(argProvider.integers, "Argument array is null");
        Assert.assertEquals(argProvider.integers.size(), 5, "Argument array is of incorrect length");
        Assert.assertEquals(argProvider.integers.get(0).intValue(), 2, "1st integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(1).intValue(), 4, "2nd integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(2).intValue(), 6, "3rd integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(3).intValue(), 8, "4th integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(4).intValue(), 10, "5th integer is incorrect");
    }

    private class IntegerListArgProvider {
        @Argument(fullName="integer_list",shortName="N",doc="integer list")
        public List<Integer> integers;
    }

    @Test
    public void untypedCollectionTest() {
        final String[] commandLine = new String[] { "-N","2","-N","4","-N","6","-N","8","-N","10" };

        parsingEngine.addArgumentSource( UntypedListArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        UntypedListArgProvider argProvider = new UntypedListArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertNotNull(argProvider.integers, "Argument array is null");
        Assert.assertEquals(argProvider.integers.size(), 5, "Argument array is of incorrect length");
        Assert.assertEquals(argProvider.integers.get(0), "2", "1st integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(1), "4", "2nd integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(2), "6", "3rd integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(3), "8", "4th integer is incorrect");
        Assert.assertEquals(argProvider.integers.get(4), "10", "5th integer is incorrect");
    }

    private class UntypedListArgProvider {
        @Argument(fullName="untyped_list",shortName="N", doc="untyped list")
        public List integers;
    }

    @Test(expectedExceptions=MissingArgumentException.class)
    public void requiredArgTest() {
        final String[] commandLine = new String[0];

        parsingEngine.addArgumentSource( RequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    private class RequiredArgProvider {
        @Argument(required=true,doc="value")
        public Integer value;
    }

    @Test
    public void defaultValueTest() {
        // First try getting the default.
        String[] commandLine = new String[0];

        parsingEngine.addArgumentSource( DefaultValueArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        DefaultValueArgProvider argProvider = new DefaultValueArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.value.intValue(), 42, "Default value is not correctly initialized");

        // Then try to override it.
        commandLine = new String[] { "--value", "27" };

        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.value.intValue(), 27, "Default value is not correctly initialized");
    }

    private class DefaultValueArgProvider {
        @Argument(doc="value",required=false)
        public Integer value = 42;
    }

    @Test
    public void disableValidationOfRequiredArgTest() {
        final String[] commandLine = new String[0];

        parsingEngine.addArgumentSource( RequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate( EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument) );

        RequiredArgProvider argProvider = new RequiredArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider );

        Assert.assertNull(argProvider.value, "Value should have remain unset");
    }

    @Test
    public void unrequiredArgTest() {
        final String[] commandLine = new String[0];

        parsingEngine.addArgumentSource( UnrequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        UnrequiredArgProvider argProvider = new UnrequiredArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertNull(argProvider.value, "Value was unrequired and unspecified; contents should be null");
    }

    private class UnrequiredArgProvider {
        @Argument(required=false,doc="unrequired value")
        public Integer value;
    }

    @Test(expectedExceptions=InvalidArgumentException.class)
    public void invalidArgTest() {
        final String[] commandLine = new String[] { "--foo" };

        parsingEngine.addArgumentSource( UnrequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test(expectedExceptions= ReviewedGATKException.class)
    public void duplicateLongNameTest() {
        parsingEngine.addArgumentSource( DuplicateLongNameProvider.class );
    }

    private class DuplicateLongNameProvider {
        @Argument(fullName="myarg",doc="my arg")
        public Integer foo;

        @Argument(fullName="myarg", doc="my arg")
        public Integer bar;
    }

    @Test(expectedExceptions= ReviewedGATKException.class)
    public void duplicateShortNameTest() {
        parsingEngine.addArgumentSource( DuplicateShortNameProvider.class );
    }


    private class DuplicateShortNameProvider {
        @Argument(shortName="myarg", doc="my arg")
        public Integer foo;

        @Argument(shortName="myarg", doc="my arg")
        public Integer bar;
    }

    @Test(expectedExceptions=UnmatchedArgumentException.class)
    public void missingArgumentNameTest() {
        final String[] commandLine = new String[] {"foo.txt"};

        parsingEngine.addArgumentSource( NoArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    private class NoArgProvider {

    }

    @Test(expectedExceptions=UnmatchedArgumentException.class)
    public void extraValueTest() {
        final String[] commandLine = new String[] {"-I", "foo.txt", "bar.txt"};

        parsingEngine.addArgumentSource( InputFileArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test(expectedExceptions=MissingArgumentException.class)
    public void multipleInvalidArgTest() {
        final String[] commandLine = new String[] {"-N1", "-N2", "-N3"};

        parsingEngine.addArgumentSource( RequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test(expectedExceptions=TooManyValuesForArgumentException.class)
    public void invalidArgCountTest() {
        final String[] commandLine = new String[] {"--value","1","--value","2","--value","3"};

        parsingEngine.addArgumentSource( RequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test
    public void packageProtectedArgTest() {
        final String[] commandLine = new String[] {"--foo", "1"};

        parsingEngine.addArgumentSource( PackageProtectedArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        PackageProtectedArgProvider argProvider = new PackageProtectedArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);

        Assert.assertEquals(argProvider.foo.intValue(), 1, "Argument is not correctly initialized");
    }

    private class PackageProtectedArgProvider {
        @Argument(doc="foo")
        Integer foo;
    }

    @Test
    public void derivedArgTest() {
        final String[] commandLine = new String[] {"--bar", "5"};

        parsingEngine.addArgumentSource( DerivedArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        DerivedArgProvider argProvider = new DerivedArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);

        Assert.assertEquals(argProvider.bar.intValue(), 5, "Argument is not correctly initialized");
    }

    private class DerivedArgProvider extends BaseArgProvider {
    }

    private class BaseArgProvider {
        @Argument(doc="bar")
        public Integer bar;
    }

    @Test
    public void correctDefaultArgNameTest() {
        parsingEngine.addArgumentSource( CamelCaseArgProvider.class );

        DefinitionMatcher matcher = ArgumentDefinitions.FullNameDefinitionMatcher;
        ArgumentDefinition definition = parsingEngine.argumentDefinitions.findArgumentDefinition("myarg", matcher);

        Assert.assertNotNull(definition, "Invalid default argument name assigned");
    }

    @SuppressWarnings("unused")
    private class CamelCaseArgProvider {
        @Argument(doc="my arg")
        Integer myArg;
    }

    @Test(expectedExceptions=UnmatchedArgumentException.class)
    public void booleanWithParameterTest() {
        final String[] commandLine = new String[] {"--mybool", "true"};

        parsingEngine.addArgumentSource( BooleanArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @SuppressWarnings("unused")
    private class BooleanArgProvider {
        @Argument(doc="my bool")
        boolean myBool;
    }

    @Test
    public void validParseForAnalysisTypeTest() {
        final String[] commandLine = new String[] {"--analysis_type", "Pileup" };

        parsingEngine.addArgumentSource( AnalysisTypeArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate( EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument) );

        AnalysisTypeArgProvider argProvider = new AnalysisTypeArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.Analysis_Name,"Pileup","Argument is not correctly initialized");
    }

    private class AnalysisTypeArgProvider {
        @Argument(fullName="analysis_type", shortName="T", doc="Type of analysis to run")
        public String Analysis_Name = null;
    }

    @Test(expectedExceptions=TooManyValuesForArgumentException.class)
    public void invalidParseForAnalysisTypeTest() {
        final String[] commandLine = new String[] {"--analysis_type", "Pileup", "-T", "CountReads" };

        parsingEngine.addArgumentSource( AnalysisTypeArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate( EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument) );
    }

    @Test(expectedExceptions=ArgumentsAreMutuallyExclusiveException.class)
    public void mutuallyExclusiveArgumentsTest() {
        // Passing only foo should work fine...
        String[] commandLine = new String[] {"--foo","5"};

        parsingEngine.addArgumentSource( MutuallyExclusiveArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        MutuallyExclusiveArgProvider argProvider = new MutuallyExclusiveArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.foo.intValue(), 5, "Argument is not correctly initialized");

        // But when foo and bar come together, danger!
        commandLine = new String[] {"--foo","5","--bar","6"};

        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @SuppressWarnings("unused")
    private class MutuallyExclusiveArgProvider {
        @Argument(doc="foo",exclusiveOf="bar")
        Integer foo;

        @Argument(doc="bar",required=false)
        Integer bar;
    }

    @Test(expectedExceptions=InvalidArgumentValueException.class)
    public void argumentValidationTest() {
        // Passing only foo should work fine...
        String[] commandLine = new String[] {"--value","521"};

        parsingEngine.addArgumentSource( ValidatingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ValidatingArgProvider argProvider = new ValidatingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.value.intValue(), 521, "Argument is not correctly initialized");

        // Try some invalid arguments
        commandLine = new String[] {"--value","foo"};
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    private class ValidatingArgProvider {
        @Argument(doc="value",validation="\\d+")
        Integer value;
    }

    @Test
    public void argumentCollectionTest() {
        String[] commandLine = new String[] { "--value", "5" };

        parsingEngine.addArgumentSource( ArgumentCollectionProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ArgumentCollectionProvider argProvider = new ArgumentCollectionProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);

        Assert.assertEquals(argProvider.rap.value.intValue(), 5, "Argument is not correctly initialized");
    }

    private class ArgumentCollectionProvider {
        @ArgumentCollection
        RequiredArgProvider rap = new RequiredArgProvider();
    }

    @Test(expectedExceptions= ReviewedGATKException.class)
    public void multipleArgumentCollectionTest() {
        parsingEngine.addArgumentSource( MultipleArgumentCollectionProvider.class );
    }

    @SuppressWarnings("unused")
    private class MultipleArgumentCollectionProvider {
        @ArgumentCollection
        RequiredArgProvider rap1 = new RequiredArgProvider();
        @ArgumentCollection
        RequiredArgProvider rap2 = new RequiredArgProvider();
    }

    // --------------------------------------------------------------------------------
    //
    // Tests of the RodBinding<T> system
    //
    // --------------------------------------------------------------------------------

    private class SingleRodBindingArgProvider {
        @Input(fullName="binding", shortName="V", required=true)
        public RodBinding<Feature> binding;
    }

    @Test
    public void basicRodBindingArgumentTest() {
        final String[] commandLine = new String[] {"-V:vcf",NON_EXISTANT_FILENAME_VCF};

        parsingEngine.addArgumentSource( SingleRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        SingleRodBindingArgProvider argProvider = new SingleRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.binding.getName(), "binding", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), NON_EXISTANT_FILENAME_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), Feature.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.isBound(), true, "Bound() isn't returning its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");
    }

    private class ShortNameOnlyRodBindingArgProvider {
        @Input(shortName="short", required=false)
        public RodBinding<Feature> binding; // = RodBinding.makeUnbound(Feature.class);
    }

    @Test
    public void shortNameOnlyRodBindingArgumentTest() {
        final String[] commandLine = new String[] {"-short:vcf",NON_EXISTANT_FILENAME_VCF};

        parsingEngine.addArgumentSource( ShortNameOnlyRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ShortNameOnlyRodBindingArgProvider argProvider = new ShortNameOnlyRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.binding.getName(), "binding", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), NON_EXISTANT_FILENAME_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), Feature.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.isBound(), true, "Bound() isn't returning its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");
    }

    private class OptionalRodBindingArgProvider {
        @Input(fullName="binding", shortName="V", required=false)
        public RodBinding<Feature> binding;

        @Input(fullName="bindingNull", shortName="VN", required=false)
        public RodBinding<VariantContext> bindingNull = null;
    }

    @Test
    public void optionalRodBindingArgumentTest() {
        final String[] commandLine = new String[] {};

        parsingEngine.addArgumentSource( OptionalRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        OptionalRodBindingArgProvider argProvider = new OptionalRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertNotNull(argProvider.binding, "Default value not applied corrected to RodBinding");
        Assert.assertEquals(argProvider.binding.getName(), RodBinding.UNBOUND_VARIABLE_NAME, "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), RodBinding.UNBOUND_SOURCE, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), Feature.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.isBound(), false, "Bound() isn't returning its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 0, "Tags aren't correctly set");

        Assert.assertNotNull(argProvider.bindingNull, "Default value not applied corrected to RodBinding");
        Assert.assertEquals(argProvider.bindingNull.getName(), RodBinding.UNBOUND_VARIABLE_NAME, "Name isn't set properly");
        Assert.assertEquals(argProvider.bindingNull.getSource(), RodBinding.UNBOUND_SOURCE, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.bindingNull.getType(), VariantContext.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.bindingNull.isBound(), false, "Bound() isn't returning its expected value");
        Assert.assertEquals(argProvider.bindingNull.getTags().getPositionalTags().size(), 0, "Tags aren't correctly set");
    }

    @Test(expectedExceptions = UserException.class)
    public void rodBindingArgumentTestMissingType() {
        final String[] commandLine = new String[] {"-V",NON_EXISTANT_FILENAME_VCF};

        parsingEngine.addArgumentSource( SingleRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        SingleRodBindingArgProvider argProvider = new SingleRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);
    }

    @Test(expectedExceptions = UserException.class)
    public void rodBindingArgumentTestTooManyTags() {
        final String[] commandLine = new String[] {"-V:x,y,z",NON_EXISTANT_FILENAME_VCF};

        parsingEngine.addArgumentSource( SingleRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        SingleRodBindingArgProvider argProvider = new SingleRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);
    }

    private class VariantContextRodBindingArgProvider {
        @Input(fullName = "binding", shortName="V")
        public RodBinding<VariantContext> binding;
    }

    @Test
    public void variantContextBindingArgumentTest() {
        final String[] commandLine = new String[] {"-V:vcf",NON_EXISTANT_FILENAME_VCF};

        parsingEngine.addArgumentSource( VariantContextRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariantContextRodBindingArgProvider argProvider = new VariantContextRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.binding.getName(), "binding", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), NON_EXISTANT_FILENAME_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), VariantContext.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");
    }

    private class ListRodBindingArgProvider {
        @Input(fullName = "binding", shortName="V", required=false)
        public List<RodBinding<Feature>> bindings;
    }

    @Test
    public void listRodBindingArgumentTest() {
        final String[] commandLine = new String[] {"-V:vcf",NON_EXISTANT_FILENAME_VCF};

        parsingEngine.addArgumentSource( ListRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ListRodBindingArgProvider argProvider = new ListRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.bindings.size(), 1, "Unexpected number of bindings");
        RodBinding<Feature> binding = argProvider.bindings.get(0);
        Assert.assertEquals(binding.getName(), "binding", "Name isn't set properly");
        Assert.assertEquals(binding.getSource(), NON_EXISTANT_FILENAME_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(binding.getType(), Feature.class, "Type isn't set to its expected value");
        Assert.assertEquals(binding.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");
    }

    @Test
    public void listRodBindingArgumentTest2Args() {
        final String[] commandLine = new String[] {"-V:vcf",NON_EXISTANT_FILENAME_VCF, "-V:vcf", "bar.vcf"};

        parsingEngine.addArgumentSource( ListRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ListRodBindingArgProvider argProvider = new ListRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.bindings.size(), 2, "Unexpected number of bindings");

        RodBinding<Feature> binding = argProvider.bindings.get(0);
        Assert.assertEquals(binding.getName(), "binding", "Name isn't set properly");
        Assert.assertEquals(binding.getSource(), NON_EXISTANT_FILENAME_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(binding.getType(), Feature.class, "Type isn't set to its expected value");
        Assert.assertEquals(binding.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");

        RodBinding<Feature> binding2 = argProvider.bindings.get(1);
        Assert.assertEquals(binding2.getName(), "binding2", "Name isn't set properly");
        Assert.assertEquals(binding2.getSource(), "bar.vcf", "Source isn't set to its expected value");
        Assert.assertEquals(binding2.getType(), Feature.class, "Type isn't set to its expected value");
        Assert.assertEquals(binding2.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");
    }

    @Test
    public void listRodBindingArgumentTest0Args() {
        final String[] commandLine = new String[] {};

        parsingEngine.addArgumentSource( ListRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ListRodBindingArgProvider argProvider = new ListRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertNull(argProvider.bindings, "Bindings were not null");
    }

    @Test
    public void listRodBindingArgumentTestExplicitlyNamed() {
        final String[] commandLine = new String[] {"-V:foo,vcf",NON_EXISTANT_FILENAME_VCF, "-V:foo,vcf", "bar.vcf"};

        parsingEngine.addArgumentSource( ListRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ListRodBindingArgProvider argProvider = new ListRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.bindings.size(), 2, "Unexpected number of bindings");
        Assert.assertEquals(argProvider.bindings.get(0).getName(), "foo", "Name isn't set properly");
        Assert.assertEquals(argProvider.bindings.get(1).getName(), "foo2", "Name isn't set properly");
    }

    private final static String HISEQ_VCF = privateTestDir + "HiSeq.10000.vcf";
    private final static String TRANCHES_FILE = privateTestDir + "tranches.6.txt";

    @Test
    public void variantContextBindingTestDynamicTyping1() {
        final String[] commandLine = new String[] {"-V", HISEQ_VCF};

        parsingEngine.addArgumentSource( VariantContextRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariantContextRodBindingArgProvider argProvider = new VariantContextRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.binding.getName(), "binding", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), HISEQ_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), VariantContext.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 0, "Tags aren't correctly set");
    }

    @Test
    public void variantContextBindingTestDynamicTypingNameAsSingleArgument() {
        final String[] commandLine = new String[] {"-V:name", HISEQ_VCF};

        parsingEngine.addArgumentSource( VariantContextRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariantContextRodBindingArgProvider argProvider = new VariantContextRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.binding.getName(), "name", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), HISEQ_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), VariantContext.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 1, "Tags aren't correctly set");
    }

    @Test()
    public void variantContextBindingTestDynamicTypingTwoTagsPassing() {
        final String[] commandLine = new String[] {"-V:name,vcf", HISEQ_VCF};

        parsingEngine.addArgumentSource( VariantContextRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariantContextRodBindingArgProvider argProvider = new VariantContextRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals(argProvider.binding.getName(), "name", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), HISEQ_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), VariantContext.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 2, "Tags aren't correctly set");
    }

    @Test()
    public void variantContextBindingTestDynamicTypingTwoTagsCausingTypeFailure() {
        final String[] commandLine = new String[] {"-V:name,beagle", HISEQ_VCF};

        parsingEngine.addArgumentSource( VariantContextRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariantContextRodBindingArgProvider argProvider = new VariantContextRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);

        Assert.assertEquals(argProvider.binding.getName(), "name", "Name isn't set properly");
        Assert.assertEquals(argProvider.binding.getSource(), HISEQ_VCF, "Source isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getType(), VariantContext.class, "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getTribbleType(), "beagle", "Type isn't set to its expected value");
        Assert.assertEquals(argProvider.binding.getTags().getPositionalTags().size(), 2, "Tags aren't correctly set");
    }

    @Test(expectedExceptions = UserException.class)
    public void variantContextBindingTestDynamicTypingUnknownTribbleType() {
        final String[] commandLine = new String[] {"-V", TRANCHES_FILE};

        parsingEngine.addArgumentSource( VariantContextRodBindingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariantContextRodBindingArgProvider argProvider = new VariantContextRodBindingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );
    }

    @Test
    public void argumentListTest() throws IOException {
        File argsFile = BaseTest.createTempListFile("args.", "-I na12878.bam");
        try {
            final String[] commandLine = new String[] {"-args", argsFile.getPath()};
            parsingEngine.addArgumentSource(InputFileArgProvider.class);
            parsingEngine.parse(commandLine);
            parsingEngine.validate();

            InputFileArgProvider argProvider = new InputFileArgProvider();
            parsingEngine.loadArgumentsIntoObject(argProvider);

            Assert.assertEquals(argProvider.inputFile, "na12878.bam", "Argument is not correctly initialized");
        } finally {
            FileUtils.deleteQuietly(argsFile);
        }
    }

    @SuppressWarnings("unused")
    private class NumericRangeArgProvider {
        @Argument(fullName = "intWithHardMinAndMax", minValue = 5, maxValue = 10)
        public int intWithHardMinAndMax;

        @Argument(fullName = "intWithHardMin", minValue = 5)
        public int intWithHardMin;

        @Argument(fullName = "intWithHardMax", maxValue = 10)
        public int intWithHardMax;

        @Argument(fullName = "intWithSoftMinAndMax", minRecommendedValue = 5, maxRecommendedValue = 10)
        public int intWithSoftMinAndMax;

        @Argument(fullName = "intWithSoftMin", minRecommendedValue = 5)
        public int intWithSoftMin;

        @Argument(fullName = "intWithSoftMax", maxRecommendedValue = 10)
        public int intWithSoftMax;

        @Argument(fullName = "intWithHardAndSoftMinAndMax", minValue = 5, minRecommendedValue = 7, maxValue = 10, maxRecommendedValue = 9)
        public int intWithHardAndSoftMinAndMax;

        @Argument(fullName = "intWithHardAndSoftMin", minValue = 5, minRecommendedValue = 7)
        public int intWithHardAndSoftMin;

        @Argument(fullName = "intWithHardAndSoftMax", maxValue = 10, maxRecommendedValue = 8)
        public int intWithHardAndSoftMax;

        @Argument(fullName = "intWithHardMinAndMaxDefaultOutsideRange", minValue = 5, maxValue = 10)
        public int intWithHardMinAndMaxDefaultOutsideRange = -1;

        @Argument(fullName = "integerWithHardMinAndMax", minValue = 5, maxValue = 10)
        public Integer integerWithHardMinAndMax;

        @Argument(fullName = "byteWithHardMinAndMax", minValue = 5, maxValue = 10)
        public byte byteWithHardMinAndMax;

        @Argument(fullName = "byteWithHardMin", minValue = 5)
        public byte byteWithHardMin;

        @Argument(fullName = "byteWithHardMax", maxValue = 10)
        public byte byteWithHardMax;

        @Argument(fullName = "doubleWithHardMinAndMax", minValue = 5.5, maxValue = 10.0)
        public double doubleWithHardMinAndMax;

        @Argument(fullName = "doubleWithHardMin", minValue = 5.5)
        public double doubleWithHardMin;

        @Argument(fullName = "doubleWithHardMax", maxValue = 10.0)
        public double doubleWithHardMax;
    }

    @DataProvider(name = "NumericRangeConstraintViolationDataProvider")
    public Object[][] numericRangeConstraintViolationDataProvider() {
        return new Object[][] {
                { new String[]{"--intWithHardMinAndMax", "11"} },
                { new String[]{"--intWithHardMinAndMax", "4"} },
                { new String[]{"--intWithHardMin", "4"} },
                { new String[]{"--intWithHardMax", "11"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "11"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "4"} },
                { new String[]{"--intWithHardAndSoftMin", "4"} },
                { new String[]{"--intWithHardAndSoftMax", "11"} },
                { new String[]{"--intWithHardMinAndMaxDefaultOutsideRange", "11"} },
                { new String[]{"--intWithHardMinAndMaxDefaultOutsideRange", "4"} },
                { new String[]{"--integerWithHardMinAndMax", "11"} },
                { new String[]{"--integerWithHardMinAndMax", "4"} },
                { new String[]{"--byteWithHardMinAndMax", "11"} },
                { new String[]{"--byteWithHardMinAndMax", "4"} },
                { new String[]{"--byteWithHardMin", "4"} },
                { new String[]{"--byteWithHardMax", "11"} },
                { new String[]{"--doubleWithHardMinAndMax", "5.4"} },
                { new String[]{"--doubleWithHardMinAndMax", "10.1"} },
                { new String[]{"--doubleWithHardMin", "5.4"} },
                { new String[]{"--doubleWithHardMax", "10.1"} }
        };
    }

    @Test(dataProvider = "NumericRangeConstraintViolationDataProvider",
          expectedExceptions = ArgumentValueOutOfRangeException.class)
    public void testNumericRangeWithConstraintViolation( final String[] commandLine ) {
        runNumericArgumentRangeTest(commandLine);
    }

    @DataProvider(name = "NumericRangeWithoutConstraintViolationDataProvider")
    public Object[][] numericRangeWithoutConstraintViolationDataProvider() {
        return new Object[][] {
                { new String[]{"--intWithHardMinAndMax", "10"} },
                { new String[]{"--intWithHardMinAndMax", "5"} },
                { new String[]{"--intWithHardMinAndMax", "7"} },
                { new String[]{"--intWithHardMin", "11"} },
                { new String[]{"--intWithHardMax", "4"} },
                { new String[]{"--intWithSoftMinAndMax", "11"} },
                { new String[]{"--intWithSoftMinAndMax", "4"} },
                { new String[]{"--intWithSoftMin", "4"} },
                { new String[]{"--intWithSoftMax", "11"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "5"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "7"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "8"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "9"} },
                { new String[]{"--intWithHardAndSoftMinAndMax", "10"} },
                { new String[]{"--intWithHardAndSoftMin", "5"} },
                { new String[]{"--intWithHardAndSoftMin", "6"} },
                { new String[]{"--intWithHardAndSoftMin", "7"} },
                { new String[]{"--intWithHardAndSoftMax", "10"} },
                { new String[]{"--intWithHardAndSoftMax", "9"} },
                { new String[]{"--intWithHardAndSoftMax", "8"} },
                { new String[]{"--intWithHardMinAndMaxDefaultOutsideRange", "10"} },
                { new String[]{"--intWithHardMinAndMaxDefaultOutsideRange", "5"} },
                { new String[]{"--intWithHardMinAndMaxDefaultOutsideRange", "7"} },
                { new String[]{"--integerWithHardMinAndMax", "10"} },
                { new String[]{"--integerWithHardMinAndMax", "5"} },
                { new String[]{"--byteWithHardMinAndMax", "10"} },
                { new String[]{"--byteWithHardMinAndMax", "5"} },
                { new String[]{"--byteWithHardMinAndMax", "7"} },
                { new String[]{"--byteWithHardMin", "5"} },
                { new String[]{"--byteWithHardMax", "10"} },
                { new String[]{"--doubleWithHardMinAndMax", "5.5"} },
                { new String[]{"--doubleWithHardMinAndMax", "10.0"} },
                { new String[]{"--doubleWithHardMinAndMax", "7.5"} },
                { new String[]{"--doubleWithHardMin", "5.5"} },
                { new String[]{"--doubleWithHardMin", "15.5"} },
                { new String[]{"--doubleWithHardMax", "10.0"} },
                { new String[]{"--doubleWithHardMax", "7.5"} }
        };
    }

    @Test(dataProvider = "NumericRangeWithoutConstraintViolationDataProvider")
    public void testNumericRangeWithoutConstraintViolation( final String[] commandLine ) {
        // These tests succeed if no exception is thrown, since no constraints have been violated
        runNumericArgumentRangeTest(commandLine);
    }

    private void runNumericArgumentRangeTest( final String[] commandLine ) {
        parsingEngine.addArgumentSource(NumericRangeArgProvider.class);
        parsingEngine.parse(commandLine);

        NumericRangeArgProvider argProvider = new NumericRangeArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);
    }

    @SuppressWarnings("unused")
    private class VariedTypeArgProvider {
        @Argument(fullName = "intVal", required=false)
        private int anInt;

        @Argument(fullName = "stringVal", required=false)
        private String aString;

        @Argument(fullName = "enumVal", required=false)
        private TestEnum anEnum;

        @Argument(fullName = "fileVal", required=false)
        private File aFile;

        @Argument(fullName = "stringSet", required=false)
        private Set<String> someStrings;

        @Argument(fullName = "intervalVal", required=false)
        private IntervalBinding<Feature> anInterval;
    }

    @DataProvider(name = "MissingArgumentValueDataProvider")
    public Object[][] missingArgumentDataProvider() {
        return new Object[][]{
                { new String[]{"--intVal"} },
                { new String[]{"--stringVal"} },
                { new String[]{"--enumVal"} },
                { new String[]{"--fileVal"} },
                { new String[]{"--stringSet"} },
                { new String[]{"--stringSet", "aha", "--stringSet"} },
                { new String[]{"--intervalVal"} }
        };
    }

    @Test(dataProvider = "MissingArgumentValueDataProvider",
          expectedExceptions = {InvalidArgumentValueException.class, MissingArgumentValueException.class})
    public void testMissingArguments( final String[] commandLine ) {
        parsingEngine.addArgumentSource( VariedTypeArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        VariedTypeArgProvider argProvider = new VariedTypeArgProvider();
        parsingEngine.loadArgumentsIntoObject(argProvider);
    }

}

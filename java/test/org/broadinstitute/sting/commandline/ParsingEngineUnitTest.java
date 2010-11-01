/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.commandline;

import org.testng.Assert;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.List;
import java.util.EnumSet;
/**
 * Test suite for the parsing engine.
 */
public class ParsingEngineUnitTest extends BaseTest {
    private ParsingEngine parsingEngine;

    @BeforeMethod
    public void setUp() {
        parsingEngine = new ParsingEngine(null);
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

    @Test(expectedExceptions= ReviewedStingException.class)
    public void duplicateLongNameTest() {
        parsingEngine.addArgumentSource( DuplicateLongNameProvider.class );
    }

    private class DuplicateLongNameProvider {
        @Argument(fullName="myarg",doc="my arg")
        public Integer foo;

        @Argument(fullName="myarg", doc="my arg")
        public Integer bar;
    }

    @Test(expectedExceptions= ReviewedStingException.class)
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

    @Test(expectedExceptions= ReviewedStingException.class)
    public void multipleArgumentCollectionTest() {
        parsingEngine.addArgumentSource( MultipleArgumentCollectionProvider.class );
    }

    private class MultipleArgumentCollectionProvider {
        @ArgumentCollection
        RequiredArgProvider rap1 = new RequiredArgProvider();
        @ArgumentCollection
        RequiredArgProvider rap2 = new RequiredArgProvider();
    }
}

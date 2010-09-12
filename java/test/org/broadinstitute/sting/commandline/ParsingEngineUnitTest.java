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

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.junit.Test;
import org.junit.Before;
import org.junit.Assert;

import java.util.List;
import java.util.EnumSet;
/**
 * Test suite for the parsing engine.
 */
public class ParsingEngineUnitTest extends BaseTest {
    private ParsingEngine parsingEngine;

    @Before
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

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }
    
    @Test
    public void multiCharShortNameArgumentTest() {
        final String[] commandLine = new String[] {"-out","out.txt"};

        parsingEngine.addArgumentSource( MultiCharShortNameArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        MultiCharShortNameArgProvider argProvider = new MultiCharShortNameArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Argument is not correctly initialized", "out.txt", argProvider.outputFile );
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

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }

    @Test
    public void extraWhitespaceTest() {
        final String[] commandLine = new String[] {"  --input_file ", "na12878.bam"};

        parsingEngine.addArgumentSource( InputFileArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }

    @Test
    public void primitiveArgumentTest() {
        final String[] commandLine = new String[] {"--foo", "5"};

        parsingEngine.addArgumentSource( PrimitiveArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        PrimitiveArgProvider argProvider = new PrimitiveArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Argument is not correctly initialized", 5, argProvider.foo );
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

        Assert.assertTrue("Argument is not correctly initialized", argProvider.allLoci );
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

        Assert.assertEquals("Argument array is of incorrect length", 2, argProvider.inputFile.length);
        Assert.assertEquals("1st filename is incorrect", "foo.txt", argProvider.inputFile[0] );
        Assert.assertEquals("2nd filename is incorrect", "bar.txt", argProvider.inputFile[1] );        
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

        Assert.assertEquals("Enum value is not correct", TestEnum.TWO, argProvider.testEnum);
    }

    @Test
    public void enumMixedCaseTest() {
        final String[] commandLine = new String[] {  "--test_enum", "oNe" };

        parsingEngine.addArgumentSource( EnumArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        EnumArgProvider argProvider = new EnumArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Enum value is not correct", TestEnum.ONE, argProvider.testEnum);
    }
    
    @Test
    public void enumDefaultTest() {
        final String[] commandLine = new String[] {};

        parsingEngine.addArgumentSource( EnumArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        EnumArgProvider argProvider = new EnumArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Enum value is not correct", TestEnum.THREE, argProvider.testEnum);
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

        Assert.assertNotNull("Argument array is null",argProvider.integers); 
        Assert.assertEquals("Argument array is of incorrect length", 5, argProvider.integers.size());
        Assert.assertEquals("1st integer is incorrect", 2, argProvider.integers.get(0).intValue() );
        Assert.assertEquals("2nd integer is incorrect", 4, argProvider.integers.get(1).intValue() );
        Assert.assertEquals("3rd integer is incorrect", 6, argProvider.integers.get(2).intValue() );
        Assert.assertEquals("4th integer is incorrect", 8, argProvider.integers.get(3).intValue() );
        Assert.assertEquals("5th integer is incorrect",10, argProvider.integers.get(4).intValue() );
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

        Assert.assertNotNull("Argument array is null",argProvider.integers);
        Assert.assertEquals("Argument array is of incorrect length", 5, argProvider.integers.size());
        Assert.assertEquals("1st integer is incorrect", "2", argProvider.integers.get(0) );
        Assert.assertEquals("2nd integer is incorrect", "4", argProvider.integers.get(1) );
        Assert.assertEquals("3rd integer is incorrect", "6", argProvider.integers.get(2) );
        Assert.assertEquals("4th integer is incorrect", "8", argProvider.integers.get(3) );
        Assert.assertEquals("5th integer is incorrect","10", argProvider.integers.get(4) );
    }

    private class UntypedListArgProvider {
        @Argument(fullName="untyped_list",shortName="N", doc="untyped list")
        public List integers;
    }

    @Test(expected=MissingArgumentException.class)
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

        Assert.assertEquals("Default value is not correctly initialized", 42, argProvider.value.intValue() );

        // Then try to override it.
        commandLine = new String[] { "--value", "27" };

        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Default value is not correctly initialized", 27, argProvider.value.intValue() );        
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

        Assert.assertNull("Value should have remain unset",argProvider.value);
    }

    @Test
    public void unrequiredArgTest() {
        final String[] commandLine = new String[0];

        parsingEngine.addArgumentSource( UnrequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        UnrequiredArgProvider argProvider = new UnrequiredArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertNull( "Value was unrequired and unspecified; contents should be null", argProvider.value );
    }

    private class UnrequiredArgProvider {
        @Argument(required=false,doc="unrequired value")
        public Integer value;
    }

    @Test(expected=InvalidArgumentException.class)
    public void invalidArgTest() {
        final String[] commandLine = new String[] { "--foo" };

        parsingEngine.addArgumentSource( UnrequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test(expected= ReviewedStingException.class)
    public void duplicateLongNameTest() {
        parsingEngine.addArgumentSource( DuplicateLongNameProvider.class );
    }

    private class DuplicateLongNameProvider {
        @Argument(fullName="myarg",doc="my arg")
        public Integer foo;

        @Argument(fullName="myarg", doc="my arg")
        public Integer bar;
    }

    @Test(expected= ReviewedStingException.class)
    public void duplicateShortNameTest() {
        parsingEngine.addArgumentSource( DuplicateShortNameProvider.class );
    }


    private class DuplicateShortNameProvider {
        @Argument(shortName="myarg", doc="my arg")
        public Integer foo;

        @Argument(shortName="myarg", doc="my arg")
        public Integer bar;
    }

    @Test(expected=UnmatchedArgumentException.class)
    public void missingArgumentNameTest() {
        final String[] commandLine = new String[] {"foo.txt"};

        parsingEngine.addArgumentSource( NoArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    private class NoArgProvider {

    }

    @Test(expected=UnmatchedArgumentException.class)
    public void extraValueTest() {
        final String[] commandLine = new String[] {"-I", "foo.txt", "bar.txt"};

        parsingEngine.addArgumentSource( InputFileArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test(expected=MissingArgumentException.class)
    public void multipleInvalidArgTest() {
        final String[] commandLine = new String[] {"-N1", "-N2", "-N3"};

        parsingEngine.addArgumentSource( RequiredArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();
    }

    @Test(expected=TooManyValuesForArgumentException.class)
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

        Assert.assertEquals("Argument is not correctly initialized", 1, argProvider.foo.intValue() );
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

        Assert.assertEquals("Argument is not correctly initialized", 5, argProvider.bar.intValue() );
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

        Assert.assertNotNull("Invalid default argument name assigned", definition );
    }

    private class CamelCaseArgProvider {
        @Argument(doc="my arg")
        Integer myArg;
    }

    @Test(expected=UnmatchedArgumentException.class)
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

        Assert.assertEquals("Argument is not correctly initialized", "Pileup", argProvider.Analysis_Name );
    }

    private class AnalysisTypeArgProvider {
        @Argument(fullName="analysis_type", shortName="T", doc="Type of analysis to run")
        public String Analysis_Name = null;
    }

    @Test(expected=TooManyValuesForArgumentException.class)
    public void invalidParseForAnalysisTypeTest() {
        final String[] commandLine = new String[] {"--analysis_type", "Pileup", "-T", "CountReads" };

        parsingEngine.addArgumentSource( AnalysisTypeArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate( EnumSet.of(ParsingEngine.ValidationType.MissingRequiredArgument) );
    }

    @Test(expected=ArgumentsAreMutuallyExclusiveException.class)
    public void mutuallyExclusiveArgumentsTest() {
        // Passing only foo should work fine...
        String[] commandLine = new String[] {"--foo","5"};

        parsingEngine.addArgumentSource( MutuallyExclusiveArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        MutuallyExclusiveArgProvider argProvider = new MutuallyExclusiveArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Argument is not correctly initialized", 5, argProvider.foo.intValue() );

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

    @Test(expected=InvalidArgumentValueException.class)
    public void argumentValidationTest() {
        // Passing only foo should work fine...
        String[] commandLine = new String[] {"--value","521"};

        parsingEngine.addArgumentSource( ValidatingArgProvider.class );
        parsingEngine.parse( commandLine );
        parsingEngine.validate();

        ValidatingArgProvider argProvider = new ValidatingArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider );

        Assert.assertEquals("Argument is not correctly initialized", 521, argProvider.value.intValue() );

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

        Assert.assertEquals("Argument is not correctly initialized", 5, argProvider.rap.value.intValue() );
    }

    private class ArgumentCollectionProvider {
        @ArgumentCollection
        RequiredArgProvider rap = new RequiredArgProvider();
    }

    @Test(expected= ReviewedStingException.class)
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

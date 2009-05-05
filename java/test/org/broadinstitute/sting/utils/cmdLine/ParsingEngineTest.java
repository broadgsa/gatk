package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.utils.StingException;
import org.junit.Test;
import org.junit.Before;
import org.junit.Assert;

import java.util.List;
/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: May 3, 2009
 * Time: 6:05:33 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
/**
 * Test suite for the parsing engine.
 */
public class ParsingEngineTest extends BaseTest {
    private ParsingEngine parsingEngine;

    @Before
    public void setUp() {
        parsingEngine = new ParsingEngine();
    }

    private class InputFileArgProvider {
        @Argument(fullName="input_file",shortName="I")
        public String inputFile;
    }

    @Test
    public void shortNameArgumentTest() {
        final String[] commandLine = new String[] {"-I","na12878.bam"};

        parsingEngine.addArgumentSources( InputFileArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }
    
    @Test
    public void shortNameCompositeArgumentTest() {
        final String[] commandLine = new String[] {"-Ina12878.bam"};

        parsingEngine.addArgumentSources( InputFileArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }

    @Test
    public void multiCharShortNameArgumentTest() {
        final String[] commandLine = new String[] {"-out","out.txt"};

        parsingEngine.addArgumentSources( MultiCharShortNameArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        MultiCharShortNameArgProvider argProvider = new MultiCharShortNameArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", "out.txt", argProvider.outputFile );
    }


    private class MultiCharShortNameArgProvider {
        @Argument(shortName="out")
        public String outputFile;
    }

    @Test
    public void longNameArgumentTest() {
        final String[] commandLine = new String[] {"--input_file", "na12878.bam"};

        parsingEngine.addArgumentSources( InputFileArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }

    @Test
    public void extraWhitespaceTest() {
        final String[] commandLine = new String[] {"  --input_file ", "na12878.bam"};

        parsingEngine.addArgumentSources( InputFileArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam", argProvider.inputFile );
    }

    @Test
    public void flagTest() {
        final String[] commandLine = new String[] {"--all_loci"};

        parsingEngine.addArgumentSources( AllLociArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        AllLociArgProvider argProvider = new AllLociArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertTrue("Argument is not correctly initialized", argProvider.allLoci );
    }

    private class AllLociArgProvider {
        @Argument(fullName="all_loci",shortName="A")
        public boolean allLoci = false;
    }

    @Test
    public void arrayTest() {
        final String[] commandLine = new String[] {"-Ifoo.txt", "--input_file", "bar.txt"};

        parsingEngine.addArgumentSources( MultiValueArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        MultiValueArgProvider argProvider = new MultiValueArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument array is of incorrect length", 2, argProvider.inputFile.length);
        Assert.assertEquals("1st filename is incorrect", "foo.txt", argProvider.inputFile[0] );
        Assert.assertEquals("2nd filename is incorrect", "bar.txt", argProvider.inputFile[1] );        
    }

    private class MultiValueArgProvider {
        @Argument(fullName="input_file",shortName="I")
        public String[] inputFile;
    }

    @Test
    public void typedCollectionTest() {
        final String[] commandLine = new String[] { "-N2", "-N4", "-N6", "-N8", "-N10" };

        parsingEngine.addArgumentSources( IntegerListArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        IntegerListArgProvider argProvider = new IntegerListArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertNotNull("Argument array is null",argProvider.integers); 
        Assert.assertEquals("Argument array is of incorrect length", 5, argProvider.integers.size());
        Assert.assertEquals("1st integer is incorrect", 2, argProvider.integers.get(0).intValue() );
        Assert.assertEquals("2nd integer is incorrect", 4, argProvider.integers.get(1).intValue() );
        Assert.assertEquals("3rd integer is incorrect", 6, argProvider.integers.get(2).intValue() );
        Assert.assertEquals("4th integer is incorrect", 8, argProvider.integers.get(3).intValue() );
        Assert.assertEquals("5th integer is incorrect",10, argProvider.integers.get(4).intValue() );
    }

    private class IntegerListArgProvider {
        @Argument(fullName="integer_list",shortName="N")
        public List<Integer> integers;
    }

    @Test
    public void untypedCollectionTest() {
        final String[] commandLine = new String[] { "-N2", "-N4", "-N6", "-N8", "-N10" };

        parsingEngine.addArgumentSources( UntypedListArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        UntypedListArgProvider argProvider = new UntypedListArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertNotNull("Argument array is null",argProvider.integers);
        Assert.assertEquals("Argument array is of incorrect length", 5, argProvider.integers.size());
        Assert.assertEquals("1st integer is incorrect", "2", argProvider.integers.get(0) );
        Assert.assertEquals("2nd integer is incorrect", "4", argProvider.integers.get(1) );
        Assert.assertEquals("3rd integer is incorrect", "6", argProvider.integers.get(2) );
        Assert.assertEquals("4th integer is incorrect", "8", argProvider.integers.get(3) );
        Assert.assertEquals("5th integer is incorrect","10", argProvider.integers.get(4) );
    }

    private class UntypedListArgProvider {
        @Argument(fullName="untyped_list",shortName="N")
        public List integers;
    }

    @Test(expected=MissingArgumentException.class)
    public void requiredArgTest() {
        final String[] commandLine = new String[0];

        parsingEngine.addArgumentSources( RequiredArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate( argumentMatches );
    }

    private class RequiredArgProvider {
        @Argument(required=true)
        public Integer value;
    }

    @Test
    public void unrequiredArgTest() {
        final String[] commandLine = new String[0];

        parsingEngine.addArgumentSources( UnrequiredArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate( argumentMatches );

        UnrequiredArgProvider argProvider = new UnrequiredArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertNull( "Value was unrequired and unspecified; contents should be null", argProvider.value );
    }

    private class UnrequiredArgProvider {
        @Argument(required=false)
        public Integer value;
    }

    @Test(expected=InvalidArgumentException.class)
    public void invalidArgTest() {
        final String[] commandLine = new String[] { "--foo" };

        parsingEngine.addArgumentSources( UnrequiredArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate( argumentMatches );        
    }

    @Test(expected=StingException.class)
    public void duplicateLongNameTest() {
        parsingEngine.addArgumentSources( DuplicateLongNameProvider.class );        
    }

    private class DuplicateLongNameProvider {
        @Argument(fullName="myarg")
        public Integer foo;

        @Argument(fullName="myarg")
        public Integer bar;
    }

    @Test(expected=StingException.class)
    public void duplicateShortNameTest() {
        parsingEngine.addArgumentSources( DuplicateShortNameProvider.class );
    }


    private class DuplicateShortNameProvider {
        @Argument(shortName="myarg")
        public Integer foo;

        @Argument(shortName="myarg")
        public Integer bar;
    }

    @Test(expected=InvalidArgumentValueException.class)
    public void missingArgumentNameTest() {
        final String[] commandLine = new String[] {"foo.txt"};

        parsingEngine.addArgumentSources( NoArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);
    }

    private class NoArgProvider {

    }

    @Test(expected=InvalidArgumentValueException.class)
    public void extraValueTest() {
        final String[] commandLine = new String[] {"-Ifoo.txt", "bar.txt"};

        parsingEngine.addArgumentSources( InputFileArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);
    }

    @Test(expected=MissingArgumentException.class)
    public void multipleInvalidArgTest() {
        final String[] commandLine = new String[] {"-N1", "-N2", "-N3"};

        parsingEngine.addArgumentSources( RequiredArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate( argumentMatches );
    }

    @Test(expected=TooManyValuesForArgumentException.class)
    public void invalidArgCountTest() {
        final String[] commandLine = new String[] {"--value","1","--value","2","--value","3"};

        parsingEngine.addArgumentSources( RequiredArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate( argumentMatches );
    }

    @Test
    public void packageProtectedArgTest() {
        final String[] commandLine = new String[] {"--foo", "1"};

        parsingEngine.addArgumentSources( PackageProtectedArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        PackageProtectedArgProvider argProvider = new PackageProtectedArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", 1, argProvider.foo.intValue() );
    }

    private class PackageProtectedArgProvider {
        @Argument
        Integer foo;
    }

    @Test
    public void correctDefaultArgNameTest() {
        parsingEngine.addArgumentSources( CamelCaseArgProvider.class );

        DefinitionMatcher matcher = ArgumentDefinitions.FullNameDefinitionMatcher;
        ArgumentDefinition definition = parsingEngine.argumentDefinitions.findArgumentDefinition("myarg", matcher);

        Assert.assertNotNull("Invalid default argument name assigned", definition );
    }

    private class CamelCaseArgProvider {
        @Argument
        Integer myArg;
    }
}

package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.BaseTest;
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


    // To test
    // misc first element
    // multiple trailing values
    // differing input types
    // spurious arguments with in conjuction with immediate setters "-Ifoo.txt bar.txt"
    // required but missing arguments
    // invalid arguments
}

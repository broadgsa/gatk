package org.broadinstitute.sting.utils.cmdLine;

import org.broadinstitute.sting.BaseTest;
import org.junit.Test;
import org.junit.Before;
import org.junit.Assert;
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

    // To test
    // 'Composite' short names
    // long names
    // flags
    // flags with arguments at every point on the line
    // flags with arguments at the end of the line

    /*
    @Test
    public void shortNameCompositeArgumentTest() {
        final String[] commandLine = new String[] {"-I na12878.bam"};

        parsingEngine.addArgumentSources( InputFileArgProvider.class );
        ArgumentMatches argumentMatches = parsingEngine.parse( commandLine );
        parsingEngine.validate(argumentMatches);

        InputFileArgProvider argProvider = new InputFileArgProvider();
        parsingEngine.loadArgumentsIntoObject( argProvider, argumentMatches);

        Assert.assertEquals("Argument is not correctly initialized", "na12878.bam" );        
    }
    */

}

// our package
package org.broadinstitute.sting.utils.cmdLine;


// the imports for unit testing.

import org.junit.*;
import static org.junit.Assert.*;
import org.apache.commons.cli.ParseException;
import org.broadinstitute.sting.BaseTest;

/**
 * The ArgumentParserTest test class.  In JUnit 4, you don't have to extent the
 * TestCase class, instead we use annotations (@Test, @Before, @After)
 * to indicate the test pattern.
 * <p/>
 * The basic idea of a test workflow in Junit 4 is:
 * <p/>
 * 1) run the method tagged @BeforeClass
 * 2) for each method tagged with @Test {
 * run all methods tagged with @Before
 * run the @Test tagged method
 * run all methods tagged with @After
 * 3) run the method tagged with @AfterClass
 * <p/>
 * You should use the methods like
 */
public class ArgumentParserTest extends BaseTest {

    // our argument parser
    private ArgumentParser m_parser = null;

    public Boolean testBool = false;
    public String testString = "";
    public Integer testInt = 0;
    public Float testFloat = 0.0f;

    /**
     * This function (because of the @BeforeClass tag) gets called only once ever,
     * before any tests are run
     */
    @BeforeClass
    public static void doBeforeAnyTests() {

    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @AfterClass
    public static void doAfterAllTests() {

    }

    /**
     * This function does the setup of our parser, before each method call.
     * <p/>
     * Called before every test case method.
     */
    @Before
    public void doForEachTest() {
        // we don't need something done to setup each test

        // setup the parser
        m_parser = new ArgumentParser("Test Program", this);
        m_parser.addOptionalFlag("testBool", "B", "our test bool", "testBool");
        m_parser.addRequiredArg("testString", "S", "our test string", "testString");
        m_parser.addRequiredArg("testInt", "I", "our test int", "testInt");
        m_parser.addRequiredArg("testFloat", "F", "our test float", "testFloat");


    }

    /**
     * Tears down the test fixture after each call.
     * <p/>
     * Called after every test case method.
     */
    @After
    public void undoForEachTest() {
        // release objects under test here, if necessary
        m_parser = null;
    }

    /**
     * Tests that we got a string parameter in correctly
     */
    @Test
    public void testStringParameter() {
        logger.warn("Executing testStringParameter");
        // setup the parameter list
        String[] params = {"-B", "-S", "String", "-I", "100", "-F", "100.0"};

        try {
            // process the arguments
            m_parser.processArgs(params,false);
            m_parser.loadArgumentsIntoObject(this);
        } catch (ParseException e) {
            fail("We received an unexpected parsing exception");
        }

        // assert that none of the parameters are still null
        org.junit.Assert.assertNotNull(testString);
        assertEquals(testString.equals("String"), true);

    }

    /**
     * Tests that we got a Boolean parameter in correctly
     */
    @Test
    public void testBooleanParameter() {
        logger.warn("Executing testBooleanParameter");
        // setup the parameter list
        String[] params = {"-B", "-S", "String", "-I", "100", "-F", "100.0"};

        try {
            // process the arguments
            m_parser.processArgs(params,false);
            m_parser.loadArgumentsIntoObject(this);
        } catch (ParseException e) {
            fail("We received an unexpected parsing exception");
        }

        assertEquals((boolean) testBool, true);
    }

    /**
     * Tests that we got a Boolean parameter in correctly
     */
    @Test
    public void testFloatParameter() {
        logger.warn("Executing testFloatParameter");

        // setup the parameter list
        String[] params = {"-B", "-S", "String", "-I", "100", "-F", "100.0"};

        try {
            // process the arguments
            m_parser.processArgs(params,false);
            m_parser.loadArgumentsIntoObject(this);
        } catch (ParseException e) {
            fail("We received an unexpected parsing exception");
        }

        assertEquals((testFloat.compareTo(100.0f)), 0);
    }

    /**
     * Tests that we got a Integer parameter in correctly
     */
    @Test
    public void testIntegerParameter() {
        logger.warn("Executing testIntegerParameter");

        // setup the parameter list
        String[] params = {"-B", "-S", "String", "-I", "100", "-F", "100.0"};

        try {
            // process the arguments
            m_parser.processArgs(params,false);
            m_parser.loadArgumentsIntoObject(this);
        } catch (ParseException e) {
            fail("We received an unexpected parsing exception");
        }


        assertEquals(testInt.compareTo(100), 0);
    }

    /**
     * Tests that if we dont pass a required parameter we get an exception
     */
    @Test
    public void testForUnpassedParameter() {
        logger.warn("Executing testForUnpassedParameter");


        // add a new required flag we won't send it
        m_parser.addRequiredArg("testNotTHere", "N", "our should be provided test", "testFloat");

        // setup the parameter list
        String[] params = {"-B", "-S", "String", "-I", "100", "-F", "100.0"};

        try {
            // process the arguments
            m_parser.processArgs(params,false);
            fail("We should have received a missing argument exception");
        } catch (ParseException e) {
            // do nothing but consume the exception
        }
    }

    /**
     * test to see if we pass a bad field name we get an runtime exception
     */
    @Test
    public void testForBadArgFieldName() {
        logger.warn("Executing testForBadArgFieldName");

        try {
            // add a new required flag we won't send it
            m_parser.addRequiredArg("testNotTHere", "N", "our should be provided test", "testDoesNotExist");
            fail("We should of recieved a runtime exception, add unavailable fields is Baaad");
        } catch (RuntimeException e) {
            // do nothing but consume the exception
        }
    }

}


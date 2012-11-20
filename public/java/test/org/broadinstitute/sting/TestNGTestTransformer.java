package org.broadinstitute.sting;

import org.apache.log4j.Logger;
import org.testng.IAnnotationTransformer;
import org.testng.annotations.ITestAnnotation;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

/**
 * Provide default @Test values for GATK testng tests.
 *
 * Currently only sets the maximum runtime to 10 minutes, if it's not been specified.
 *
 * See http://beust.com/weblog/2006/10/18/annotation-transformers-in-java/
 *
 * @author depristo
 * @since 10/31/12
 * @version 0.1
 */
public class TestNGTestTransformer implements IAnnotationTransformer {
    public static final long DEFAULT_TIMEOUT = 1000 * 60 * 10; // 10 minutes max per test

    final static Logger logger = Logger.getLogger(TestNGTestTransformer.class);

    public void transform(ITestAnnotation annotation,
                          Class testClass,
                          Constructor testConstructor,
                          Method testMethod)
    {
        if ( annotation.getTimeOut() == 0 ) {
            logger.warn("test " + testMethod.toString() + " has no specified timeout, adding default timeout " + DEFAULT_TIMEOUT / 1000 / 60 + " minutes");
            annotation.setTimeOut(DEFAULT_TIMEOUT);
        }
    }
}


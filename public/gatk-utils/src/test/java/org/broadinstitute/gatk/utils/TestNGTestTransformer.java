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

package org.broadinstitute.gatk.utils;

import org.apache.log4j.Logger;
import org.testng.IAnnotationTransformer;
import org.testng.annotations.ITestAnnotation;

import java.lang.reflect.Constructor;
import java.lang.reflect.Method;

/**
 * Provide default @Test values for GATK testng tests.
 *
 * Currently only sets the maximum runtime to 40 minutes, if it's not been specified.
 *
 * See http://beust.com/weblog/2006/10/18/annotation-transformers-in-java/
 *
 * @author depristo
 * @since 10/31/12
 * @version 0.1
 */
public class TestNGTestTransformer implements IAnnotationTransformer {
    public static final long DEFAULT_TIMEOUT = 1000 * 60 * 40; // 40 minutes max per test

    final static Logger logger = Logger.getLogger(TestNGTestTransformer.class);

    public void transform(ITestAnnotation annotation,
                          Class testClass,
                          Constructor testConstructor,
                          Method testMethod)
    {
        if ( annotation.getTimeOut() == 0 ) {
            logger.warn("test " + (testMethod == null ? "<null>" : testMethod.toString()) + " has no specified timeout, adding default timeout " + DEFAULT_TIMEOUT / 1000 / 60 + " minutes");
            annotation.setTimeOut(DEFAULT_TIMEOUT);
        }
    }
}


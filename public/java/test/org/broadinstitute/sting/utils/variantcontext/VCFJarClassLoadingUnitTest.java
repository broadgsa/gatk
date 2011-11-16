/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.variantcontext;

import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;

/**
 * Test to ensure that, given only the VCF jar and its expected dependencies, core VCF classes will load.
 */
public class VCFJarClassLoadingUnitTest {
    @Test
    public void testVCFJarClassLoading() throws ClassNotFoundException, MalformedURLException {
        URL[] jarURLs;

        try {
            jarURLs = new URL[] { getVCFJarFile().toURI().toURL(), getTribbleJarFile().toURI().toURL() };
        }
        catch ( FileNotFoundException e ) {
            throw new ReviewedStingException("Could not find the VCF jar and/or its dependencies", e);
        }

        ClassLoader classLoader = new URLClassLoader(jarURLs, null);
        classLoader.loadClass("org.broadinstitute.sting.utils.variantcontext.VariantContext");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.VCFCodec");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.VCFWriter");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter");
    }

    /**
     * Locates the tribble jar within the dist directory.
     *
     * Makes the horrible assumption that tests will always be run from the root of a Sting clone,
     * but this is much less problematic than using the classpath to locate tribble, since
     * the classpath won't explicitly contain tribble when we're testing the fully-packaged
     * GATK jar.
     *
     * @return The tribble jar file, if found
     * @throws FileNotFoundException If we couldn't locate a tribble jar within the dist directory
     */
    private File getTribbleJarFile() throws FileNotFoundException {
        File distDir = new File("dist");
        if ( ! distDir.isDirectory() ) {
            throw new FileNotFoundException("The dist directory does not exist");
        }

        for ( File distDirEntry : distDir.listFiles() ) {
            if ( distDirEntry.getName().startsWith("tribble") && distDirEntry.getName().endsWith(".jar") ) {
                return distDirEntry;
            }
        }

        throw new FileNotFoundException("Could not find a tribble jar file in the dist directory.");
    }

    /**
     * Locates the vcf jar within the dist directory.
     *
     * Makes the horrible assumption that tests will always be run from the root of a Sting clone,
     * but this is much less problematic than using the classpath to locate vcf.jar, since
     * the classpath won't explicitly contain vcf.jar when we're testing the fully-packaged
     * GATK jar.
     *
     * @return The vcf jar file, if found
     * @throws FileNotFoundException If we couldn't locate a vcf jar within the dist directory
     */
    private File getVCFJarFile() throws FileNotFoundException {
        File vcfJar = new File("dist/vcf.jar");

        if ( ! vcfJar.exists() ) {
            throw new FileNotFoundException("Could not find dist/vcf.jar");
        }

        return vcfJar;
    }
}

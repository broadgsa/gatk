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

import org.testng.annotations.Test;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.URL;
import java.net.URLClassLoader;

/**
 * Test to ensure that, given only the VCF jar and its expected dependencies, core VCF classes will load.
 */
public class VCFJarClassLoadingUnitTest {
    @Test
    public void testVCFJarClassLoading() throws ClassNotFoundException, MalformedURLException {
        URI vcfURI = new File("dist/vcf.jar").toURI();
        URI tribbleURI = new File("lib/tribble-24.jar").toURI();

        ClassLoader classLoader = new URLClassLoader(new URL[] {vcfURI.toURL(),tribbleURI.toURL()}, null);
        classLoader.loadClass("org.broadinstitute.sting.utils.variantcontext.VariantContext");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.VCFCodec");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.VCF3Codec");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.VCFWriter");
        classLoader.loadClass("org.broadinstitute.sting.utils.codecs.vcf.StandardVCFWriter");
    }
}

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

package org.broadinstitute.gatk.utils.report;

import org.broadinstitute.gatk.utils.BaseTest;
import org.testng.annotations.Test;


/**
 * @author aaron
 *         <p/>
 *         Class ReportMarshallerUnitTest
 *         <p/>
 *         test out the marshaller
 */
public class ReportMarshallerUnitTest extends BaseTest {
    @Test
    public void testMarshalling() {
        /*Configuration cfg = new Configuration();
        try {
            cfg.setDirectoryForTemplateLoading(new File("templates"));
        } catch (IOException e) {
            e.printStackTrace(); 
        }
        cfg.setObjectWrapper(new DefaultObjectWrapper());
        Template temp = null;
        try {
            temp = cfg.createMarhsaller("myTestTemp.ftl");
        } catch (IOException e) {
            e.printStackTrace();
        }
        FakeAnalysis fa = new FakeAnalysis();
        File fl = new File("testFile.out");
        fl.deleteOnExit();
        ReportMarshaller marsh = new ReportMarshaller("report",fl,temp);
        marsh.write(fa);
        marsh.write(fa);
        marsh.write(fa);
        marsh.close();*/
    }
}

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
package org.broadinstitute.sting.gatk.executive;

/**
 * Created by IntelliJ IDEA.
 * User: mhanna
 * Date: Jan 12, 2011
 * Time: 9:19:27 PM
 * To change this template use File | Settings | File Templates.
 */
public interface MicroSchedulerMBean {
    /**
     * Gets the filename to which performance data is currently being written.
     * @return Filename to which performance data is currently being written.
     */
    public String getPerformanceLogFileName();

    /**
     * Set the filename of the log for performance.  If set,
     * @param fileName filename to use when writing performance data.
     */
    public void setPerformanceLogFileName(String fileName);

    /**
     * Gets the frequency with which performance data is written.
     * @return Frequency, in seconds, of performance log writes.
     */
    public long getPerformanceProgressPrintFrequencySeconds();    

    /**
     * How often should the performance log message be written?
     * @param seconds number of seconds between messages indicating performance frequency.
     */
    public void setPerformanceProgressPrintFrequencySeconds(long seconds);
}

/*
* The Broad Institute
* SOFTWARE COPYRIGHT NOTICE AGREEMENT
* This software and its documentation are copyright 2009 by the
* Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
*
* This software is supplied without any warranty or guaranteed support whatsoever. Neither
* the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
*/

package edu.mit.broad.picard.sam;

import edu.mit.broad.picard.metrics.MetricBase;

/**
 * Metrics class for insert size statistics
 *
 * @author Doug Voet
 */
public class InsertSizeMetrics extends MetricBase {
    public double MEDIAN_INSERT_SIZE;
	public int MIN_INSERT_SIZE;
	public int MAX_INSERT_SIZE;
	public double MEAN_INSERT_SIZE;
	public double STANDARD_DEVIATION;
    public long READ_PAIRS;

    public int WIDTH_OF_10_PERCENT;
    public int WIDTH_OF_20_PERCENT;
    public int WIDTH_OF_30_PERCENT;
    public int WIDTH_OF_40_PERCENT;
    public int WIDTH_OF_50_PERCENT;
    public int WIDTH_OF_60_PERCENT;
    public int WIDTH_OF_70_PERCENT;
    public int WIDTH_OF_80_PERCENT;
    public int WIDTH_OF_90_PERCENT;
    public int WIDTH_OF_99_PERCENT;
}

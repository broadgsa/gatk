///*
// * Copyright (c) 2012, The Broad Institute
// *
// * Permission is hereby granted, free of charge, to any person
// * obtaining a copy of this software and associated documentation
// * files (the "Software"), to deal in the Software without
// * restriction, including without limitation the rights to use,
// * copy, modify, merge, publish, distribute, sublicense, and/or sell
// * copies of the Software, and to permit persons to whom the
// * Software is furnished to do so, subject to the following
// * conditions:
// *
// * The above copyright notice and this permission notice shall be
// * included in all copies or substantial portions of the Software.
// * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// * OTHER DEALINGS IN THE SOFTWARE.
// */
//
//package org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.genotypePhasingEvaluator;
//
//import org.broadinstitute.sting.gatk.walkers.varianteval.util.TableType;
//
//import java.util.HashMap;
//import java.util.Map;
//
///**
// * a table of sample names to genotype phasing statistics
// */
//class SamplePhasingStatistics extends TableType {
//    private HashMap<String, PhaseStats> sampleStats = null;
//    private double minPhaseQuality;
//
//    public SamplePhasingStatistics(double minPhaseQuality) {
//        this.sampleStats = new HashMap<String, PhaseStats>();
//        this.minPhaseQuality = minPhaseQuality;
//    }
//
//    public PhaseStats ensureSampleStats(String samp) {
//        PhaseStats ps = sampleStats.get(samp);
//        if (ps == null) {
//            ps = new PhaseStats();
//            sampleStats.put(samp, ps);
//        }
//        return ps;
//    }
//
//    /**
//     * @return one row per sample
//     */
//    public String[] getRowKeys() {
//        return sampleStats.keySet().toArray(new String[sampleStats.size()]);
//    }
//
//    /**
//     * get the column keys
//     *
//     * @return a list of objects, in this case strings, that are the column names
//     */
//    public String[] getColumnKeys() {
//        return PhaseStats.getFieldNamesArray();
//    }
//
//    public Object getCell(int x, int y) {
//        String[] rowKeys = getRowKeys();
//        PhaseStats ps = sampleStats.get(rowKeys[x]);
//        return ps.getField(y);
//    }
//
//    public String getName() {
//        return "Sample Phasing Statistics (for PQ >= " + minPhaseQuality + ")";
//    }
//
//    public String toString() {
//        StringBuilder sb = new StringBuilder();
//        for (Map.Entry<String, PhaseStats> sampPhaseStatsEnt : sampleStats.entrySet()) {
//            String sample = sampPhaseStatsEnt.getKey();
//            PhaseStats ps = sampPhaseStatsEnt.getValue();
//
//            sb.append(sample + "\t" + ps);
//        }
//        return sb.toString();
//    }
//}

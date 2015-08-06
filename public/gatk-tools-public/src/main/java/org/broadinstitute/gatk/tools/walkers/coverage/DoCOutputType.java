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

package org.broadinstitute.gatk.tools.walkers.coverage;

/**
 * Models a single output file in the DoC walker.
 *
 * @author mhanna
 * @version 0.1
 */
public class DoCOutputType {
    public enum Partition { readgroup, sample, library, platform, center, sample_by_platform, sample_by_center, sample_by_platform_by_center }
    public enum Aggregation { locus, interval, gene, cumulative }
    public enum FileType { summary, statistics, coverage_counts, coverage_proportions }

    private final Partition partition;
    private final Aggregation aggregation;
    private final FileType fileType;

    public DoCOutputType(final Partition partition,
                         final Aggregation aggregation,
                         final FileType fileType) {
        this.partition = partition;
        this.aggregation = aggregation;
        this.fileType = fileType;
    }

    public String getFileName(final String baseName) {
        // main output
        if(partition == null)
            return baseName;

        if(baseName.trim().equals("/dev/null"))
            return "/dev/null";

        // mhanna 22 Aug 2010 - Deliberately force this header replacement to make sure integration tests pass.
        // TODO: Update integration tests and get rid of this.
        String partitionType = (partition == DoCOutputType.Partition.readgroup ? "read_group" : partition.toString());        

        if(fileType == FileType.coverage_counts || fileType == FileType.coverage_proportions) {
            // coverage counts / proportions files always include aggregation.
            return baseName + "." +
                    partitionType + "_" +
                    aggregation + "_" +
                    fileType;
        }

        return  baseName + "." +
                partitionType + "_" +
                (aggregation == Aggregation.interval || aggregation == Aggregation.gene ? aggregation + "_" : "") +
                fileType;
    }

    public int hashCode() {
        return (partition!=null?partition.ordinal()+1:0) * aggregation.ordinal() * fileType.ordinal();
    }

    public boolean equals(Object other) {
        if(!(other instanceof DoCOutputType))
            return false;
        DoCOutputType otherOutputType = (DoCOutputType)other;
        return partition == otherOutputType.partition &&
                aggregation == otherOutputType.aggregation &&
                fileType == otherOutputType.fileType;
    }
}

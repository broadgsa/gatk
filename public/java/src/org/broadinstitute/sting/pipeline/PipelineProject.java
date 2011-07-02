/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.pipeline;

import java.io.File;
import java.util.Map;
import java.util.TreeMap;

/**
 * Java bean defining the project for a pipeline.
 */
public class PipelineProject {
    private String name;
    private File referenceFile;
    private File intervalList;
    private File genotypeDbsnp;
    private File evalDbsnp;
    private File refseqTable;
    private Map<String, String> tags = new TreeMap<String, String>();

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public File getIntervalList() {
        return intervalList;
    }

    public void setIntervalList(File intervalList) {
        this.intervalList = intervalList;
    }

    public File getReferenceFile() {
        return referenceFile;
    }

    public void setReferenceFile(File referenceFile) {
        this.referenceFile = referenceFile;
    }

    public File getGenotypeDbsnp() {
        return genotypeDbsnp;
    }

    public void setGenotypeDbsnp(File genotypeDbsnp) {
        this.genotypeDbsnp = genotypeDbsnp;
    }

    public String getGenotypeDbsnpType() {
        return getDbsnpType(genotypeDbsnp);
    }

    public File getEvalDbsnp() {
        return evalDbsnp;
    }

    public void setEvalDbsnp(File evalDbsnp) {
        this.evalDbsnp = evalDbsnp;
    }

    public String getEvalDbsnpType() {
        return getDbsnpType(evalDbsnp);
    }

    public File getRefseqTable() {
        return refseqTable;
    }

    public void setRefseqTable(File refseqTable) {
        this.refseqTable = refseqTable;
    }

    public Map<String, String> getTags() {
        return tags;
    }

    public void setTags(Map<String, String> tags) {
        this.tags = tags;
    }

    private String getDbsnpType(File file) {
        if (file == null)
            return null;
        else if (file.getName().toLowerCase().endsWith(".vcf"))
            return "vcf";
        else
            return "dbsnp";
    }
}

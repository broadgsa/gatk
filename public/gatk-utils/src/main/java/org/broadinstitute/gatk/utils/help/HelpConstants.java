/*
* Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.gatk.utils.help;

public class HelpConstants {

    public final static String BASE_GATK_URL = "https://www.broadinstitute.org/gatk";
    public final static String GATK_DOCS_URL = BASE_GATK_URL + "/guide/tooldocs/";
    public final static String GATK_ARTICLE_URL = BASE_GATK_URL + "/guide/article";
    public final static String GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";
    public final static String GATK_FORUM_API_URL = "https://gatkforums.broadinstitute.org/api/v1/";

    /**
     * Arguments for parallelism options
     */
    public final static String ARG_TREEREDUCIBLE = "-nt";
    public final static String ARG_NANOSCHEDULABLE = "-nct";

    /**
     * Definition of the group names / categories of tools.
     * The names get parsed to make supercategories in the doc index,
     * so be careful when making big changes -- see GATKDoclet.java toMap()
     */
    public final static String DOCS_CAT_ENGINE = "Engine Parameters (available to all tools)";
    public final static String DOCS_CAT_QC = "Diagnostics and Quality Control Tools";
    public final static String DOCS_CAT_DATA = "Sequence Data Processing Tools";
    public final static String DOCS_CAT_VARDISC = "Variant Discovery Tools";
    public final static String DOCS_CAT_VAREVAL = "Variant Evaluation Tools";
    public final static String DOCS_CAT_VARMANIP = "Variant Manipulation Tools";
    public final static String DOCS_CAT_ANNOT = "Annotation Modules";
    public final static String DOCS_CAT_RF = "Read Filters";
    public final static String DOCS_CAT_RODCODECS = "Resource File Codecs";
    public final static String DOCS_CAT_REFUTILS = "Reference Utilities";
    public final static String DOCS_CAT_USRERR = "User Exceptions (Exclude)";
    public final static String DOCS_CAT_TOY = "Toy Examples (Exclude)";

    public static String articlePost(Integer id) {
        return GATK_ARTICLE_URL + "?id=" + id.toString();
    }

}
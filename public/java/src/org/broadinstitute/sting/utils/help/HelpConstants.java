/*
* Copyright (c) 2012 The Broad Institute
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

package org.broadinstitute.sting.utils.help;

public class HelpConstants {

    public final static String BASE_GATK_URL = "http://www.broadinstitute.org/gatk";
    public final static String GATK_DOCS_URL = BASE_GATK_URL + "/gatkdocs/";
    public final static String GATK_FORUM_URL = "http://gatkforums.broadinstitute.org/";
    public final static String GATK_FORUM_API_URL = "https://gatkforums.broadinstitute.org/api/v1/";

    /**
     * Arguments for parallelism options
     */
    public final static String ARG_TREEREDUCIBLE = "-nt";
    public final static String ARG_NANOSCHEDULABLE = "-nct";
    public final static String CMDLINE_GATK_URL = GATK_DOCS_URL + "org_broadinstitute_sting_gatk_CommandLineGATK.html";

    /**
     * Definition of the group names / categories of tools.
     * The names get parsed to make supercategories in the doc index,
     * so be careful when making big changes -- see GATKDoclet.java toMap()
     */
    public final static String DOCS_CAT_DATA = "Sequence Data Processing Tools";
    public final static String DOCS_CAT_QC = "Diagnostics and Quality Control Tools";
    public final static String DOCS_CAT_ENGINE = "Engine Parameters (available to all tools)";
    public final static String DOCS_CAT_RF = "Read Filters";
    public final static String DOCS_CAT_REFUTILS = "Reference Utilities";
    public final static String DOCS_CAT_RODCODECS = "ROD Codecs";
    public final static String DOCS_CAT_USRERR = "User Exceptions";
    public final static String DOCS_CAT_VALIDATION = "Validation Utilities";
    public final static String DOCS_CAT_ANNOT = "Variant Annotations";
    public final static String DOCS_CAT_VARDISC = "Variant Discovery Tools";
    public final static String DOCS_CAT_VARMANIP = "Variant Evaluation and Manipulation Tools";
    public final static String DOCS_CAT_TEST = "Testing Tools";
    public final static String DOCS_CAT_HELPUTILS = "Help Utilities";

    public static String forumPost(String post) {
    	return GATK_FORUM_URL + post;
    }
}
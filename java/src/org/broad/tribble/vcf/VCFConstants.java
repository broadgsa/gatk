/*
 * Copyright (c) 2010.
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

package org.broad.tribble.vcf;

public final class VCFConstants {

    // standard INFO/FORMAT field keys
    public static final String ANCESTRAL_ALLELE_KEY = "AA";
    public static final String ALLELE_COUNT_KEY = "AC";
    public static final String ALLELE_FREQUENCY_KEY = "AF";
    public static final String ALLELE_NUMBER_KEY = "AN";
    public static final String RMS_BASE_QUALITY_KEY = "BQ";
    public static final String CIGAR_KEY = "CIGAR";
    public static final String DBSNP_KEY = "DB";
    public static final String DEPTH_KEY = "DP";
    public static final String END_KEY = "END";
    public static final String GENOTYPE_FILTER_KEY = "FT";
    public static final String GENOTYPE_KEY = "GT";
    public static final String GENOTYPE_LIKELIHOODS_KEY = "GL";
    public static final String GENOTYPE_QUALITY_KEY = "GQ";
    public static final String HAPMAP2_KEY = "H2";
    public static final String HAPMAP3_KEY = "H3";
    public static final String HAPLOTYPE_QUALITY_KEY = "HQ";
    public static final String RMS_MAPPING_QUALITY_KEY = "MQ";
    public static final String MAPPING_QUALITY_ZERO_KEY = "MQ0";
    public static final String SAMPLE_NUMBER_KEY = "NS";
    public static final String OLD_DEPTH_KEY = "RD";
    public static final String STRAND_BIAS_KEY = "SB";
    public static final String SOMATIC_KEY = "SOMATIC";
    public static final String VALIDATED_KEY = "VALIDATED";

    // separators
    public static final String FORMAT_FIELD_SEPARATOR = ":";
    public static final String GENOTYPE_FIELD_SEPARATOR = ":";
    public static final String FIELD_SEPARATOR = "\t";
    public static final String FILTER_CODE_SEPARATOR = ";";
    public static final String INFO_FIELD_SEPARATOR = ";";

    // missing/default values
    public static final String UNFILTERED = ".";
    public static final String PASSES_FILTERS_v3 = "0";
    public static final String PASSES_FILTERS_v4 = "PASS";
    public static final String EMPTY_ID_FIELD = ".";
    public static final String EMPTY_INFO_FIELD = ".";
    public static final String EMPTY_ALTERNATE_ALLELE_FIELD = ".";
    public static final String MISSING_VALUE_v4 = ".";
    public static final String MISSING_QUALITY_v3 = "-1";
    public static final String MISSING_GENOTYPE_QUALITY_v3 = "-1";
    public static final String MISSING_HAPLOTYPE_QUALITY_v3 = "-1";
    public static final String MISSING_DEPTH_v3 = "-1";
    public static final String UNBOUNDED_ENCODING_v4 = ".";
    public static final String UNBOUNDED_ENCODING_v3 = "-1";
    public static final String EMPTY_ALLELE = ".";
    public static final String EMPTY_GENOTYPE = "./.";
    public static final double MAX_GENOTYPE_QUAL = 99.0;

    public static final String DOUBLE_PRECISION_FORMAT_STRING = "%.2f";
}
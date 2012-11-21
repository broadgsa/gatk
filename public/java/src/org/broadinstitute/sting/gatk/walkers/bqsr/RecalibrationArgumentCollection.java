/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.bqsr;

import org.broad.tribble.Feature;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.utils.recalibration.RecalUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.Collections;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 27, 2009
 *
 * A collection of the arguments that are common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 */

public class RecalibrationArgumentCollection {

    /**
     * This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference,
     * so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of RodBindings (VCF, Bed, etc.)
     * for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites.
     * Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
     */
    @Input(fullName = "knownSites", shortName = "knownSites", doc = "A database of known polymorphic sites to skip over in the recalibration algorithm", required = false)
    public List<RodBinding<Feature>> knownSites = Collections.emptyList();

    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.   Use '/dev/stdout' to print to standard out.
     */
    @Gather(BQSRGatherer.class)
    @Output(doc = "The output recalibration table file to create", required = true)
    public File RECAL_TABLE_FILE = null;
    public PrintStream RECAL_TABLE;

    /**
     * If not provided, then no plots will be generated (useful for queue scatter/gathering).
     * However, we *highly* recommend that users generate these plots whenever possible for QC checking.
     */
    @Output(fullName = "plot_pdf_file", shortName = "plots", doc = "The output recalibration pdf file to create", required = false)
    public File RECAL_PDF_FILE = null;

    /**
     * If not provided, then a temporary file is created and then deleted upon completion.
     * For advanced users only.
     */
    @Advanced
    @Argument(fullName = "intermediate_csv_file", shortName = "intermediate", doc = "The intermediate csv file to create", required = false)
    public File RECAL_CSV_FILE = null;

    /**
     * Note that the --list argument requires a fully resolved and correct command-line to work.
     */
    @Argument(fullName = "list", shortName = "ls", doc = "List the available covariates and exit", required = false)
    public boolean LIST_ONLY = false;

    /**
     * Note that the ReadGroup and QualityScore covariates are required and do not need to be specified.
     * Also, unless --no_standard_covs is specified, the Cycle and Context covariates are standard and are included by default.
     * Use the --list argument to see the available covariates.
     */
    @Argument(fullName = "covariate", shortName = "cov", doc = "One or more covariates to be used in the recalibration. Can be specified multiple times", required = false)
    public String[] COVARIATES = null;

    /*
     * The Cycle and Context covariates are standard and are included by default unless this argument is provided.
     * Note that the ReadGroup and QualityScore covariates are required and cannot be excluded.
     */
    @Argument(fullName = "no_standard_covs", shortName = "noStandard", doc = "Do not use the standard set of covariates, but rather just the ones listed using the -cov argument", required = false)
    public boolean DO_NOT_USE_STANDARD_COVARIATES = false;

    /**
     * This calculation is critically dependent on being able to skip over known polymorphic sites. Please be sure that you know what you are doing if you use this option.
     */
    @Advanced
    @Argument(fullName = "run_without_dbsnp_potentially_ruining_quality", shortName = "run_without_dbsnp_potentially_ruining_quality", required = false, doc = "If specified, allows the recalibrator to be used without a dbsnp rod. Very unsafe and for expert users only.")
    public boolean RUN_WITHOUT_DBSNP = false;

    /**
     * CountCovariates and TableRecalibration accept a --solid_recal_mode <MODE> flag which governs how the recalibrator handles the
     * reads which have had the reference inserted because of color space inconsistencies.
     */
    @Argument(fullName = "solid_recal_mode", shortName = "sMode", required = false, doc = "How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS")
    public RecalUtils.SOLID_RECAL_MODE SOLID_RECAL_MODE = RecalUtils.SOLID_RECAL_MODE.SET_Q_ZERO;

    /**
     * CountCovariates and TableRecalibration accept a --solid_nocall_strategy <MODE> flag which governs how the recalibrator handles
     * no calls in the color space tag. Unfortunately because of the reference inserted bases mentioned above, reads with no calls in
     * their color space tag can not be recalibrated.
     */
    @Argument(fullName = "solid_nocall_strategy", shortName = "solid_nocall_strategy", doc = "Defines the behavior of the recalibrator when it encounters no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ", required = false)
    public RecalUtils.SOLID_NOCALL_STRATEGY SOLID_NOCALL_STRATEGY = RecalUtils.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION;

    /**
     * The context covariate will use a context of this size to calculate it's covariate value for base mismatches
     */
    @Argument(fullName = "mismatches_context_size", shortName = "mcs", doc = "size of the k-mer context to be used for base mismatches", required = false)
    public int MISMATCHES_CONTEXT_SIZE = 2;

    /**
     * The context covariate will use a context of this size to calculate it's covariate value for base insertions and deletions
     */
    @Argument(fullName = "indels_context_size", shortName = "ics", doc = "size of the k-mer context to be used for base insertions and deletions", required = false)
    public int INDELS_CONTEXT_SIZE = 3;

    /**
     * The cycle covariate will generate an error if it encounters a cycle greater than this value.
     * This argument is ignored if the Cycle covariate is not used.
     */
    @Argument(fullName = "maximum_cycle_value", shortName = "maxCycle", doc = "the maximum cycle value permitted for the Cycle covariate", required = false)
    public int MAXIMUM_CYCLE_VALUE = 500;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off (default is off)
     */
    @Argument(fullName = "mismatches_default_quality", shortName = "mdq", doc = "default quality for the base mismatches covariate", required = false)
    public byte MISMATCHES_DEFAULT_QUALITY = -1;

    /**
     * A default base qualities to use as a prior (reported quality) in the insertion covariate model. This parameter is used for all reads without insertion quality scores for each base. (default is on)
     */
    @Argument(fullName = "insertions_default_quality", shortName = "idq", doc = "default quality for the base insertions covariate", required = false)
    public byte INSERTIONS_DEFAULT_QUALITY = 45;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off (default is off)
     */
    @Argument(fullName = "deletions_default_quality", shortName = "ddq", doc = "default quality for the base deletions covariate", required = false)
    public byte DELETIONS_DEFAULT_QUALITY = 45;

    /**
     * Reads with low quality bases on either tail (beginning or end) will not be considered in the context. This parameter defines the quality below which (inclusive) a tail is considered low quality
     */
    @Argument(fullName = "low_quality_tail", shortName = "lqt", doc = "minimum quality for the bases in the tail of the reads to be considered", required = false)
    public byte LOW_QUAL_TAIL = 2;

    /**
     * BQSR generates a quantization table for quick quantization later by subsequent tools. BQSR does not quantize the base qualities, this is done by the engine with the -qq or -BQSR options.
     * This parameter tells BQSR the number of levels of quantization to use to build the quantization table.
     */
    @Argument(fullName = "quantizing_levels", shortName = "ql", required = false, doc = "number of distinct quality scores in the quantized output")
    public int QUANTIZING_LEVELS = 16;

    /**
     * The tag name for the binary tag covariate (if using it)
     */
    @Argument(fullName = "binary_tag_name", shortName = "bintag", required = false, doc = "the binary tag covariate name if using it")
    public String BINARY_TAG_NAME = null;


    /////////////////////////////
    // Debugging-only Arguments
    /////////////////////////////

    @Hidden
    @Argument(fullName = "default_platform", shortName = "dP", required = false, doc = "If a read has no platform then default to the provided String. Valid options are illumina, 454, and solid.")
    public String DEFAULT_PLATFORM = null;

    @Hidden
    @Argument(fullName = "force_platform", shortName = "fP", required = false, doc = "If provided, the platform of EVERY read will be forced to be the provided String. Valid options are illumina, 454, and solid.")
    public String FORCE_PLATFORM = null;

    @Hidden
    @Output(fullName = "recal_table_update_log", shortName = "recal_table_update_log", required = false, doc = "If provided, log all updates to the recalibration tables to the given file. For debugging/testing purposes only")
    public PrintStream RECAL_TABLE_UPDATE_LOG = null;

    public File existingRecalibrationReport = null;

    public GATKReportTable generateReportTable(final String covariateNames) {
        GATKReportTable argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run", 2);
        argumentsTable.addColumn("Argument");
        argumentsTable.addColumn(RecalUtils.ARGUMENT_VALUE_COLUMN_NAME);
        argumentsTable.addRowID("covariate", true);
        argumentsTable.set("covariate", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, covariateNames);
        argumentsTable.addRowID("no_standard_covs", true);
        argumentsTable.set("no_standard_covs", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DO_NOT_USE_STANDARD_COVARIATES);
        argumentsTable.addRowID("run_without_dbsnp", true);
        argumentsTable.set("run_without_dbsnp", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, RUN_WITHOUT_DBSNP);
        argumentsTable.addRowID("solid_recal_mode", true);
        argumentsTable.set("solid_recal_mode", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, SOLID_RECAL_MODE);
        argumentsTable.addRowID("solid_nocall_strategy", true);
        argumentsTable.set("solid_nocall_strategy", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, SOLID_NOCALL_STRATEGY);
        argumentsTable.addRowID("mismatches_context_size", true);
        argumentsTable.set("mismatches_context_size", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MISMATCHES_CONTEXT_SIZE);
        argumentsTable.addRowID("indels_context_size", true);
        argumentsTable.set("indels_context_size", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INDELS_CONTEXT_SIZE);
        argumentsTable.addRowID("mismatches_default_quality", true);
        argumentsTable.set("mismatches_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MISMATCHES_DEFAULT_QUALITY);
        argumentsTable.addRowID("insertions_default_quality", true);
        argumentsTable.set("insertions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INSERTIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("low_quality_tail", true);
        argumentsTable.set("low_quality_tail", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, LOW_QUAL_TAIL);
        argumentsTable.addRowID("default_platform", true);
        argumentsTable.set("default_platform", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DEFAULT_PLATFORM);
        argumentsTable.addRowID("force_platform", true);
        argumentsTable.set("force_platform", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, FORCE_PLATFORM);
        argumentsTable.addRowID("quantizing_levels", true);
        argumentsTable.set("quantizing_levels", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, QUANTIZING_LEVELS);
        argumentsTable.addRowID("recalibration_report", true);
        argumentsTable.set("recalibration_report", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, existingRecalibrationReport == null ? "null" : existingRecalibrationReport.getAbsolutePath());
        argumentsTable.addRowID("plot_pdf_file", true);
        argumentsTable.set("plot_pdf_file", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, RECAL_PDF_FILE == null ? "null" : RECAL_PDF_FILE.getAbsolutePath());
        argumentsTable.addRowID("binary_tag_name", true);
        argumentsTable.set("binary_tag_name", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, BINARY_TAG_NAME == null ? "null" : BINARY_TAG_NAME);
        return argumentsTable;
    }

}

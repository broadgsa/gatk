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

package org.broadinstitute.gatk.engine.recalibration;

import com.google.java.contract.Requires;
import htsjdk.tribble.Feature;
import org.broadinstitute.gatk.utils.commandline.*;
import org.broadinstitute.gatk.utils.report.GATKReportTable;
import org.broadinstitute.gatk.utils.Utils;
import org.broadinstitute.gatk.utils.exceptions.GATKException;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: rpoplin
 * Date: Nov 27, 2009
 *
 * A collection of the arguments that are used for BQSR. Used to be common to both CovariateCounterWalker and TableRecalibrationWalker.
 * This set of arguments will also be passed to the constructor of every Covariate when it is instantiated.
 */

public class RecalibrationArgumentCollection implements Cloneable {

    /**
     * This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference,
     * so it is critical that a database of known polymorphic sites (e.g. dbSNP) is given to the tool in order to mask out those sites.
     */
    @Input(fullName = "knownSites", shortName = "knownSites", doc = "A database of known polymorphic sites", required = false)
    public List<RodBinding<Feature>> knownSites = Collections.emptyList();

    /**
     * After the header, data records occur one per line until the end of the file. The first several items on a line are the
     * values of the individual covariates and will change depending on which covariates were specified at runtime. The last
     * three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches,
     * and the raw empirical quality score calculated by phred-scaling the mismatch rate.
     */
    @Gather(BQSRGatherer.class)
    @Output(doc = "The output recalibration table file to create", required = true)
    public File RECAL_TABLE_FILE = null;
    public PrintStream RECAL_TABLE;

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

    /**
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
     * BaseRecalibrator accepts a --solid_recal_mode <MODE> flag which governs how the recalibrator handles the
     * reads which have had the reference inserted because of color space inconsistencies.
     */
    @Argument(fullName = "solid_recal_mode", shortName = "sMode", required = false, doc = "How should we recalibrate solid bases in which the reference was inserted? Options = DO_NOTHING, SET_Q_ZERO, SET_Q_ZERO_BASE_N, or REMOVE_REF_BIAS")
    public RecalUtils.SOLID_RECAL_MODE SOLID_RECAL_MODE = RecalUtils.SOLID_RECAL_MODE.SET_Q_ZERO;

    /**
     * BaseRecalibrator accepts a --solid_nocall_strategy <MODE> flag which governs how the recalibrator handles
     * no calls in the color space tag. Unfortunately because of the reference inserted bases mentioned above, reads with no calls in
     * their color space tag can not be recalibrated.
     */
    @Argument(fullName = "solid_nocall_strategy", shortName = "solid_nocall_strategy", doc = "Defines the behavior of the recalibrator when it encounters no calls in the color space. Options = THROW_EXCEPTION, LEAVE_READ_UNRECALIBRATED, or PURGE_READ", required = false)
    public RecalUtils.SOLID_NOCALL_STRATEGY SOLID_NOCALL_STRATEGY = RecalUtils.SOLID_NOCALL_STRATEGY.THROW_EXCEPTION;

    /**
     * The context covariate will use a context of this size to calculate its covariate value for base mismatches. Must be between 1 and 13 (inclusive). Note that higher values will increase runtime and required java heap size.
     */
    @Argument(fullName = "mismatches_context_size", shortName = "mcs", doc = "Size of the k-mer context to be used for base mismatches", required = false)
    public int MISMATCHES_CONTEXT_SIZE = 2;

    /**
     * The context covariate will use a context of this size to calculate its covariate value for base insertions and deletions. Must be between 1 and 13 (inclusive). Note that higher values will increase runtime and required java heap size.
     */
    @Argument(fullName = "indels_context_size", shortName = "ics", doc = "Size of the k-mer context to be used for base insertions and deletions", required = false)
    public int INDELS_CONTEXT_SIZE = 3;

    /**
     * The cycle covariate will generate an error if it encounters a cycle greater than this value.
     * This argument is ignored if the Cycle covariate is not used.
     */
    @Argument(fullName = "maximum_cycle_value", shortName = "maxCycle", doc = "The maximum cycle value permitted for the Cycle covariate", required = false)
    public int MAXIMUM_CYCLE_VALUE = 500;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off. [default is off]
     */
    @Advanced
    @Argument(fullName = "mismatches_default_quality", shortName = "mdq", doc = "default quality for the base mismatches covariate", required = false)
    public byte MISMATCHES_DEFAULT_QUALITY = -1;

    /**
     * A default base qualities to use as a prior (reported quality) in the insertion covariate model. This parameter is used for all reads without insertion quality scores for each base. [default is on]
     */
    @Advanced
    @Argument(fullName = "insertions_default_quality", shortName = "idq", doc = "default quality for the base insertions covariate", required = false)
    public byte INSERTIONS_DEFAULT_QUALITY = 45;

    /**
     * A default base qualities to use as a prior (reported quality) in the mismatch covariate model. This value will replace all base qualities in the read for this default value. Negative value turns it off. [default is on]
     */
    @Advanced
    @Argument(fullName = "deletions_default_quality", shortName = "ddq", doc = "default quality for the base deletions covariate", required = false)
    public byte DELETIONS_DEFAULT_QUALITY = 45;

    /**
     * Reads with low quality bases on either tail (beginning or end) will not be considered in the context. This parameter defines the quality below which (inclusive) a tail is considered low quality
     */
    @Advanced
    @Argument(fullName = "low_quality_tail", shortName = "lqt", doc = "minimum quality for the bases in the tail of the reads to be considered", required = false)
    public byte LOW_QUAL_TAIL = 2;

    /**
     * BQSR generates a quantization table for quick quantization later by subsequent tools. BQSR does not quantize the base qualities, this is done by the engine with the -qq or -BQSR options.
     * This parameter tells BQSR the number of levels of quantization to use to build the quantization table.
     */
    @Advanced
    @Argument(fullName = "quantizing_levels", shortName = "ql", required = false, doc = "number of distinct quality scores in the quantized output")
    public int QUANTIZING_LEVELS = 16;

    /**
     * The tag name for the binary tag covariate (if using it)
     */
    @Advanced
    @Argument(fullName = "binary_tag_name", shortName = "bintag", required = false, doc = "the binary tag covariate name if using it")
    public String BINARY_TAG_NAME = null;

    /**
     * Whether GATK report tables should have rows in sorted order, starting from leftmost column
     */
    @Argument(fullName = "sort_by_all_columns", shortName = "sortAllCols", doc = "Sort the rows in the tables of reports", required = false)
    public Boolean SORT_BY_ALL_COLUMNS  = false;

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
    @Argument(fullName = "force_readgroup", shortName = "fRG", required = false, doc = "If provided, the read group of EVERY read will be forced to be the provided String.")
    public String FORCE_READGROUP = null;

    @Hidden
    @Output(fullName = "recal_table_update_log", shortName = "recal_table_update_log", required = false, doc = "If provided, log all updates to the recalibration tables to the given file. For debugging/testing purposes only", defaultToStdout = false)
    public PrintStream RECAL_TABLE_UPDATE_LOG = null;

    /**
     * The repeat covariate will use a context of this size to calculate its covariate value for base insertions and deletions
     */
    @Hidden
    @Argument(fullName = "max_str_unit_length", shortName = "maxstr", doc = "Max size of the k-mer context to be used for repeat covariates", required = false)
    public int MAX_STR_UNIT_LENGTH = 8;

    @Hidden
    @Argument(fullName = "max_repeat_length", shortName = "maxrep", doc = "Max number of repetitions to be used for repeat covariates", required = false)
    public int MAX_REPEAT_LENGTH = 20;


    public File existingRecalibrationReport = null;

    public GATKReportTable generateReportTable(final String covariateNames) {
        GATKReportTable argumentsTable;
        if(SORT_BY_ALL_COLUMNS) {
            argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run", 2, GATKReportTable.TableSortingWay.SORT_BY_COLUMN);
        } else {
            argumentsTable = new GATKReportTable("Arguments", "Recalibration argument collection values used in this run", 2);
        }
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
        argumentsTable.addRowID("deletions_default_quality", true);
        argumentsTable.set("deletions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, DELETIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("insertions_default_quality", true);
        argumentsTable.set("insertions_default_quality", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, INSERTIONS_DEFAULT_QUALITY);
        argumentsTable.addRowID("maximum_cycle_value", true);
        argumentsTable.set("maximum_cycle_value", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, MAXIMUM_CYCLE_VALUE);
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
        argumentsTable.addRowID("binary_tag_name", true);
        argumentsTable.set("binary_tag_name", RecalUtils.ARGUMENT_VALUE_COLUMN_NAME, BINARY_TAG_NAME == null ? "null" : BINARY_TAG_NAME);
        return argumentsTable;
    }

    /**
     * Returns a map with the arguments that differ between this an
     * another {@link RecalibrationArgumentCollection} instance.
     * <p/>
     * The key is the name of that argument in the report file. The value is a message
     * that explains the difference to the end user.
     * <p/>
     * Thus, a empty map indicates that there is no differences between both argument collection that
     * is relevant to report comparison.
     * <p/>
     * This method should not throw any exception.
     *
     * @param other the argument-collection to compare against.
     * @param thisRole the name used to refer to this RAC report that makes sense to the end user.
     * @param otherRole the name used to refer to the other RAC report that makes sense to the end user.
     *
     * @return never <code>null</code>, but a zero-size collection if there are no differences.
     */
    @Requires("other != null && thisRole != null && otherRole != null && !thisRole.equalsIgnoreCase(otherRole)")
    public Map<String,? extends CharSequence> compareReportArguments(final RecalibrationArgumentCollection other,final String thisRole, final String otherRole) {
        final Map<String,String> result = new LinkedHashMap<>(15);
        compareRequestedCovariates(result, other, thisRole, otherRole);
        compareSimpleReportArgument(result,"no_standard_covs", DO_NOT_USE_STANDARD_COVARIATES, other.DO_NOT_USE_STANDARD_COVARIATES, thisRole, otherRole);
        compareSimpleReportArgument(result,"run_without_dbsnp",RUN_WITHOUT_DBSNP,other.RUN_WITHOUT_DBSNP,thisRole,otherRole);
        compareSimpleReportArgument(result,"solid_recal_mode", SOLID_RECAL_MODE, other.SOLID_RECAL_MODE,thisRole,otherRole);
        compareSimpleReportArgument(result,"solid_nocall_strategy", SOLID_NOCALL_STRATEGY, other.SOLID_NOCALL_STRATEGY,thisRole,otherRole);
        compareSimpleReportArgument(result,"mismatches_context_size", MISMATCHES_CONTEXT_SIZE,other.MISMATCHES_CONTEXT_SIZE,thisRole,otherRole);
        compareSimpleReportArgument(result,"mismatches_default_quality", MISMATCHES_DEFAULT_QUALITY, other.MISMATCHES_DEFAULT_QUALITY,thisRole,otherRole);
        compareSimpleReportArgument(result,"deletions_default_quality", DELETIONS_DEFAULT_QUALITY, other.DELETIONS_DEFAULT_QUALITY,thisRole,otherRole);
        compareSimpleReportArgument(result,"insertions_default_quality", INSERTIONS_DEFAULT_QUALITY, other.INSERTIONS_DEFAULT_QUALITY,thisRole,otherRole);
        compareSimpleReportArgument(result,"maximum_cycle_value", MAXIMUM_CYCLE_VALUE, other.MAXIMUM_CYCLE_VALUE,thisRole,otherRole);
        compareSimpleReportArgument(result,"low_quality_tail", LOW_QUAL_TAIL, other.LOW_QUAL_TAIL,thisRole,otherRole);
        compareSimpleReportArgument(result,"default_platform", DEFAULT_PLATFORM, other.DEFAULT_PLATFORM,thisRole,otherRole);
        compareSimpleReportArgument(result,"force_platform", FORCE_PLATFORM, other.FORCE_PLATFORM,thisRole,otherRole);
        compareSimpleReportArgument(result,"quantizing_levels", QUANTIZING_LEVELS, other.QUANTIZING_LEVELS,thisRole,otherRole);
        compareSimpleReportArgument(result,"binary_tag_name", BINARY_TAG_NAME, other.BINARY_TAG_NAME,thisRole,otherRole);
        return result;
    }


    /**
     * Compares the covariate report lists.
     *
     * @param diffs map where to annotate the difference.
     * @param other the argument collection to compare against.
     * @param thisRole the name for this argument collection that makes sense to the user.
     * @param otherRole  the name for the other argument collection that makes sense to the end user.
     *
     * @return <code>true</code> if a difference was found.
     */
    @Requires("diffs != null && other != null && thisRole != null && otherRole != null")
    private boolean compareRequestedCovariates(final Map<String,String> diffs,
            final RecalibrationArgumentCollection other, final String thisRole, final String otherRole) {

        final Set<String> beforeNames = new HashSet<>(this.COVARIATES.length);
        final Set<String> afterNames = new HashSet<>(other.COVARIATES.length);
        Utils.addAll(beforeNames, this.COVARIATES);
        Utils.addAll(afterNames,other.COVARIATES);
        final Set<String> intersect = new HashSet<>(Math.min(beforeNames.size(),afterNames.size()));
        intersect.addAll(beforeNames);
        intersect.retainAll(afterNames);

        String diffMessage = null;
        if (intersect.size() == 0) { // In practice this is not possible due to required covariates but...
            diffMessage = String.format("There are no common covariates between '%s' and '%s'"
                    + " recalibrator reports. Covariates in '%s': {%s}. Covariates in '%s': {%s}.",thisRole,otherRole,
                    thisRole,Utils.join(", ",this.COVARIATES),
                    otherRole,Utils.join(",",other.COVARIATES));
        } else if (intersect.size() != beforeNames.size() || intersect.size() != afterNames.size()) {
            beforeNames.removeAll(intersect);
            afterNames.removeAll(intersect);
            diffMessage = String.format("There are differences in the set of covariates requested in the"
                    + " '%s' and '%s' recalibrator reports. "
                    + " Exclusive to '%s': {%s}. Exclusive to '%s': {%s}.",thisRole,otherRole,
                    thisRole,Utils.join(", ",beforeNames),
                    otherRole,Utils.join(", ",afterNames));
        }
        if (diffMessage != null) {
            diffs.put("covariate",diffMessage);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Annotates a map with any difference encountered in a simple value report argument that differs between this an
     * another {@link RecalibrationArgumentCollection} instance.
     * <p/>
     * The key of the new entry would be the name of that argument in the report file. The value is a message
     * that explains the difference to the end user.
     * <p/>
     *
     * <p/>
     * This method should not return any exception.
     *
     * @param diffs where to annotate the differences.
     * @param name the name of the report argument to compare.
     * @param thisValue this argument collection value for that argument.
     * @param otherValue the other collection value for that argument.
     * @param thisRole the name used to refer to this RAC report that makes sense to the end user.
     * @param otherRole the name used to refer to the other RAC report that makes sense to the end user.
     *
     * @type T the argument Object value type.
     *
     * @return <code>true</code> if a difference has been spotted, thus <code>diff</code> has been modified.
     */
    private <T> boolean compareSimpleReportArgument(final Map<String,String> diffs,
            final String name, final T thisValue, final T otherValue, final String thisRole, final String otherRole) {
        if (thisValue == null && otherValue == null) {
            return false;
        } else if (thisValue != null && thisValue.equals(otherValue)) {
            return false;
        } else {
            diffs.put(name,
                    String.format("differences between '%s' {%s} and '%s' {%s}.",
                            thisRole,thisValue == null ? "" : thisValue,
                            otherRole,otherValue == null ? "" : otherValue));
            return true;
        }

    }

    /**
     * Create a shallow copy of this argument collection.
     *
     * @return never <code>null</code>.
     */
    @Override
    public RecalibrationArgumentCollection clone() {
        try {
            return (RecalibrationArgumentCollection) super.clone();
        } catch (CloneNotSupportedException e) {
            throw new GATKException("Unreachable code clone not supported thrown when the class "
                    + this.getClass().getName() + " is cloneable ",e);
        }
    }

}

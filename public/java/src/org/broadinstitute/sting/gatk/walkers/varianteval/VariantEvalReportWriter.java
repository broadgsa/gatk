/*
 * Copyright (c) 2012, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTableV2;
import org.broadinstitute.sting.gatk.walkers.varianteval.evaluators.VariantEvaluator;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.manager.StratificationManager;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.Analysis;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.AnalysisModuleScanner;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.DataPoint;
import org.broadinstitute.sting.gatk.walkers.varianteval.util.EvaluationContext;
import org.broadinstitute.sting.utils.collections.Pair;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.StingException;

import java.io.PrintStream;
import java.lang.reflect.Field;
import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Class for writing the GATKReport for VariantEval
 *
 * Accepts a fulled evaluated (i.e., there's no more data coming) set of stratifications and evaluators
 * and supports writing out the data in these evaluators to a GATKReport.
 */
public class VariantEvalReportWriter {

    protected VariantEvalReportWriter() {}  // no public access

    /**
     * The business end of the class.  Writes out the data in the provided stratManager
     * to the PrintStream out
     *
     * @param out            the output stream
     * @param stratManager   the stratification manager
     * @param stratifiers    the stratifiers
     * @param evaluators     the evaluators
     */
    public static void writeReport(final PrintStream out,
                                   final StratificationManager<VariantStratifier, EvaluationContext> stratManager,
                                   final Collection<VariantStratifier> stratifiers,
                                   final Collection<VariantEvaluator> evaluators) {

        final GATKReport report = initializeGATKReport(stratifiers, evaluators);

        for ( int key = 0; key < stratManager.size(); key++ ) {
            final String stratStateString = stratManager.getStratsAndStatesStringForKey(key);
            final List<Pair<VariantStratifier, Object>> stratsAndStates = stratManager.getStratsAndStatesForKey(key);
            final EvaluationContext nec = stratManager.get(key);

            for ( final VariantEvaluator ve : nec.getVariantEvaluators() ) {
                final GATKReportTableV2 table = report.getTable(ve.getSimpleName());

                final AnalysisModuleScanner scanner = new AnalysisModuleScanner(ve);
                final Map<Field, DataPoint> datamap = scanner.getData();
                try {
                    if ( scanner.hasMoltenField() ) {
                        final Field field = scanner.getMoltenField();
                        final Object fieldValue = field.get(ve);

                        if ( fieldValue == null || ! (fieldValue instanceof Map) )
                            throw new ReviewedStingException("BUG field " + field.getName() + " must be a non-null instance of Map in " + scanner.getAnalysis().name());
                        final Map<Object, Object> map = (Map<Object, Object>)fieldValue;
                        if ( map.isEmpty() )
                            throw new ReviewedStingException("BUG: map is null or empty in analysis " + scanner.getAnalysis());
                        
                        int counter = 0; // counter is used to ensure printing order is as defined by entrySet
                        for ( Map.Entry<Object, Object> keyValue : map.entrySet() ) {
                            // "%05d" is a terrible hack to ensure sort order
                            final String moltenStratStateString = stratStateString + String.format("%05d", counter++);
                            setStratificationColumns(table, moltenStratStateString, stratsAndStates);
                            table.set(moltenStratStateString, scanner.getMoltenAnnotation().variableName(), keyValue.getKey());
                            table.set(moltenStratStateString, scanner.getMoltenAnnotation().valueName(), keyValue.getValue());
                        }
                    } else {
                        setStratificationColumns(table, stratStateString, stratsAndStates);
                        for ( final Field field : datamap.keySet()) {
                            table.set(stratStateString, field.getName(), field.get(ve));
                        }
                    }
                } catch (IllegalAccessException e) {
                    throw new ReviewedStingException("BUG: analysis field not public: " + e);
                }
            }
        }

        report.print(out);
    }

    /**
     * Common utility to configure a GATKReportTable columns
     *
     * Sets the column names to the strat names in stratsAndStates for the primary key in table
     *
     * @param table
     * @param primaryKey
     * @param stratsAndStates
     */
    private static void setStratificationColumns(final GATKReportTableV2 table,
                                                 final String primaryKey,
                                                 final List<Pair<VariantStratifier, Object>> stratsAndStates) {
        table.set(primaryKey, table.getTableName(), table.getTableName());
        for ( final Pair<VariantStratifier, Object> stratAndState : stratsAndStates ) {
            final VariantStratifier vs = stratAndState.getFirst();
            final String columnName = vs.getName();
            final Object strat = stratAndState.getSecond();
            if ( columnName == null || strat == null )
                throw new ReviewedStingException("Unexpected null variant stratifier state at " + table + " key = " + primaryKey);
            table.set(primaryKey, columnName, strat);
        }
    }

    /**
     * Initialize the output report
     *
     * We have a set of stratifiers and evaluation objects.  We need to create tables that look like:
     *
     * strat1 strat2 ... stratN eval1.field1 eval1.field2 ... eval1.fieldM
     *
     * for each eval.
     *
     * Note that this procedure doesn't support the creation of the old TableType system.  As the
     * VariantEvaluators are effectively tables themselves, we require authors to just create new
     * evaluation modules externally instead of allow them to embed them in other evaluation modules
     *
     * @return an initialized report object
     */
    private static GATKReport initializeGATKReport(final Collection<VariantStratifier> stratifiers,
                                                   final Collection<VariantEvaluator> evaluators) {
        final GATKReport report = new GATKReport();

        for (final VariantEvaluator ve : evaluators) {
            final AnalysisModuleScanner scanner = new AnalysisModuleScanner(ve);
            final Map<Field, DataPoint> datamap = scanner.getData();

            // create the table
            final String tableName = ve.getSimpleName();
            final String tableDesc = ve.getClass().getAnnotation(Analysis.class).description();
            report.addTable(tableName, tableDesc, 1 + stratifiers.size() + (scanner.hasMoltenField() ? 2 : datamap.size()), true);

            // grab the table, and add the columns we need to it
            final GATKReportTableV2 table = report.getTable(tableName);
            table.addColumn(tableName, tableName);

            // first create a column to hold each stratifier state
            for (final VariantStratifier vs : stratifiers) {
                final String columnName = vs.getName();
                table.addColumn(columnName, vs.getFormat());
            }

            if ( scanner.hasMoltenField() ) {
                // deal with molten data
                table.addColumn(scanner.getMoltenAnnotation().variableName(), scanner.getMoltenAnnotation().variableFormat());
                table.addColumn(scanner.getMoltenAnnotation().valueName(), scanner.getMoltenAnnotation().valueFormat());
            } else {
                if ( datamap.isEmpty() )
                    throw new ReviewedStingException("Datamap is empty for analysis " + scanner.getAnalysis());
                
                // add DataPoint's for each field marked as such
                for (final Field field : datamap.keySet()) {
                    try {
                        field.setAccessible(true);

                        // this is an atomic value, add a column for it
                        final String format = datamap.get(field).format();
                        table.addColumn(field.getName(), format);
                    } catch (SecurityException e) {
                        throw new StingException("SecurityException: " + e);
                    }
                }
            }
        }

        return report;
    }
}

/*
*  By downloading the PROGRAM you agree to the following terms of use:
*  
*  BROAD INSTITUTE - SOFTWARE LICENSE AGREEMENT - FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
*  
*  This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 (BROAD) and the LICENSEE and is effective at the date the downloading is completed (EFFECTIVE DATE).
*  
*  WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
*  WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
*  NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
*  
*  1. DEFINITIONS
*  1.1 PROGRAM shall mean copyright in the object code and source code known as GATK2 and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute/GATK on the EFFECTIVE DATE.
*  
*  2. LICENSE
*  2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM. 
*  The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
*  2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
*  2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.  
*  
*  3. OWNERSHIP OF INTELLECTUAL PROPERTY 
*  LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
*  Copyright 2012 Broad Institute, Inc.
*  Notice of attribution:  The GATK2 program was made available through the generosity of Medical and Population Genetics program at the Broad Institute, Inc.
*  LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
*  
*  4. INDEMNIFICATION
*  LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, (Indemnitees), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
*  
*  5. NO REPRESENTATIONS OR WARRANTIES
*  THE PROGRAM IS DELIVERED AS IS.  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
*  IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
*  
*  6. ASSIGNMENT
*  This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
*  
*  7. MISCELLANEOUS
*  7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
*  7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
*  7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
*  7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt. 
*  7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter. 
*  7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
*  7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
*/

package org.broadinstitute.sting.gatk.walkers.varianteval;

import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
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
                final GATKReportTable table = report.getTable(ve.getSimpleName());

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
    private static void setStratificationColumns(final GATKReportTable table,
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
            report.addTable(tableName, tableDesc, 1 + stratifiers.size() + (scanner.hasMoltenField() ? 2 : datamap.size()), GATKReportTable.TableSortingWay.SORT_BY_ROW);

            // grab the table, and add the columns we need to it
            final GATKReportTable table = report.getTable(tableName);
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
                for (final Map.Entry<Field, DataPoint> field : datamap.entrySet()) {
                    try {
                        field.getKey().setAccessible(true);

                        // this is an atomic value, add a column for it
                        final String format = field.getValue().format();
                        table.addColumn(field.getKey().getName(), format);
                    } catch (SecurityException e) {
                        throw new StingException("SecurityException: " + e);
                    }
                }
            }
        }

        return report;
    }
}

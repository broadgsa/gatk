/*
 * Copyright (c) 2011, The Broad Institute
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

package org.broadinstitute.sting.gatk.walkers.diffengine;

import com.google.java.contract.Requires;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.gatk.report.GATKReport;
import org.broadinstitute.sting.gatk.report.GATKReportTable;
import org.broadinstitute.sting.gatk.walkers.varianteval.stratifications.VariantStratifier;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.PluginManager;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: 7/4/11
 * Time: 12:51 PM
 * A generic engine for comparing tree-structured objects
 */
public class DiffEngine {
    final protected static Logger logger = Logger.getLogger(DiffEngine.class);

    private final Map<String, DiffableReader> readers = new HashMap<String, DiffableReader>();

    public DiffEngine() {
        loadDiffableReaders();
    }

    // --------------------------------------------------------------------------------
    //
    // difference calculation
    //
    // --------------------------------------------------------------------------------

    public List<Difference> diff(DiffElement master, DiffElement test) {
        DiffValue masterValue = master.getValue();
        DiffValue testValue = test.getValue();

        if ( masterValue.isCompound() && masterValue.isCompound() ) {
            return diff(master.getValueAsNode(), test.getValueAsNode());
        } else if ( masterValue.isAtomic() && testValue.isAtomic() ) {
            return diff(masterValue, testValue);
        } else {
            // structural difference in types.  one is node, other is leaf
            return Arrays.asList(new Difference(master, test));
        }
    }

    public List<Difference> diff(DiffNode master, DiffNode test) {
        Set<String> allNames = new HashSet<String>(master.getElementNames());
        allNames.addAll(test.getElementNames());
        List<Difference> diffs = new ArrayList<Difference>();

        for ( String name : allNames ) {
            DiffElement masterElt = master.getElement(name);
            DiffElement testElt = test.getElement(name);
            if ( masterElt == null && testElt == null ) {
                throw new ReviewedStingException("BUG: unexceptedly got two null elements for field: " + name);
            } else if ( masterElt == null || testElt == null ) { // if either is null, we are missing a value
                // todo -- should one of these be a special MISSING item?
                diffs.add(new Difference(masterElt, testElt));
            } else {
                diffs.addAll(diff(masterElt, testElt));
            }
        }

        return diffs;
    }

    public List<Difference> diff(DiffValue master, DiffValue test) {
        if ( master.getValue().equals(test.getValue()) ) {
            return Collections.emptyList();
        } else {
            return Arrays.asList(new Difference(master.getBinding(), test.getBinding()));
        }
    }

    // --------------------------------------------------------------------------------
    //
    // Summarizing differences
    //
    // --------------------------------------------------------------------------------

    /**
     * Emits a summary of the diffs to out.  Suppose you have the following three differences:
     *
     *   A.X.Z:1!=2
     *   A.Y.Z:3!=4
     *   B.X.Z:5!=6
     *
     * The above is the itemized list of the differences.  The summary looks for common differences
     * in the name hierarchy, counts those shared elements, and emits the differences that occur
     * in order of decreasing counts.
     *
     * So, in the above example, what are the shared elements?
     *
     * A.X.Z and B.X.Z share X.Z, so there's a *.X.Z with count 2
     * A.X.Z, A.Y.Z, and B.X.Z all share *.*.Z, with count 3
     * Each of A.X.Z, A.Y.Z, and B.X.Z are individually unique, with count 1
     *
     * So we would emit the following summary:
     *
     * *.*.Z: 3
     * *.X.Z: 2
     * A.X.Z: 1 [specific difference: 1!=2]
     * A.Y.Z: 1 [specific difference: 3!=4]
     * B.X.Z: 1 [specific difference: 5!=6]
     *
     * The algorithm to accomplish this calculation is relatively simple. Start with all of the
     * concrete differences.  For each pair of differences A1.A2....AN and B1.B2....BN:
     *
     * find the longest common subsequence Si.Si+1...SN where Ai = Bi = Si
     * If i == 0, then there's no shared substructure
     * If i > 0, then generate the summarized value X = *.*...Si.Si+1...SN
     * if X is a known summary, increment it's count, otherwise set its count to 1
     *
     * Not that only pairs of the same length are considered as potentially equivalent
     *
     * @param params determines how we display the items
     * @param diffs
     */
    public void reportSummarizedDifferences(List<Difference> diffs, SummaryReportParams params ) {
        printSummaryReport(summarizeDifferences(diffs), params );
    }

    public List<SummarizedDifference> summarizeDifferences(List<Difference> diffs) {
        List<String[]> diffPaths = new ArrayList<String[]>(diffs.size());

        for ( Difference diff1 : diffs ) {
            diffPaths.add(diffNameToPath(diff1.getFullyQualifiedName()));
        }

        return summarizedDifferencesOfPaths(diffPaths);
    }

    final protected static String[] diffNameToPath(String diffName) {
        return diffName.split("\\.");
    }

    protected List<SummarizedDifference> summarizedDifferencesOfPaths(List<String[]> diffPaths) {
        Map<String, SummarizedDifference> summaries = new HashMap<String, SummarizedDifference>();

        // create the initial set of differences
        for ( int i = 0; i < diffPaths.size(); i++ ) {
            for ( int j = 0; j <= i; j++ ) {
                String[] diffPath1 = diffPaths.get(i);
                String[] diffPath2 = diffPaths.get(j);
                if ( diffPath1.length == diffPath2.length ) {
                    int lcp = longestCommonPostfix(diffPath1, diffPath2);
                    String path = lcp > 0 ? summarizedPath(diffPath2, lcp) : Utils.join(".", diffPath2);
                    addSummary(summaries, path, true);
                }
            }
        }

        // count differences
        for ( String[] diffPath : diffPaths ) {
            for ( SummarizedDifference sumDiff : summaries.values() ) {
                if ( sumDiff.matches(diffPath) )
                    addSummary(summaries, sumDiff.getPath(), false);
            }
        }

        List<SummarizedDifference> sortedSummaries = new ArrayList<SummarizedDifference>(summaries.values());
        Collections.sort(sortedSummaries);
        return sortedSummaries;
    }

    private static void addSummary(Map<String, SummarizedDifference> summaries, String path, boolean onlyCatalog) {
        if ( summaries.containsKey(path) ) {
            if ( ! onlyCatalog )
                summaries.get(path).incCount();
        } else {
            SummarizedDifference sumDiff = new SummarizedDifference(path);
            summaries.put(sumDiff.getPath(), sumDiff);
        }
    }

    protected void printSummaryReport(List<SummarizedDifference> sortedSummaries, SummaryReportParams params ) {
        GATKReport report = new GATKReport();
        final String tableName = "diffences";
        report.addTable(tableName, "Summarized differences between the master and test files.\nSee http://www.broadinstitute.org/gsa/wiki/index.php/DiffObjectsWalker_and_SummarizedDifferences for more information");
        GATKReportTable table = report.getTable(tableName);
        table.addPrimaryKey("Difference", true);
        table.addColumn("NumberOfOccurrences", 0);

        int count = 0, count1 = 0;
        for ( SummarizedDifference diff : sortedSummaries ) {
            if ( diff.getCount() < params.minSumDiffToShow )
                // in order, so break as soon as the count is too low
                break;

            if ( params.maxItemsToDisplay != 0 && count++ > params.maxItemsToDisplay )
                break;

            if ( diff.getCount() == 1 ) {
                count1++;
                if ( params.maxCountOneItems != 0 && count1 > params.maxCountOneItems )
                    break;
            }

            table.set(diff.getPath(), "NumberOfOccurrences", diff.getCount());
        }

        table.write(params.out);
    }

    protected static int longestCommonPostfix(String[] diffPath1, String[] diffPath2) {
        int i = 0;
        for ( ; i < diffPath1.length; i++ ) {
            int j = diffPath1.length - i - 1;
            if ( ! diffPath1[j].equals(diffPath2[j]) )
                break;
        }
        return i;
    }

    /**
     * parts is [A B C D]
     * commonPostfixLength: how many parts are shared at the end, suppose its 2
     * We want to create a string *.*.C.D
     *
     * @param parts
     * @param commonPostfixLength
     * @return
     */
    protected static String summarizedPath(String[] parts, int commonPostfixLength) {
        int stop = parts.length - commonPostfixLength;
        if ( stop > 0 ) parts = parts.clone();
        for ( int i = 0; i < stop; i++ ) {
            parts[i] = "*";
        }
        return Utils.join(".", parts);
    }

    /**
     * TODO -- all of the algorithms above should use SummarizedDifference instead
     * TODO -- of some SummarizedDifferences and some low-level String[]
     */
    public static class SummarizedDifference implements Comparable<SummarizedDifference> {
        final String path; // X.Y.Z
        final String[] parts;
        int count = 0;

        public SummarizedDifference(String path) {
            this.path = path;
            this.parts = diffNameToPath(path);
        }

        public void incCount() { count++; }

        public int getCount() {
            return count;
        }

        /**
         * The fully qualified path object A.B.C etc
         * @return
         */
        public String getPath() {
            return path;
        }

        /**
         * @return the length of the parts of this summary
         */
        public int length() {
            return this.parts.length;
        }

        /**
         * Returns true if the string parts matches this summary.  Matches are
         * must be equal() everywhere where this summary isn't *.
         * @param otherParts
         * @return
         */
        public boolean matches(String[] otherParts) {
            if ( otherParts.length != length() )
                return false;

            // TODO optimization: can start at right most non-star element
            for ( int i = 0; i < length(); i++ ) {
                String part = parts[i];
                if ( ! part.equals("*") && ! part.equals(otherParts[i]) )
                    return false;
            }

            return true;
        }

        @Override
        public String toString() {
            return String.format("%s:%d", getPath(), getCount());
        }

        @Override
        public int compareTo(SummarizedDifference other) {
            // sort first highest to lowest count, then by lowest to highest path
            int countCmp = Integer.valueOf(count).compareTo(other.count);
            return countCmp != 0 ? -1 * countCmp : path.compareTo(other.path);
        }


    }

    // --------------------------------------------------------------------------------
    //
    // plugin manager
    //
    // --------------------------------------------------------------------------------

    public void loadDiffableReaders() {
        List<Class<? extends DiffableReader>> drClasses = new PluginManager<DiffableReader>( DiffableReader.class ).getPlugins();

        logger.info("Loading diffable modules:");
        for (Class<? extends DiffableReader> drClass : drClasses ) {
            logger.info("\t" + drClass.getSimpleName());

            try {
                DiffableReader dr = drClass.newInstance();
                readers.put(dr.getName(), dr);
            } catch (InstantiationException e) {
                throw new ReviewedStingException("Unable to instantiate module '" + drClass.getSimpleName() + "'");
            } catch (IllegalAccessException e) {
                throw new ReviewedStingException("Illegal access error when trying to instantiate '" + drClass.getSimpleName() + "'");
            }
        }
    }

    protected Map<String, DiffableReader> getReaders() {
        return readers;
    }

    protected DiffableReader getReader(String name) {
        return readers.get(name);
    }

    /**
     * Returns a reader appropriate for this file, or null if no such reader exists
     * @param file
     * @return
     */
    public DiffableReader findReaderForFile(File file) {
        for ( DiffableReader reader : readers.values() )
            if (reader.canRead(file) )
                return reader;

        return null;
    }

    /**
     * Returns true if reader appropriate for this file, or false if no such reader exists
     * @param file
     * @return
     */
    public boolean canRead(File file) {
        return findReaderForFile(file) != null;
    }


    public DiffElement createDiffableFromFile(File file) {
        return createDiffableFromFile(file, -1);
    }

    public DiffElement createDiffableFromFile(File file, int maxElementsToRead) {
        DiffableReader reader = findReaderForFile(file);
        if ( reader == null )
            throw new UserException("Unsupported file type: " + file);
        else
            return reader.readFromFile(file, maxElementsToRead);
    }

    public static boolean simpleDiffFiles(File masterFile, File testFile, DiffEngine.SummaryReportParams params) {
        DiffEngine diffEngine = new DiffEngine();

        if ( diffEngine.canRead(masterFile) && diffEngine.canRead(testFile) ) {
            DiffElement master = diffEngine.createDiffableFromFile(masterFile);
            DiffElement test = diffEngine.createDiffableFromFile(testFile);
            List<Difference> diffs = diffEngine.diff(master, test);
            diffEngine.reportSummarizedDifferences(diffs, params);
            return true;
        } else {
            return false;
        }
    }

    public static class SummaryReportParams {
        PrintStream out = System.out;
        int maxItemsToDisplay = 0;
        int maxCountOneItems = 0;
        int minSumDiffToShow = 0;

        public SummaryReportParams(PrintStream out, int maxItemsToDisplay, int maxCountOneItems, int minSumDiffToShow) {
            this.out = out;
            this.maxItemsToDisplay = maxItemsToDisplay;
            this.maxCountOneItems = maxCountOneItems;
            this.minSumDiffToShow = minSumDiffToShow;
        }
    }
}

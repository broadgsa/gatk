/*
 * Copyright (c) 2010 The Broad Institute
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the ”Software”), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED ”AS IS”, WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 * THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.playground.analyzeconcordance;

import org.broadinstitute.sting.commandline.CommandLineProgram;
import org.broadinstitute.sting.commandline.Argument;
import org.broadinstitute.sting.utils.text.XReadLines;
import org.broadinstitute.sting.utils.StingException;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.PathUtils;
import org.broadinstitute.sting.playground.utils.ProcessUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.util.*;

/**
 * Compares results of VariantEval across a population or a case/control group.
 */
public class AnalyzeConcordance extends CommandLineProgram {

	@Argument(fullName = "group_name", shortName = "groupName", doc = "The name of the group which will be prefixed output files", required = false)
	private String baseName = "analyze_concordance";
	@Argument(fullName = "eval_list", shortName = "evalList", doc = "The input list of unfiltered eval files to analyze", required = true)
	private String evalListFile = null;
    @Argument(fullName = "filtered_eval_list", shortName = "filteredEvalList", doc = "The input list of filtered eval files to analyze", required = false)
    private String filteredEvalListFile = null;
    @Argument(fullName = "output_dir", shortName = "outputDir", doc = "The directory in which to output all the plots and intermediate data files", required = false)
    private String outputDir = "analyzeConcordance";
    @Argument(fullName = "path_to_Rscript", shortName = "Rscript", doc = "The path to your implementation of Rscript", required = false)
    private String pathToRscript = "env Rscript";
    @Argument(fullName = "path_to_resources", shortName = "resources", doc = "Path to resources folder holding the Sting analyze concordance R scripts", required = false)
    private String pathToResources = "R" + File.separator + "analyzeConcordance";

    private enum EvalFilterType {
        UNFILTERED, FILTERED
    }

    private static final AnalyzeConcordanceField[] ANALYZE_CONCORDANCE_FIELDS = AnalyzeConcordanceField.values();

	private String evalDataFile;
    private List<String[]> data = new ArrayList<String[]>();

    private static Logger logger = Logger.getLogger(AnalyzeConcordance.class);

    protected int execute() {
        int result;

        try {
            createOutputDirectory();

            // initialize all the data from the csv file and allocate the list of covariates
            logger.info("Reading in input csv file...");
            initializeData();
            logger.info("...Done!");

            // output data tables for Rscript to read in
            logger.info("Writing out intermediate tables for R...");
            writeDataTables();
            logger.info("...Done!");

			// perform the analysis using Rscript and output the plots
			logger.info("Calling analysis R scripts and writing out figures...");
			result = callRScripts();
			logger.info("...Done!");

			// perform the analysis using Rscript and output the plots
			logger.info("Generating html report...");
			generateHtmlReport();
			logger.info("...Done!");
			
        } catch (StingException se) {
            throw se;
        } catch (Exception e) {
            throw new StingException("Error analyzing concordance", e);
        }

        return result;
    }

    private void createOutputDirectory() {
        // create the output directory where all the data tables and plots will go
        File outputDir = new File(this.outputDir);
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new StingException("Couldn't create directory: " + this.outputDir);
        }
    }

    private void initializeData() throws FileNotFoundException {
        // add the column headers to the data
        addHeader();

        // read the list of unfiltered eval files
        addEvalListFile(EvalFilterType.UNFILTERED, new File(evalListFile));

        // if provided, read the list of filtered eval files
        if (filteredEvalListFile != null) {
            addEvalListFile(EvalFilterType.FILTERED, new File(filteredEvalListFile));
        }
    }

    private void addHeader() {
        String[] headers = new String[ANALYZE_CONCORDANCE_FIELDS.length + 2];
        int column = 0;
        headers[column++] = "eval_id";
        headers[column++] = "filter_type";

        for (AnalyzeConcordanceField field : ANALYZE_CONCORDANCE_FIELDS) {
            headers[column++] = field.getColumnHeader();
        }

        data.add(headers);
    }

    private void addEvalListFile(EvalFilterType filterType, File evalListFile) throws FileNotFoundException {
        for (String line : new XReadLines(evalListFile)) {
            String[] parts = line.split("\t");
            addEvalFile(parts[0], filterType, new File(parts[1]));
        }
    }

    private void addEvalFile(String evalID, EvalFilterType filterType, File evalFile) throws FileNotFoundException {
        SortedMap<AnalyzeConcordanceField, String> fieldValues = new TreeMap<AnalyzeConcordanceField, String>();

        for (String line : new XReadLines(evalFile)) {
            for (AnalyzeConcordanceField field : ANALYZE_CONCORDANCE_FIELDS) {
                String value = field.parseLine(line);
                if (value != null) {
                    fieldValues.put(field, value);
                    break;  // continue to the next line.
                }
            }
        }

        String[] values = new String[ANALYZE_CONCORDANCE_FIELDS.length + 2];
        int column = 0;
        values[column++] = evalID;
        values[column++] = filterType.toString().toLowerCase();

        // get all the values, including null if for some reason a value wasn't found
        for (AnalyzeConcordanceField field : ANALYZE_CONCORDANCE_FIELDS) {
            values[column++] = fieldValues.get(field);
        }

        data.add(values);
    }

    private void writeDataTables() throws FileNotFoundException {
		evalDataFile = baseName + ".eval_data.tsv";
        // Create a PrintStream
        PrintStream output = new PrintStream(new File(outputDir, evalDataFile));
        for (String[] line : data) {
            output.println(Utils.join("\t", line));
        }
        output.close();
    }

    private int callRScripts() {
        String command = pathToRscript + " "
                + new File(pathToResources, "analyzeConcordance.R") + " "
				+ new File(outputDir, baseName) + " "
                + new File(outputDir, evalDataFile);

        return ProcessUtils.runCommandAndWait(command);
    }

	private void generateHtmlReport() throws FileNotFoundException {
		// TODO: Enhance the reports
		PrintStream output = new PrintStream(new File(outputDir, "report.html"));
		output.println("<html><body>");
		for (File pngFile : new File(outputDir).listFiles(new PathUtils.ExtensionFilter("png"))) {
			output.println("<div><img src=\"" + pngFile.getName() + "\"/></div>");
		}
		output.println("</body></html>");
		output.close();
	}

    public static void main(String[] argv) {
        try {
            AnalyzeConcordance instance = new AnalyzeConcordance();
            start(instance, argv);
            System.exit(CommandLineProgram.result);
        } catch (Exception e) {
            exitSystemWithError(e);
        }
    }
}

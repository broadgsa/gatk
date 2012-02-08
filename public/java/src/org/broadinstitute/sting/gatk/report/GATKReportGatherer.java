package org.broadinstitute.sting.gatk.report;

import org.broadinstitute.sting.commandline.Gatherer;
import org.broadinstitute.sting.utils.exceptions.UserException;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: roger
 * Date: 1/9/12
 * Time: 11:17 PM
 * To change this template use File | Settings | File Templates.
 */
public class GATKReportGatherer extends Gatherer {
    @Override
    public void gather(List<File> inputs, File output) {
        //Combines inputs GATKReport to one output

        PrintStream o;
        try {
            o = new PrintStream(output);
        } catch (FileNotFoundException e) {
            throw new UserException("File to be output by CoverageByRG Gather function was not found");
        }

        GATKReport current = new GATKReport();
        boolean isFirst = true;
        for (File input : inputs) {

            // If the table is empty
            if (isFirst) {
                current = new GATKReport(input);
                isFirst = false;
            } else {
                GATKReport toAdd = new GATKReport(input);
                current.combineWith(toAdd);
            }
        }

        current.print(o);
    }
}

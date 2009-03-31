package org.broadinstitute.sting.gatk.dataSources;

import org.broadinstitute.sting.gatk.walkers.Walker;
import org.broadinstitute.sting.gatk.walkers.ReadWalker;

import java.util.ArrayList;
import java.io.File;

/**
 * User: aaron
 * Date: Mar 25, 2009
 * Time: 4:51:39 PM
 * <p/>
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2009 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 * <p/>
 * This software is supplied without any warranty or guaranteed support whatsoever. Neither
 * the Broad Institute nor MIT can be responsible for its use, misuse, or functionality.
 */
public class DataSourceBuilder {

    // storage for the passed file
    ArrayList<File> passFiles = new ArrayList<File>();

    public DataSourceBuilder() {

    }

    /**
     * add a file used to generate the data sources
     *
     * @param fileName the filename that should be used
     */
    public void addDataFile(String fileName) {
        // for now, just add it to the internal file list
        passFiles.add(new File(fileName));
    }

    /**
     * add a file used to generate the data sources
     *
     * @param file the filename that should be used
     */
    public void addDataFile(File file) {
        // for now, just add it to the internal file list
        passFiles.add(file);
    }

    public DataSource build(Walker inputWalker) {
        if (inputWalker instanceof ReadWalker) {

        }

        return null;
    }


    /**
     * this section contains the private methods to create data sources
     * based on the type of walker we're passed in.
     */


    /**
     * we know we have a read data source, let's get the 
     * @return
     */
    //private ReadDataSource generateReadDataSource() {
        // 
    //}

}

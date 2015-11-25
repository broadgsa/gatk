/*
* Copyright 2012-2015 Broad Institute, Inc.
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

package org.broadinstitute.gatk.engine.filters;

import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.text.XReadLines;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.regex.Pattern;

/**
 * This is a utility class, its sole purpose is to populate PlatformUnitFilter with data. When a command line argument
 * (@Argument) of the type PlatformUnitFilterHelper is declared in an application (walker), its constuctor
 * PlatformUnitFilterHelper(String) automatically called by the argument system will parse its String argument
 * and set up static fields of PlatformUnitFilter object.
 *
 * The String argument can be either a name of existing file, or a list of comma-separated lane (Platform Unit) names.
 * First, the constructor will check if a file with specified name exists. If it does, then it is assumed that each line
 * in the file contains one name of a lane (Platfor Unit) to filter out. If such file does not exist, then the argument is
 * interpreted as a comma-separated list. Blank spaces around lane names are allowed in both cases and will be trimmed out.
 *
 * In other words, all it takes to request filtering out reads from specific lane(s) is
 *
 * 1) declare filter usage in the walker
 *
 * @ReadFilters({PlatformUnitFilter.class,...})
 *
 * 2) specify the argument that will take the list of lanes to filter:
 *
 * @Argument(fullName="filterLanes", shortName="FL", doc="all specified lanes will be ignored", required=false)
 *   PlatformUnitFilterHelper dummy;
 *
 * After that, the walker can be invoked with "--filterLanes 302UBAAXX090508.8,302YAAAXX090427.8" argument.
 *
 * Created by IntelliJ IDEA.
 * User: asivache
 * Date: Sep 22, 2009
 * Time: 11:11:48 AM
 * To change this template use File | Settings | File Templates.
 */
public class PlatformUnitFilterHelper {
    final public static Pattern EMPTYLINE_PATTERN = Pattern.compile("^\\s*$");

    public PlatformUnitFilterHelper(String arg) {
         File f = new File(arg);

         if ( f.exists() ) {
             try {
                 XReadLines reader = new XReadLines(f);
                 for ( String line : reader ) {
                     if ( EMPTYLINE_PATTERN.matcher(line).matches() ) continue; // skip empty lines
                     PlatformUnitFilter.addBlackListedLane(line); // PlatformUnitFilter will trim the line as needed
                 }
             } catch ( FileNotFoundException e) { throw new UserException.CouldNotReadInputFile(f, e); } // this should NEVER happen
             return;
         }

        // no such file, must be a comma-separated list:

        PlatformUnitFilter.setBlackListedLanes(arg); // PlatformUnitFilter will split on commas and trim as needed

    }
}

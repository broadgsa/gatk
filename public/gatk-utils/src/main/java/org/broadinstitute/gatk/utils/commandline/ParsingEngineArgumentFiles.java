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

package org.broadinstitute.gatk.utils.commandline;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.SortedMap;

/**
 * Container class to store the list of argument files.
 * The files will be parsed after the command line arguments.
 */
public class ParsingEngineArgumentFiles extends ParsingEngineArgumentProvider {
    @Argument(fullName = "arg_file", shortName = "args", doc = "Reads arguments from the specified file", required = false)
    public List<File> files = new ArrayList<File>();

    @Override
    public void parse(ParsingEngine parsingEngine, SortedMap<ArgumentMatchSource, ParsedArgs> parsedArgs) {
        ArgumentMatches argumentMatches = parsingEngine.getArgumentMatches();
        for (File file: this.files) {
            List<String> fileTokens = parsingEngine.getArguments(file);
            parsingEngine.parse(new ArgumentMatchFileSource(file), fileTokens, argumentMatches, parsedArgs);
        }
    }
}

class ArgumentMatchFileSource extends ArgumentMatchSource {
    ArgumentMatchFileSource(File file) {
        super("file " + file.getAbsolutePath());
    }
}

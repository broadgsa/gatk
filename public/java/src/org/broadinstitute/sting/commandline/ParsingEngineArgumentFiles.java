package org.broadinstitute.sting.commandline;

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

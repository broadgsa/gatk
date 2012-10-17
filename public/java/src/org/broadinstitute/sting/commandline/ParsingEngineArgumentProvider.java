package org.broadinstitute.sting.commandline;

import java.util.List;
import java.util.SortedMap;

/**
 * A class that can parse arguments for the engine
 */
public abstract class ParsingEngineArgumentProvider {
    public abstract void parse(ParsingEngine parsingEngine, SortedMap<ArgumentMatchSource, ParsedArgs> parsedArgs);
}


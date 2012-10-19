package org.broadinstitute.sting.commandline;

/**
 * Represents a collection of parsed arguments for an argument source.
 *
 * Useful for printing out help documents.
 */
public abstract class ParsedArgs {
    /**
     * @return A compact description of the arguments from an provider/source.
     */
    public abstract String getDescription();
}

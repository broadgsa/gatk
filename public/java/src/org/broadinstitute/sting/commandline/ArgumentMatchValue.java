package org.broadinstitute.sting.commandline;

import java.io.File;

/**
 * Returns argument values as either strings or values.
 */
public abstract class ArgumentMatchValue {
    /**
     * @return the value of this argument as a String object.
     */
    public abstract String asString();

    /**
     * @return the value of this argument as a File object.
     */
    public abstract File asFile();
}

package org.broadinstitute.sting.commandline;

import java.io.File;

/**
 * Argument values that originated from a string.
 */
public class ArgumentMatchStringValue extends ArgumentMatchValue {
    private final String value;

    public ArgumentMatchStringValue(String value) {
        this.value = value;
    }

    @Override
    public String asString() {
        return value;
    }

    @Override
    public File asFile() {
        return value == null ? null : new File(value);
    }
}

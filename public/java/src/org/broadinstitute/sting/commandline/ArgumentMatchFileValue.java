package org.broadinstitute.sting.commandline;

import java.io.File;

/**
 * Holds a reference to a file as an argument match value.
 *
 * This is useful when the type of the stored file may be a subclass of java.io.File,
 * for example a Queue RemoteFile.
 */
public class ArgumentMatchFileValue extends ArgumentMatchValue {
    private final File file;

    public ArgumentMatchFileValue(File file) {
        this.file = file;
    }

    @Override
    public String asString() {
        return file == null ? null : file.getAbsolutePath();
    }

    @Override
    public File asFile() {
        return file;
    }
}

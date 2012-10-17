package org.broadinstitute.sting.commandline;

import org.apache.commons.lang.StringUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A list of string arguments, usually from the command line or an args list file.
 */
public class ParsedListArgs extends ParsedArgs {
    private final List<String> args = new ArrayList<String>();

    public ParsedListArgs() {
    }

    public ParsedListArgs(List<String> args) {
        this.args.addAll(args);
    }

    public void add(String... args) {
        this.args.addAll(Arrays.asList(args));
    }

    @Override
    public String getDescription() {
        return StringUtils.join(this.args, " ");
    }
}

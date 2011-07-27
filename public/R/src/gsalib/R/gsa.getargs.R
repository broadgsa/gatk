.gsa.getargs.usage <- function(argspec, doc) {
    cargs = commandArgs();

    usage = "Usage:";

    fileIndex = grep("--file=", cargs);
    if (length(fileIndex) > 0) {
        progname = gsub("--file=", "", cargs[fileIndex[1]]);

        usage = sprintf("Usage: Rscript %s [arguments]", progname);

        if (!is.na(doc)) {
            message(sprintf("%s: %s\n", progname, doc));
        }
    }

    message(usage);

    for (argname in names(argspec)) {
        key = argname;
        defaultValue = 0;
        doc = "";

        if (is.list(argspec[[argname]])) {
            defaultValue = argspec[[argname]]$value;
            doc = argspec[[argname]]$doc;
        }

        message(sprintf(" -%-10s\t[default: %s]\t%s", key, defaultValue, doc));
    }

    message("");

    stop(call. = FALSE);
}

gsa.getargs <- function(argspec, doc = NA) {
    argsenv = new.env();

    for (argname in names(argspec)) {
        value = 0;
        if (is.list(argspec[[argname]])) {
            value = argspec[[argname]]$value;
        } else {
            value = argspec[[argname]];
        }

        assign(argname, value, envir=argsenv);
    }

    if (interactive()) {
        for (argname in names(argspec)) {
            value = get(argname, envir=argsenv);

            if (is.na(value) | is.null(value)) {
                if (exists("cmdargs")) {
                    assign(argname, cmdargs[[argname]], envir=argsenv);
                } else {
                    assign(argname, readline(sprintf("Please enter a value for '%s': ", argname)), envir=argsenv);
                }
            } else {
                assign(argname, value, envir=argsenv);
            }
        }
    } else {
        cargs = commandArgs(TRUE);

        if (length(cargs) == 0) {
            .gsa.getargs.usage(argspec, doc);
        }

        for (i in 1:length(cargs)) {
            if (length(grep("^-", cargs[i], ignore.case=TRUE)) > 0) {
                key = gsub("-", "", cargs[i]);
                value = cargs[i+1];

                if (key == "h" | key == "help") {
                    .gsa.getargs.usage(argspec, doc);
                }

                if (length(grep("^[\\d\\.e\\+\\-]+$", value, perl=TRUE, ignore.case=TRUE)) > 0) {
                    value = as.numeric(value);
                }
                
                assign(key, value, envir=argsenv);
            }
        }
    }

    args = as.list(argsenv);

    isMissingArgs = 0;
    missingArgs = c();
    
    for (arg in names(argspec)) {
        if (is.na(args[[arg]]) | is.null(args[[arg]])) {
            gsa.warn(sprintf("Value for required argument '-%s' was not specified", arg));

            isMissingArgs = 1;
            missingArgs = c(missingArgs, arg);
        }
    }

    if (isMissingArgs) {
        gsa.error(
            paste(
                "Missing required arguments: -",
                paste(missingArgs, collapse=" -"),
                ".  Specify -h or -help to this script for a list of available arguments.",
                sep=""
            )
        );
    }

    args;
}

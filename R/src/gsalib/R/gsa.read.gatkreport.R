# Load a table into the specified environment.  Make sure that each new table gets a unique name (this allows one to cat a bunch of tables with the same name together and load them into R without each table overwriting the last.
.gsa.assignGATKTableToEnvironment <- function(tableName, tableHeader, tableRows, tableEnv) {
    d = data.frame(tableRows, row.names=NULL, stringsAsFactors=FALSE);
    colnames(d) = tableHeader;

    for (i in 1:ncol(d)) {
        v = suppressWarnings(as.numeric(d[,i]));

        if (length(na.omit(as.numeric(v))) == length(d[,i])) {
            d[,i] = v;
        }
    }

    usedNames = ls(envir=tableEnv, pattern=tableName);

    if (length(usedNames) > 0) {
        tableName = paste(tableName, ".", length(usedNames), sep="");
    }

    assign(tableName, d, envir=tableEnv);
}

# Load all GATKReport tables from a file
gsa.read.gatkreport <- function(filename) {
    con = file(filename, "r", blocking = TRUE);
    lines = readLines(con);
    close(con);

    tableEnv = new.env();

    tableName = NA;
    tableHeader = c();
    tableRows = c();

    for (line in lines) {
        if (length(grep("^##:GATKReport.v0.1[[:space:]]+", line, ignore.case=TRUE)) > 0) {
            headerFields = unlist(strsplit(line, "[[:space:]]+"));

            if (!is.na(tableName)) {
                .gsa.assignGATKTableToEnvironment(tableName, tableHeader, tableRows, tableEnv);
            }

            tableName = headerFields[2];
            tableHeader = c();
            tableRows = c();
        } else if (length(grep("^[[:space:]]*$", line)) > 0 | length(grep("^[[:space:]]*#", line)) > 0) {
            # do nothing
        } else if (!is.na(tableName)) {
            row = unlist(strsplit(line, "[[:space:]]+"));

            if (length(tableHeader) == 0) {
                tableHeader = row;
            } else {
                tableRows = rbind(tableRows, row);
            }
        }
    }

    if (!is.na(tableName)) {
        .gsa.assignGATKTableToEnvironment(tableName, tableHeader, tableRows, tableEnv);
    }

    gatkreport = as.list(tableEnv);
}

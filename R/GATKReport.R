read.gatkreport <- function(filename) {
    con = file(filename, "r", blocking = TRUE);
    lines = readLines(con);
    close(con);

    tableEnv = new.env();

    tableName = NA;
    tableHeader = c();
    tableRows = c();

    for (line in lines) {
        if (length(grep("^#:table[[:space:]]+", line, ignore.case=TRUE)) > 0) {
            headerFields = unlist(strsplit(line, "[[:space:]]+"));

            if (!is.na(tableName)) {
                d = data.frame(tableRows, row.names=NULL, stringsAsFactors=FALSE);
                colnames(d) = tableHeader;

                assign(tableName, d, envir=tableEnv);
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
        d = data.frame(tableRows, row.names=NULL);
        colnames(d) = tableHeader;

        assign(tableName, d, envir=tableEnv);
    }

    gatkreport = as.list(tableEnv);
}

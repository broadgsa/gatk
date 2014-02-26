gsa.read.squidmetrics = function(project, bylane = FALSE) {
    suppressMessages(library(ROracle));

    drv = dbDriver("Oracle");
    con = dbConnect(drv, "REPORTING/REPORTING@ora01:1521/SEQPROD");

    if (bylane) {
        statement = paste("SELECT * FROM ILLUMINA_PICARD_METRICS WHERE \"Project\" = '", project, "'", sep="");
        print(statement);

        rs  = dbSendQuery(con, statement = statement);
        d = fetch(rs, n=-1);
        dbHasCompleted(rs);
        dbClearResult(rs);
    } else {
        statement = paste("SELECT * FROM ILLUMINA_SAMPLE_STATUS_AGG WHERE \"Project\" = '", project, "'", sep="");
        print(statement);

        rs = dbSendQuery(con, statement = statement);
        d = fetch(rs, n=-1);
        dbHasCompleted(rs);
        dbClearResult(rs);
    }

    oraCloseDriver(drv);

    subset(d, Project == project);
}

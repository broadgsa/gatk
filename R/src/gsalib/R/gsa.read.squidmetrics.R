gsa.read.squidmetrics = function(project, bylane = FALSE) {
    suppressMessages(library(ROracle));

    drv = dbDriver("Oracle");
    con = dbConnect(drv, "REPORTING/REPORTING@ora01:1521/SEQPROD");

    if (bylane) {
        rs  = dbSendQuery(con, statement = paste("SELECT * FROM ILLUMINA_PICARD_METRICS"));
        d = fetch(rs, n=-1);
        dbHasCompleted(rs);
        dbClearResult(rs);
    } else {
        rs = dbSendQuery(con, statement = paste("SELECT * FROM ILLUMINA_SAMPLE_STATUS_AGG"));
        d = fetch(rs, n=-1);
        dbHasCompleted(rs);
        dbClearResult(rs);
    }

    oraCloseDriver(drv);

    subset(d, Project == project);
}

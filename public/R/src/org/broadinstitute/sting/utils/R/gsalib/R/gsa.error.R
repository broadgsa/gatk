gsa.error <- function(message) {
    message("");
    gsa.message("Error: **********");
    gsa.message(sprintf("Error: %s", message));
    gsa.message("Error: **********");
    message("");

    traceback();

    message("");
    stop(message, call. = FALSE);
}

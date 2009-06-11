package org.broadinstitute.sting.gatk.datasources.providers;

/**
 * Created by IntelliJ IDEA.
 * User: hanna
 * Date: Apr 16, 2009
 * Time: 4:11:40 PM
 *
 * Thrown to indicate invalid positions passed to the providers.
 * Extend from RuntimeException to make it easier on our walker writers; don't make
 * them catch every exception that comes their way.
 */
public class InvalidPositionException extends RuntimeException {
    public InvalidPositionException(String message) {
        super(message);
    }

    public InvalidPositionException(String message, Throwable throwable) {
        super(message,throwable);
    }
}

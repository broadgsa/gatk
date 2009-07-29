package org.broadinstitute.sting.utils.genotype;


/**
 * 
 * @author aaron 
 * 
 * Class GenotypeException
 *
 * This exception is thrown when a genotype call is passed in that cannot be processed, i.e. invalid.
 */
public class InvalidGenotypeException extends Exception {
    public InvalidGenotypeException(String msg) {
        super(msg);
    }

    public InvalidGenotypeException(String message, Throwable throwable) {
        super(message, throwable);
    }
}

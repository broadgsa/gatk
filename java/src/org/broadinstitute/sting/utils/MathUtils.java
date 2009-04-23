package org.broadinstitute.sting.utils;

/**
 * MathUtils is a static class (no instantiation allowed!) with some useful math methods.
 *
 * @author Kiran Garimella
 */
public class MathUtils {
    /** Private constructor.  No instantiating this class! */
    private MathUtils() {}

    /**
     * Compares double values for equality (within 1e-6), or inequality.
     *
     * @param a  the first double value
     * @param b  the second double value
     * @return   -1 if a is greater than b, 0 if a is equal to be within 1e-6, 1 if b is greater than a.
     */
    public static byte compareDoubles(double a, double b) { return compareDoubles(a, b, 1e-6); }

    /**
     * Compares double values for equality (within epsilon), or inequality.
     *
     * @param a       the first double value
     * @param b       the second double value
     * @param epsilon the precision within which two double values will be considered equal
     * @return        -1 if a is greater than b, 0 if a is equal to be within epsilon, 1 if b is greater than a.
     */
    public static byte compareDoubles(double a, double b, double epsilon)
    {
        if (Math.abs(a - b) < epsilon) { return 0; }
        if (a > b) { return -1; }
        return 1;
    }

    /**
     * Compares float values for equality (within 1e-6), or inequality.
     *
     * @param a  the first float value
     * @param b  the second float value
     * @return   -1 if a is greater than b, 0 if a is equal to be within 1e-6, 1 if b is greater than a.
     */
    public static byte compareFloats(float a, float b) { return compareFloats(a, b, 1e-6f); }

    /**
     * Compares float values for equality (within epsilon), or inequality.
     *
     * @param a       the first float value
     * @param b       the second float value
     * @param epsilon the precision within which two float values will be considered equal
     * @return        -1 if a is greater than b, 0 if a is equal to be within epsilon, 1 if b is greater than a.
     */
    public static byte compareFloats(float a, float b, float epsilon)
    {
        if (Math.abs(a - b) < epsilon) { return 0; }
        if (a > b) { return -1; }
        return 1;
    }
}

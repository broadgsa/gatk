package edu.mit.broad.picard.metrics;

import edu.mit.broad.picard.PicardException;
import edu.mit.broad.picard.util.FormatUtil;

import java.lang.reflect.Field;

/**
 * A base class from which all Metric classes should inherit.
 *
 * @author Tim Fennell
 */
public class MetricBase {
    /**
     * An equals method that checks equality by asserting that the classes are of the exact
     * same type and that all public fields are equal.
     *
     * @param o an instance to compare to
     * @return true if they are equal, false otherwise
     */
    public boolean equals(Object o) {
        if (o == null) return false;
        if (o.getClass() != getClass()) return false;

        // Loop through all the fields and check that they are either
        // null in both objects or equal in both objects
        for (Field f : getClass().getFields()) {
            try {
                Object lhs = f.get(this);
                Object rhs = f.get(o);

                if (lhs == null) {
                    if (rhs == null) {
                        // keep going
                    }
                    else if (rhs != null) {
                        return false;
                    }
                }
                else {
                    if (lhs.equals(rhs)) {
                        // keep going
                    }
                    else {
                        return false;
                    }
                }
            }
            catch (IllegalAccessException iae) {
                throw new PicardException("Could not read field " + f.getName() + " from a " + getClass().getSimpleName());
            }
        }

        // If we got this far all the fields are equal
        return true;
    }

    /** Converts the metric class to a human readable string. */
    public String toString() {
        StringBuilder buffer = new StringBuilder();
        FormatUtil formatter = new FormatUtil();

        for (Field f : getClass().getFields()) {
            try {
                buffer.append(f.getName());
                buffer.append("\t");
                buffer.append(formatter.format(f.get(this)));
                buffer.append("\n");
            }
            catch (IllegalAccessException iae) {
                throw new PicardException("Could not read field " + f.getName() + " from a " + getClass().getSimpleName());
            }
        }

        return buffer.toString();
    }
}

package org.broadinstitute.sting.utils.windowmaker;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.Pair;

import java.util.Iterator;

/**
 * Created by IntelliJ IDEA.
 * User: depristo
 * Date: Apr 14, 2009
 * Time: 9:33:00 AM
 * To change this template use File | Settings | File Templates.
 */
public interface PositionalDataGenerator<Pos, T> {
    public Pos lastPos();
}
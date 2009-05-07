package org.broadinstitute.sting.utils.windowmaker;

import org.broadinstitute.sting.utils.Pair;

import java.util.ArrayList;

/**
 * System for allowing "windowed" access into a data stream
 *
 * User: depristo
 * Date: May 1, 2009
 * Time: 3:03:20 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class WindowMaker<Pos, T> {
    private PositionalDataGenerator<Pos, T> dataSource;

    public WindowMaker(PositionalDataGenerator pdg)
    {
        dataSource = pdg;        
    }

    /**
     * 
     */
    public abstract boolean hasNext();

    /**
     * Get the next object in the stream. Does not advance the stream pointer.  Successive peek()
     * calls return the same object.
     * 
     * @return
     */
    public abstract Pair<Pos, T> peek();

    /**
     * Return the next N objects in the stream.  Does not advance the stream pointer.  Successive peek()
     * calls return the same Array of Objects.
     *
     * @param N
     * @return
     */
    public abstract ArrayList<Pair<Pos, T>> peek(int N);

    /**
     * Pop the leftmost object off the window
     */
    public abstract void pop();

    /**
     * Pop the leftmost N objects off the window
     * @param N
     */
    public abstract void pop(int N);

    /**
     * Advance the internal stream pointer to from, collect all data until to, and return the data
     * from -> to as an ArrayList
     * 
     * @param from
     * @param to
     * @return
     */
    public abstract ArrayList<Pair<Pos, T>> speek(Pos from, Pos to);

    /**
     * Advance the internal stream pointer to from, collect N data units, and return the data
     * as an ArrayList
     *
     * @param from
     * @param length
     * @return
     */
    public abstract ArrayList<Pair<Pos, T>> speek(Pos from, int N);


    /**
     * Advance the internal stream pointer to from, collect all remaining data units from the
     * data stream, and return them all as one (potentially large) arraylist
     *
     * @param from
     * @return
     */
    public abstract ArrayList<Pair<Pos, T>> speek(Pos from);
}

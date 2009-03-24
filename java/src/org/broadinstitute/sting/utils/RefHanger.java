package org.broadinstitute.sting.utils;

import org.broadinstitute.sting.utils.GenomeLoc;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Hanging data off the reference sequence
 *
 * Supports in effect the following data structure
 *
 * <-- reference bases: A T G C -->
 *                      d d d d
 *                      d d d d
 *                      d d d d
 *                      d d   d
 *                        d   d
 *                        d
 *
 * Where the little d's are data associated with each position in the reference.
 * Supports adding and removing data to either side of the data structure, as well as
 * randomly accessing data anywhere within window.  
 */
public class RefHanger<T> {

    // -----------------------------------------------------------------------------------------------------------------
    //
    // member fields
    //
    // -----------------------------------------------------------------------------------------------------------------
    ArrayList<Hanger> hangers;

    // -----------------------------------------------------------------------------------------------------------------
    //
    // Info structure
    //
    // -----------------------------------------------------------------------------------------------------------------
    public class Hanger {
        public GenomeLoc loc = null;
        public ArrayList<T> data = null;
        
        public Hanger(GenomeLoc loc, ArrayList<T> data) {
            this.loc = loc;
            this.data = data;
        }

        public final ArrayList<T> getData() { return data; }
        public final int size() { return this.data.size(); }
        public final T get(int i) { return this.data.get(i); }
    }


    // -----------------------------------------------------------------------------------------------------------------
    //
    // constructors and other basic operations
    //
    // -----------------------------------------------------------------------------------------------------------------
    public RefHanger() {
        hangers = new ArrayList<Hanger>();
    }

    public void clear() {
        hangers.clear();
        //System.out.printf("leftLoc is %s%n", getLeftLoc());
    }

    protected int getLeftOffset() { return 0; }
    protected int getRightOffset() { return hangers.size() - 1; }
    protected int getOffset(GenomeLoc loc) {
        //System.out.printf("Loc: %s vs %s%n", loc, getLeftLoc());
        return loc.minus(getLeftLoc());
    }
    
    public GenomeLoc getLeftLoc() { return hangers.get(getLeftOffset()).loc; }
    public GenomeLoc getRightLoc() { return hangers.get(getRightOffset()).loc; }

    public boolean hasLocation(GenomeLoc loc) {
        return ! isEmpty() && loc.isBetween(getLeftLoc(), getRightLoc());
    }

    public boolean isEmpty() {
        return hangers.isEmpty();
    }
    public boolean hasHangers() {
        return ! isEmpty();
    }

    /**
     * Pops off the left most data from the structure
     *
     * @return
     */
    public Hanger popLeft() {
        assert hasHangers();
        return hangers.remove(0);
    }

    public void dropLeft() {
        popLeft();
    }

    /**
     * Looks at the left most data from the structure
     *
     * @return
     */
    public Hanger getLeft() {
        assert hasHangers();
        return getHanger(0);
    }

    /**
     * Returns data at offset relativePos
     *
     * @return
     */
    public Hanger getHanger(int relativePos) {
        //assert hangers.contains(relativePos) : hangers + " " + relativePos;
        return hangers.get(relativePos);
    }

    /**
     * Returns data at GenomeLoc
     *
     * @return
     */
    public Hanger getHanger(GenomeLoc pos) {
        return getHanger(getOffset(pos));
    }

    public int size() {
        return hangers.size();
    }

    // -----------------------------------------------------------------------------------------------------------------
    //
    // Adding data to the left and right
    //
    // -----------------------------------------------------------------------------------------------------------------

    public void pushLeft(GenomeLoc pos) {
        pushLeft(pos, new ArrayList<T>());
    }

    public void pushLeft(GenomeLoc pos, T datum) {
        pushLeft(pos, new ArrayList<T>(Arrays.asList(datum)));
    }

    public void pushLeft(GenomeLoc pos, ArrayList<T> data) {
        hangers.add(0, new Hanger(pos, data));
    }
    
    public void pushRight(GenomeLoc pos) {
        pushRight(pos, new ArrayList<T>());
    }

    public void pushRight(GenomeLoc pos, T datum) {
        pushRight(pos, new ArrayList<T>(Arrays.asList(datum)));
    }

    public void pushRight(GenomeLoc pos, ArrayList<T> data) {
        hangers.add(new Hanger(pos, data));
    }

    public boolean ensurePos(GenomeLoc pos) {
        if ( hasLocation(pos) )
            return true;
        else {
            pushRight(pos);
            return false;
        }
    }

    public void extendData(GenomeLoc pos, T datum) {
        getHanger(pos).data.add(datum);
    }

    public void addData(List<GenomeLoc> positions, List<T> dataByPos) {
        assert( positions.size() == dataByPos.size() );

        for ( int i = 0; i < positions.size(); i++ ) {
            GenomeLoc pos = positions.get(i);
            T datum = dataByPos.get(i);
            expandingPut1(pos, datum);
        }
    }

    public void expandingPut1(final GenomeLoc loc, T datum) {
        ensurePos(loc);
        extendData(loc, datum);
    }

    public void printState() {
        System.out.printf("Hanger:%n");
        for ( Hanger hanger : hangers ) {
            System.out.printf("  -> %s => %s:%n", hanger.loc, Utils.join("/", hanger.data) );
       }
    }

    /**
     * Pushes locations on the right until we reach the expected position for pos.
     *
     * For example, if we have chr1:1 and 2 in the hanger, and we push 4 into the hangers
     * this function will add 3 -> {} to the hanger too
     *
     * @param pos
     * @param datum
     */
    public void expandingPut(GenomeLoc pos, T datum) {
        //System.out.printf("expandingPut(%s, %s)%n", pos, datum);
        //printState();
        if ( isEmpty() )
            // we have nothing, just push right
            pushRight(pos, datum);
        else {
            //assert pos.compareTo(getRightLoc()) == 1 : pos + " " + getRightLoc() + " => " + pos.compareTo(getRightLoc());

            GenomeLoc nextRight = getRightLoc().nextLoc();
            while ( pos.compareTo(nextRight) == 1 ) {
                //printState();
                //System.out.printf("    *** Extending %s, heading for %s%n", nextRight, pos);
                ensurePos(nextRight);
                nextRight = nextRight.nextLoc();
            }

            ensurePos(pos);
            extendData(pos, datum);
        }
        //printState();
    }
}
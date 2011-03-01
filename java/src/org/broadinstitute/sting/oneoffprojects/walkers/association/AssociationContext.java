package org.broadinstitute.sting.oneoffprojects.walkers.association;

import java.lang.reflect.InvocationTargetException;
import java.util.*;

/**
 * @Author chartl
 * @Date 2011-02-23
 * A general context for windowed association
 */
public abstract class AssociationContext<X extends AssociationContextAtom> {
    protected Class<? extends AssociationContextAtom> clazz;
    protected List<X> window;

    public AssociationContext( Class<X> zclaz ) {
        window = new ArrayList<X>(getWindowSize());
        clazz = zclaz;
    }

    public X map(MapExtender e) throws NoSuchMethodException, InstantiationException, IllegalAccessException, InvocationTargetException, ClassCastException {
        return (X) clazz.getConstructor(new Class[] {MapExtender.class}).newInstance(new Object[] {e});
    }

    public boolean filter(MapExtender m) { return true; }

    public void reduce(X context) {
        window.add(context);
    }

    public abstract int getWindowSize();
    public abstract int slideByValue();

    public boolean isFull() {
        return window.size() >= getWindowSize();
    }

    public void slide() {
        ArrayList<X> newWindow = new ArrayList<X>((window.subList(slideByValue(),window.size())));
        newWindow.ensureCapacity(getWindowSize());
        window = newWindow;
    }
}
/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.utils.instrumentation;

import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;

import java.lang.instrument.Instrumentation;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.lang.reflect.Modifier;
import java.util.IdentityHashMap;

/**
 * A sizeof implementation for Java.  Relies on the Java instrumentation API, so
 * it must be added as an agent to function properly.
 *
 * To run, add -javaagent:$STING_HOME/dist/StingUtils.jar as a command-line
 * JVM argument.
 *
 * @author mhanna
 * @version 0.1
 */
public class Sizeof {
    /**
     * Instrumentation object.  Registered by the JVM via the premain() method.
     */
    private static Instrumentation instrumentation;

    /**
     * Called by the JVM before the agent is started.
     * @param args Arguments?
     * @param inst Instrumentation object, used to perform instrumentation in the JVM.
     */
    public static void premain(String args, Instrumentation inst) {
        instrumentation = inst;
    }

    /**
     * Is this Sizeof operator enabled?  To enable, add the -javaagent directive listed in the class-level javadoc.
     * @return True if sizeof() is enabled.  If false, any calls to utility methods of this class will throw an exception.
     */
    public static boolean isEnabled() {
        return instrumentation != null;
    }

    /**
     * Gets the size of the given object.  Retrieves the size for only this object; any reference fields in the object will only be
     * counted as single pointers.
     * @param o The object to sizeof().
     * @return Gets the best possible approximation we can get of the size of the object in memory.  On Sun JVM, includes some object padding.
     */
    public static long getObjectSize(Object o) {
        if(!isEnabled())
            throw new ReviewedGATKException("Sizeof operator is currently disabled!  To enable, review the documentation in Sizeof.java");
        return instrumentation.getObjectSize(o);
    }

    /**
     * Gets the size of the given object, including the size of the objects to which this object refers.
     * @param o The object to sizeof().
     * @return Gets the best possible approximation we can get of the size of the object in memory, including all references within each object.
     */
    public static long getObjectGraphSize(Object o) {
        if(!isEnabled())
            throw new ReviewedGATKException("Sizeof operator is currently disabled!  To enable, review the documentation in Sizeof.java");
        IdentityHashMap<Object,Object> objectsSeen = new IdentityHashMap<Object,Object>();
        return getObjectGraphSize(o,objectsSeen);
    }

    /**
     * The engine for walking the graph of all objects and their children.
     * @param o The object to traverse.
     * @param objectsSeen A list of all objects already seen.
     * @return Gets the best possible approximation we can get of the size of the object in memory, including all references within each object.
     */
    private static long getObjectGraphSize(Object o,IdentityHashMap<Object,Object> objectsSeen) {
        // Size of a null object itself (as opposed to the reference to the null object) is 0.
        if(o == null)
            return 0;
        
        // Don't allow repeated traversals of the same object.
        if(objectsSeen.containsKey(o))
            return 0;
        objectsSeen.put(o,o);

        // Get the size of the object itself, plus all contained primitives.
        long totalSize = instrumentation.getObjectSize(o);

        // Get the size of (non-primitive) array elements.
        Class<?> classToInspect = o.getClass();
        if(classToInspect.isArray()) {
            if(!classToInspect.getComponentType().isPrimitive()) {
                for(int i = 0; i < Array.getLength(o); i++)
                    totalSize += getObjectGraphSize(Array.get(o,i),objectsSeen);
            }
        }

        // Walk the descendents of each field of this class.  Be sure to avoid synthetic fields like this$0 -- these
        // are back references to the parent of the object contained in the inner class.
        // Potential BUG: Are there other types of synthetic fields we should be tracking?
        while(classToInspect != null) {
            for(Field field: classToInspect.getDeclaredFields()) {
                if(field.getType().isPrimitive())
                    continue;
                if(Modifier.isStatic(field.getModifiers()))
                    continue;
                if(field.isSynthetic())
                    continue;
                field.setAccessible(true);
                Object fieldValue;
                try {
                    fieldValue = field.get(o);
                }
                catch(IllegalAccessException ex) {
                    throw new ReviewedGATKException("Unable to access field " + field.getName(),ex);
                }
                totalSize += getObjectGraphSize(fieldValue,objectsSeen);
            }
            classToInspect = classToInspect.getSuperclass();
        }
        return totalSize;
    }
}

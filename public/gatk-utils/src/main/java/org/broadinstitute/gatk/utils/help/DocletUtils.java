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

package org.broadinstitute.gatk.utils.help;

import com.sun.javadoc.FieldDoc;
import com.sun.javadoc.PackageDoc;
import com.sun.javadoc.ProgramElementDoc;
import org.broadinstitute.gatk.utils.classloader.JVMUtils;

import java.lang.reflect.Field;

/**
 * Methods in the class must ONLY be used by doclets, since the com.sun.javadoc.* classes are not
 * available on all systems, and we don't want the GATK proper to depend on them.
 */
public class DocletUtils {

    protected static boolean assignableToClass(ProgramElementDoc classDoc, Class lhsClass, boolean requireConcrete) {
        try {
            Class type = getClassForDoc(classDoc);
            return lhsClass.isAssignableFrom(type) && (!requireConcrete || JVMUtils.isConcrete(type));
        } catch (Throwable t) {
            // Ignore errors.
            return false;
        }
    }

    protected static Class getClassForDoc(ProgramElementDoc doc) throws ClassNotFoundException {
        return Class.forName(getClassName(doc, true));
    }

    protected static Field getFieldForFieldDoc(FieldDoc fieldDoc) {
        try {
            Class clazz = getClassForDoc(fieldDoc.containingClass());
            return JVMUtils.findField(clazz, fieldDoc.name());
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     * Reconstitute the class name from the given class JavaDoc object.
     *
     * @param doc the Javadoc model for the given class.
     * @return The (string) class name of the given class.
     */
    protected static String getClassName(ProgramElementDoc doc, boolean binaryName) {
        PackageDoc containingPackage = doc.containingPackage();
        String className = doc.name();
        if (binaryName) {
            className = className.replaceAll("\\.", "\\$");
        }
        return containingPackage.name().length() > 0 ?
                String.format("%s.%s", containingPackage.name(), className) :
                String.format("%s", className);
    }
}
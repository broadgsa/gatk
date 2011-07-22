package org.broadinstitute.sting.utils.help;

import com.sun.javadoc.FieldDoc;
import com.sun.javadoc.PackageDoc;
import com.sun.javadoc.ProgramElementDoc;

import java.lang.reflect.Field;

public class HelpUtils {
    protected static boolean implementsInterface(ProgramElementDoc classDoc, Class... interfaceClasses) {
        for (Class interfaceClass : interfaceClasses)
            if (assignableToClass(classDoc, interfaceClass, false))
                return true;
        return false;
    }

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
        return Class.forName(getClassName(doc));
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
    protected static String getClassName(ProgramElementDoc doc) {
        PackageDoc containingPackage = doc.containingPackage();
        return containingPackage.name().length() > 0 ?
                String.format("%s.%s", containingPackage.name(), doc.name()) :
                String.format("%s", doc.name());
    }
}
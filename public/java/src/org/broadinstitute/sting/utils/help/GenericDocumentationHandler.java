/*
 * Copyright (c) 2011, The Broad Institute
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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

package org.broadinstitute.sting.utils.help;

import com.google.java.contract.Requires;
import com.sun.javadoc.ClassDoc;
import com.sun.javadoc.FieldDoc;
import com.sun.javadoc.RootDoc;
import com.sun.javadoc.Tag;
import org.apache.log4j.Logger;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.CommandLineGATK;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.classloader.JVMUtils;
import org.broadinstitute.sting.utils.exceptions.ReviewedStingException;

import java.io.*;
import java.lang.reflect.Field;
import java.util.*;

/**
 *
 */
public class GenericDocumentationHandler extends DocumentedGATKFeatureHandler {
    private static Logger logger = Logger.getLogger(GenericDocumentationHandler.class);
    GATKDocWorkUnit toProcess;
    ClassDoc classdoc;
    Set<GATKDocWorkUnit> all;
    RootDoc rootDoc;

    @Override
    public boolean includeInDocs(ClassDoc doc) {
        return true;
//        try {
//            Class type = HelpUtils.getClassForDoc(doc);
//            return JVMUtils.isConcrete(type);
//        } catch ( ClassNotFoundException e ) {
//            return false;
//        }
    }


    @Override
    public String getTemplateName(ClassDoc doc) throws IOException {
        return "generic.template.html";
    }

    @Override
    public void processOne(RootDoc rootDoc, GATKDocWorkUnit toProcessArg, Set<GATKDocWorkUnit> allArg) {
        this.rootDoc = rootDoc;
        this.toProcess = toProcessArg;
        this.all = allArg;
        this.classdoc = toProcess.classDoc;

        //System.out.printf("%s class %s%n", toProcess.group, toProcess.classDoc);
        Map<String, Object> root = new HashMap<String, Object>();

        addHighLevelBindings(root);
        addArgumentBindings(root);
        addRelatedBindings(root);

        toProcess.setHandlerContent((String)root.get("summary"), root);
    }

    protected void addHighLevelBindings(Map<String, Object> root) {
        root.put("name", classdoc.name());

        // Extract overrides from the doc tags.
        StringBuilder summaryBuilder = new StringBuilder();
        for(Tag tag: classdoc.firstSentenceTags())
            summaryBuilder.append(tag.text());
        root.put("summary", summaryBuilder.toString());
        root.put("description", classdoc.commentText().substring(summaryBuilder.toString().length()));
        root.put("timestamp", toProcess.buildTimestamp);
        root.put("version", toProcess.absoluteVersion);

        for(Tag tag: classdoc.tags()) {
            root.put(tag.name(), tag.text());
        }
    }

    protected void addArgumentBindings(Map<String, Object> root) {
        ParsingEngine parsingEngine = createStandardGATKParsingEngine();

        // attempt to instantiate the class
        Object instance = makeInstanceIfPossible(toProcess.clazz);

        Map<String, List<Object>> args = new HashMap<String, List<Object>>();
        root.put("arguments", args);
        args.put("all", new ArrayList<Object>());
        args.put("required", new ArrayList<Object>());
        args.put("optional", new ArrayList<Object>());
        args.put("hidden", new ArrayList<Object>());
        args.put("depreciated", new ArrayList<Object>());
        try {
            for ( ArgumentSource argumentSource : parsingEngine.extractArgumentSources(HelpUtils.getClassForDoc(classdoc)) ) {
                ArgumentDefinition argDef = argumentSource.createArgumentDefinitions().get(0);
                FieldDoc fieldDoc = getFieldDoc(classdoc, argumentSource.field.getName());
                Map<String, Object> argBindings = docForArgument(fieldDoc, argumentSource, argDef); // todo -- why can you have multiple ones?
                if ( ! argumentSource.isHidden() || getDoclet().showHiddenFeatures() ) {
                    logger.debug(String.format("Processing %s", argumentSource));
                    String kind = "optional";
                    if ( argumentSource.isRequired() ) kind = "required";
                    else if ( argumentSource.isHidden() ) kind = "hidden";
                    else if ( argumentSource.isDeprecated() ) kind = "depreciated";

                    // get the value of the field
                    if ( instance != null ) {
                        Object value = getFieldValue(toProcess.clazz, instance, fieldDoc.name());
                        if ( value != null )
                            argBindings.put("defaultValue", prettyPrintValueString(value));
                    }

                    args.get(kind).add(argBindings);
                    args.get("all").add(argBindings);
                } else {
                    logger.debug(String.format("Skipping hidden feature %s", argumentSource));
                }
            }
        } catch ( ClassNotFoundException e ) {
            throw new RuntimeException(e);
        }
    }

    private Object getFieldValue(Class c, Object instance, String fieldName) {
        Field field = JVMUtils.findField(c, fieldName);
        if ( field != null ) {
            Object value = JVMUtils.getFieldValue(field, instance);
            //System.out.printf("Fetched value of field %s in class %s: %s%n", fieldName, c, value);
            return value;
        } else {
            return findFieldValueInArgumentCollections(c, instance, fieldName);
        }
    }

    private Object findFieldValueInArgumentCollections(Class c, Object instance, String fieldName) {
        for ( Field field : JVMUtils.getAllFields(c) ) {
            if ( field.isAnnotationPresent(ArgumentCollection.class) ) {
                //System.out.printf("Searching for %s in argument collection field %s%n", fieldName, field);
                Object fieldValue = JVMUtils.getFieldValue(field, instance);
                Object value = getFieldValue(fieldValue.getClass(), fieldValue, fieldName);
                if ( value != null )
                    return value;
            }
        }

        return null;
    }

    /**
     * Assumes value != null
     * @param value
     * @return
     */
    private Object prettyPrintValueString(Object value) {
        if ( value.getClass().isArray() ) {
            Class type = value.getClass().getComponentType();
            if ( boolean.class.isAssignableFrom(type) )
                return Arrays.toString((boolean[])value);
            if ( byte.class.isAssignableFrom(type) )
                return Arrays.toString((byte[])value);
            if ( char.class.isAssignableFrom(type) )
                return Arrays.toString((char[])value);
            if ( double.class.isAssignableFrom(type) )
                return Arrays.toString((double[])value);
            if ( float.class.isAssignableFrom(type) )
                return Arrays.toString((float[])value);
            if ( int.class.isAssignableFrom(type) )
                return Arrays.toString((int[])value);
            if ( long.class.isAssignableFrom(type) )
                return Arrays.toString((long[])value);
            if ( short.class.isAssignableFrom(type) )
                return Arrays.toString((short[])value);
            if ( Object.class.isAssignableFrom(type) )
                return Arrays.toString((Object[])value);
            else
                throw new RuntimeException("Unexpected array type in prettyPrintValue.  Value was " + value + " type is " + type);
        } else
            return value.toString();
    }

    private Object makeInstanceIfPossible(Class c) {
        Object instance = null;
        try {
            // don't try to make something where we will obviously fail
            if (! c.isEnum() && ! c.isAnnotation() && ! c.isAnonymousClass() &&
                    ! c.isArray() && ! c.isPrimitive() & JVMUtils.isConcrete(c) ) {
                instance = c.newInstance();
                //System.out.printf("Created object of class %s => %s%n", c, instance);
                return instance;
            } else
                return null;
        }
        catch (IllegalAccessException e ) { }
        catch (InstantiationException e ) { }
        catch (ExceptionInInitializerError e ) { }
        catch (SecurityException e ) { }
        // this last one is super dangerous, but some of these methods catch ClassNotFoundExceptions
        // and rethrow then as RuntimeExceptions
        catch (RuntimeException e) {}
        finally {
            if ( instance == null )
                logger.warn(String.format("Unable to create instance of class %s => %s", c, instance));
        }

        return instance;
    }

    protected void addRelatedBindings(Map<String, Object> root) {
        List<Map<String, Object>> extraDocsData = new ArrayList<Map<String, Object>>();

        // add in all of the explicitly related items
        for ( final Class extraDocClass : toProcess.annotation.extraDocs() ) {
            final GATKDocWorkUnit otherUnit = GATKDoclet.findWorkUnitForClass(extraDocClass, all);
            if ( otherUnit == null )
                throw new ReviewedStingException("Requested extraDocs for class without any documentation: " + extraDocClass);
            extraDocsData.add(
                    new HashMap<String, Object>(){{
                        put("filename", otherUnit.filename);
                        put("name", otherUnit.name);}});

        }

        List<Map<String, Object>> hierarchyDocs = new ArrayList<Map<String, Object>>();
        for (final GATKDocWorkUnit other : all ) {
            final String relation = classRelationship(toProcess.clazz, other.clazz);
            if ( relation != null )
                hierarchyDocs.add(
                        new HashMap<String, Object>(){{
                            put("filename", other.filename);
                            put("relation", relation);
                            put("name", other.name);}});

        }

        root.put("relatedDocs", hierarchyDocs);
        root.put("extradocs", extraDocsData);
    }

    private static final String classRelationship(Class me, Class other) {
        if ( other.equals(me) )
            // no circular references
            return null;
        else if ( other.isAssignableFrom(me) )
            // toProcess is a superclass of other.clazz
            return "superclass";
        else if ( me.isAssignableFrom(other) )
            // toProcess inherits from other.clazz
            return "subclass";
        else
            return null;

    }

    protected ParsingEngine createStandardGATKParsingEngine() {
        CommandLineProgram clp = new CommandLineGATK();
        try {
            CommandLineProgram.start(clp, new String[]{}, true);
            return clp.parser;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    private FieldDoc getFieldDoc(ClassDoc classDoc, String name) {
        return getFieldDoc(classDoc, name, true);
    }

    private FieldDoc getFieldDoc(ClassDoc classDoc, String name, boolean primary) {
        //System.out.printf("Looking for %s in %s%n", name, classDoc.name());
        for ( FieldDoc fieldDoc : classDoc.fields(false) ) {
            //System.out.printf("fieldDoc " + fieldDoc + " name " + fieldDoc.name());
            if ( fieldDoc.name().equals(name) )
                return fieldDoc;

            Field field = HelpUtils.getFieldForFieldDoc(fieldDoc);
            if ( field.isAnnotationPresent(ArgumentCollection.class) ) {
                ClassDoc typeDoc = getRootDoc().classNamed(fieldDoc.type().qualifiedTypeName());
                if ( typeDoc == null )
                    throw new ReviewedStingException("Tried to get javadocs for ArgumentCollection field " + fieldDoc + " but could't find the class in the RootDoc");
                else {
                    FieldDoc result = getFieldDoc(typeDoc, name, false);
                    if ( result != null )
                        return result;
                    // else keep searching
                }
            }
        }

        // if we didn't find it here, wander up to the superclass to find the field
        if ( classDoc.superclass() != null ) {
            return getFieldDoc(classDoc.superclass(), name, false);
        }

        if ( primary )
            throw new RuntimeException("No field found for expected field " + name);
        else
            return null;
    }

    protected Map<String, Object> docForArgument(FieldDoc fieldDoc, ArgumentSource source, ArgumentDefinition def) {
        Map<String, Object> root = new HashMap<String, Object>();
        root.put("name", def.shortName != null ? "-" + def.shortName : "--" + def.fullName );

        if ( def.shortName != null && def.fullName != null )
            root.put("synonyms", "--" + def.fullName);

        root.put("required", def.required ? "yes" : "no");
        root.put("type", def.argumentType.getSimpleName());

        // summary and fulltext
        root.put("summary", def.doc != null ? def.doc : "");
        root.put("fulltext", fieldDoc.commentText());

        List<String> attributes = new ArrayList<String>();
        // this one below is just too much.
        //attributes.add(def.ioType.annotationClass.getSimpleName());
        if ( def.required ) attributes.add("required");
        // flag is just boolean, not interesting
        //if ( def.isFlag ) attributes.add("flag");
        if ( def.isHidden ) attributes.add("hidden");
        if ( source.isDeprecated() ) attributes.add("depreciated");
        if ( attributes.size() > 0 )
            root.put("attributes", Utils.join(", ", attributes));

        if ( def.validOptions != null ) {
            root.put("options", docForEnumArgument(source.field.getType()));
        }

        return root;
    }

    @Requires("enumClass.isEnum()")
    private List<Map<String, Object>> docForEnumArgument(Class enumClass) {
        ClassDoc doc = GATKDoclet.getClassDocForClass(rootDoc, enumClass);
        if ( doc == null ) //  || ! doc.isEnum() )
            throw new RuntimeException("Tried to get docs for enum " + enumClass + " but got instead: " + doc);

        List<Map<String, Object>> bindings = new ArrayList<Map<String, Object>>();
        for (final FieldDoc field : doc.fields(false) ) {
            bindings.add(
                    new HashMap<String, Object>(){{
                        put("name", field.name());
                        put("summary", field.commentText());}});
        }

        return bindings;
    }

}

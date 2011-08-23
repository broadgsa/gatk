/*
 * Copyright (c) 2010, The Broad Institute
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

package org.broadinstitute.sting.queue.extensions.gatk;

import net.sf.samtools.BAMIndex;
import net.sf.samtools.SAMFileWriter;
import org.broad.tribble.Tribble;
import org.broadinstitute.sting.commandline.*;
import org.broadinstitute.sting.gatk.io.stubs.SAMFileWriterArgumentTypeDescriptor;
import org.broadinstitute.sting.utils.codecs.vcf.VCFWriter;

import java.io.File;
import java.lang.annotation.Annotation;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public abstract class ArgumentDefinitionField extends ArgumentField {

    protected final ArgumentDefinition argumentDefinition;
    protected ArgumentDefinitionField(ArgumentDefinition argumentDefinition) {
        this.argumentDefinition = argumentDefinition;
    }

    @Override protected String getRawFieldName() { return argumentDefinition.fullName; }
    @Override protected Class<? extends Annotation> getAnnotationIOClass() { return argumentDefinition.ioType.annotationClass; }
    @Override protected String getDoc() { return escape(argumentDefinition.doc); }
    @Override protected String getFullName() { return escape(argumentDefinition.fullName); }
    @Override protected String getShortName() { return escape(argumentDefinition.shortName); }
    @Override protected boolean isRequired() { return argumentDefinition.required; }
    @Override protected String getExclusiveOf() { return escape(argumentDefinition.exclusiveOf); }
    @Override protected String getValidation() { return escape(argumentDefinition.validation); }
    protected boolean isFlag() { return argumentDefinition.isFlag; }
    protected boolean isMultiValued() { return argumentDefinition.isMultiValued; }

    protected final String getShortFieldGetter() { return getFieldName(getRawShortFieldName()); }
    protected final String getShortFieldSetter() { return getFieldName(getRawShortFieldName() + "_="); }
    protected String getRawShortFieldName() { return argumentDefinition.shortName; }
    @Override protected String getDefineAddition() {
        if (argumentDefinition.shortName == null)
            return "";
        else if(getShortFieldGetter().equals(getFieldName()))
            return "";
        else
            return String.format("%n" +
                    "/**%n" +
                    " * Short name of %1$s%n" +
                    " * @return Short name of %1$s%n" +
                    " */%n" +
                    "def %3$s = this.%1$s%n" +
                    "%n" +
                    "/**%n" +
                    " * Short name of %1$s%n" +
                    " * @param value Short name of %1$s%n" +
                    " */%n" +
                    "def %4$s(value: %2$s) { this.%1$s = value }%n",
                    getFieldName(),
                    getFieldType(),
                    getShortFieldGetter(),
                    getShortFieldSetter());
    }

    protected static final String REQUIRED_TEMPLATE = " + \" %1$s \" + %2$s.format(%3$s)";
    protected static final String REPEAT_TEMPLATE = " + repeat(\" %1$s \", %3$s, format=formatValue(%2$s))";
    protected static final String OPTIONAL_TEMPLATE = " + optional(\" %1$s \", %3$s, format=formatValue(%2$s))";
    protected static final String FLAG_TEMPLATE = " + (if (%3$s) \" %1$s\" else \"\")";

    public final String getCommandLineAddition() {
        return String.format(getCommandLineTemplate(), getCommandLineParam(), getCommandLineFormat(), getFieldName());
    }

    protected String getCommandLineParam() {
        return (argumentDefinition.shortName != null)
                ? "-" + argumentDefinition.shortName
                : "--" + argumentDefinition.fullName;
    }

    protected String getCommandLineFormat() {
        return "\"%s\"";
    }

    @Override
    protected String getGatherAnnotation() {
        return "";
    }

    protected String getCommandLineTemplate() {
        if (isFlag()) return FLAG_TEMPLATE;
        if (isMultiValued()) return REPEAT_TEMPLATE;
        if (isRequired()) return REQUIRED_TEMPLATE;
        return OPTIONAL_TEMPLATE;
    }

    public static List<? extends ArgumentField> getArgumentFields(ParsingEngine parsingEngine,Class<?> classType) {
        List<ArgumentField> argumentFields = new ArrayList<ArgumentField>();
        for (ArgumentSource argumentSource: parsingEngine.extractArgumentSources(classType))
            if (!argumentSource.isDeprecated()) {
                Class<?> gatherer = null;
                if (argumentSource.field.isAnnotationPresent(Gather.class))
                    gatherer = argumentSource.field.getAnnotation(Gather.class).value();
                for (ArgumentDefinition argumentDefinition: argumentSource.createArgumentDefinitions())
                    argumentFields.addAll(getArgumentFields(argumentDefinition, gatherer));
            }
        return argumentFields;
    }

    private static final List<String> intervalFields = Arrays.asList("intervals", "excludeIntervals", "targetIntervals");

    private static List<? extends ArgumentField> getArgumentFields(ArgumentDefinition argumentDefinition, Class<?> gatherer) {
        if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(
                    new IntervalFileArgumentField(argumentDefinition),
                    new IntervalStringArgumentField(argumentDefinition));

        // ROD Bindings are set by the RodBindField
        } else if (RodBindField.ROD_BIND_FIELD.equals(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            // TODO: Once everyone is using @Allows and @Requires correctly, we can stop blindly allowing Triplets
            return Arrays.asList(new RodBindArgumentField(argumentDefinition), new InputIndexesArgumentField(argumentDefinition, Tribble.STANDARD_INDEX_EXTENSION));
            //return Collections.<ArgumentField>emptyList();

        } else if ("input_file".equals(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(new InputTaggedFileDefinitionField(argumentDefinition), new InputIndexesArgumentField(argumentDefinition, BAMIndex.BAMIndexSuffix, ".bam"));

        } else if ((RodBinding.class.equals(argumentDefinition.argumentType) || RodBinding.class.equals(argumentDefinition.componentType)) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(new InputTaggedFileDefinitionField(argumentDefinition), new InputIndexesArgumentField(argumentDefinition, Tribble.STANDARD_INDEX_EXTENSION));

        } else if (argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Collections.singletonList(new InputArgumentField(argumentDefinition));

        } else if (argumentDefinition.ioType == ArgumentIOType.OUTPUT) {

            List<ArgumentField> fields = new ArrayList<ArgumentField>();

            String gatherClass;
            if (gatherer != null)
                gatherClass = gatherer.getName();
            else if (SAMFileWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                gatherClass = "BamGatherFunction";
            else if (VCFWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                gatherClass = "VcfGatherFunction";
            else
                gatherClass = "org.broadinstitute.sting.queue.function.scattergather.SimpleTextGatherFunction";

            fields.add(new OutputArgumentField(argumentDefinition, gatherClass));

            if (SAMFileWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                fields.add(new SAMFileWriterIndexArgumentField(argumentDefinition));
            else if (VCFWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                fields.add(new VCFWriterIndexArgumentField(argumentDefinition));

            return fields;

        } else if (argumentDefinition.isFlag) {
            return Collections.singletonList(new FlagArgumentField(argumentDefinition));

        } else if (argumentDefinition.isMultiValued) {
            return Collections.singletonList(new MultiValuedArgumentField(argumentDefinition));

        } else if (!argumentDefinition.required && useOption(argumentDefinition.argumentType)) {
            boolean useFormat = useFormatter(argumentDefinition.argumentType);
            List<ArgumentField> fields = new ArrayList<ArgumentField>();
            ArgumentField field = new OptionedArgumentField(argumentDefinition, useFormat);
            fields.add(field);
            if (useFormat) fields.add(new FormatterArgumentField(field));
            return fields;

        } else {
            boolean useFormat = useFormatter(argumentDefinition.argumentType);
            List<ArgumentField> fields = new ArrayList<ArgumentField>();
            ArgumentField field = new DefaultArgumentField(argumentDefinition, useFormat);
            fields.add(field);
            if (useFormat) fields.add(new FormatterArgumentField(field));
            return fields;

        }
    }

    // if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT)
    // Change intervals exclusive of intervalsString.
    private static class IntervalFileArgumentField extends InputArgumentField {
        public IntervalFileArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override
        protected String getExclusiveOf() {
            StringBuilder exclusiveOf = new StringBuilder(super.getExclusiveOf());
            if (exclusiveOf.length() > 0)
                exclusiveOf.append(",");
            exclusiveOf.append(escape(argumentDefinition.fullName)).append("String");
            return exclusiveOf.toString();
        }
    }

    // if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT)
    // Change intervals to a string but as an argument.
    private static class IntervalStringArgumentField extends ArgumentDefinitionField {
        public IntervalStringArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @SuppressWarnings("unchecked")
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Argument.class; }
        @Override protected Class<?> getInnerType() { return String.class; }
        @Override protected String getRawFieldName() { return super.getRawFieldName() + "String"; }
        @Override protected String getFullName() { return super.getFullName() + "String"; }
        @Override protected String getRawShortFieldName() { return super.getRawShortFieldName() + "String"; }
        @Override protected String getFieldType() { return "List[String]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override public String getCommandLineTemplate() { return REPEAT_TEMPLATE; }

        @Override
        protected String getExclusiveOf() {
            StringBuilder exclusiveOf = new StringBuilder(super.getExclusiveOf());
            if (exclusiveOf.length() > 0)
                exclusiveOf.append(",");
            exclusiveOf.append(escape(argumentDefinition.fullName));
            return exclusiveOf.toString();
        }
    }

    // if (argumentDefinition.ioType == ArgumentIOType.INPUT)
    // Map all inputs to files.  Handles multi valued files.
    private static class InputArgumentField extends ArgumentDefinitionField {
        public InputArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getFieldType() { return isMultiValued() ? "List[File]" : "File"; }
        @Override protected String getDefaultValue() { return isMultiValued() ? "Nil" : "_"; }
    }

    // if (argumentDefinition.ioType == ArgumentIOType.OUTPUT)
    // Map all outputs to files.
    private static class OutputArgumentField extends ArgumentDefinitionField {
        private final String gatherClass;
        public OutputArgumentField(ArgumentDefinition argumentDefinition, String gatherClass) {
            super(argumentDefinition);
            this.gatherClass = gatherClass;
        }

        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getFieldType() { return "File"; }
        @Override protected String getDefaultValue() { return "_"; }

        @Override public boolean isGather() { return true; }
        @Override protected String getGatherAnnotation() {
            return String.format("@Gather(classOf[%s])%n", gatherClass);
        }
    }

    // if (argumentDefinition.isFlag)
    // Booleans should be set on the commandline only if they are true.
    private static class FlagArgumentField extends ArgumentDefinitionField {
        public FlagArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return boolean.class; }
        @Override protected String getFieldType() { return "Boolean"; }
        @Override protected String getDefaultValue() { return "_"; }
        @Override protected String getCommandLineTemplate() { return FLAG_TEMPLATE; }
    }

    // if (argumentDefinition.isMultiValued)
    // Multi value arguments are mapped to List[] and use repeat.
    private static class MultiValuedArgumentField extends ArgumentDefinitionField {
        public MultiValuedArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }

        @Override protected Class<?> getInnerType() { return mapType(argumentDefinition.componentType); }
        @Override protected String getFieldType() { return String.format("List[%s]", getType(getInnerType())); }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected String getCommandLineTemplate() { return REPEAT_TEMPLATE; }
    }

    // if (!argumentDefinition.required && useOption(argumentDefinition.argumentType))
    // Any optional arguments that are primitives are wrapped in options.
    private static class OptionedArgumentField extends ArgumentDefinitionField {
        private final boolean useFormatter;

        public OptionedArgumentField(ArgumentDefinition argumentDefinition, boolean useFormatter) {
            super(argumentDefinition);
            this.useFormatter = useFormatter;
        }

        @Override protected Class<?> getInnerType() { return mapType(argumentDefinition.argumentType); }
        @Override protected String getFieldType() { return String.format("Option[%s]", getType(getInnerType())); }
        @Override protected String getDefaultValue() { return "None"; }
        @Override protected String getCommandLineTemplate() { return OPTIONAL_TEMPLATE; }
        @Override protected String getCommandLineFormat() {
            return this.useFormatter ? getFieldName(this.getRawFieldName() + "Format") : super.getCommandLineFormat();
        }
    }

    // Any other @Arguments
    private static class DefaultArgumentField extends ArgumentDefinitionField {
        private final boolean useFormatter;

        public DefaultArgumentField(ArgumentDefinition argumentDefinition, boolean useFormatter) {
            super(argumentDefinition);
            this.useFormatter = useFormatter;
        }

        @Override protected Class<?> getInnerType() { return mapType(argumentDefinition.argumentType); }
        @Override protected String getFieldType() { return getType(getInnerType()); }
        @Override protected String getDefaultValue() { return "_"; }
        @Override protected String getCommandLineFormat() {
            return this.useFormatter ? getFieldName(this.getRawFieldName() + "Format") : super.getCommandLineFormat();
        }
    }

    // Allows the user to specify the track name, track type, and the file.
    public static class RodBindArgumentField extends ArgumentDefinitionField {
        public RodBindArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }
        @Override protected Class<?> getInnerType() { return null; } // RodBind does not need to be imported.
        @Override protected String getFieldType() { return "List[RodBind]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected String getCommandLineTemplate() {
            return " + repeat(\"\", %3$s, format=RodBind.formatCommandLine(\"%1$s\"))";
        }
    }

    // Tagged input_files or other rods.
    public static class InputTaggedFileDefinitionField extends ArgumentDefinitionField {
        public InputTaggedFileDefinitionField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }
        @Override protected Class<?> getInnerType() { return null; } // TaggedFile does not need to be imported.
        @Override protected String getFieldType() { return argumentDefinition.isMultiValued ? "List[File]" :  "File"; }
        @Override protected String getDefaultValue() { return argumentDefinition.isMultiValued ? "Nil" : "_"; }
        @Override protected String getCommandLineTemplate() {
            if (argumentDefinition.isMultiValued) {
                return " + repeat(\"\", %3$s, format=TaggedFile.formatCommandLine(\"%1$s\"))";
            } else if (!argumentDefinition.required) {
                return " + optional(\"\", %3$s, format=TaggedFile.formatCommandLine(\"%1$s\"))";
            } else {
                return " + TaggedFile.formatCommandLine(\"%1$s\")(\"\", %3$s, \"\")";
            }
        }
    }

    // Adds optional inputs for the indexes of any rods added to this function.
    private static class InputIndexesArgumentField extends ArgumentField {
        private final boolean originalIsMultiValued;
        private final String indexFieldName;
        private final String originalFieldName;
        private final String indexSuffix;
        private final String originalSuffix;
        public InputIndexesArgumentField(ArgumentDefinition originalArgumentDefinition, String indexSuffix) {
            this(originalArgumentDefinition, indexSuffix, null);
        }
        public InputIndexesArgumentField(ArgumentDefinition originalArgumentDefinition, String indexSuffix, String originalSuffix) {
            this.originalIsMultiValued = originalArgumentDefinition.isMultiValued;
            this.indexFieldName = originalArgumentDefinition.fullName + "Index" + (originalIsMultiValued ? "es" : "");
            this.originalFieldName = originalArgumentDefinition.fullName;
            this.indexSuffix = indexSuffix;
            this.originalSuffix = originalSuffix;
        }
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Input.class; }
        @Override public String getCommandLineAddition() { return ""; }
        @Override protected String getDoc() {
            return originalIsMultiValued
                    ? "Dependencies on any indexes of " + this.originalFieldName
                    : "Dependencies on the index of " + this.originalFieldName;
        }
        @Override protected String getFullName() { return this.indexFieldName; }
        @Override protected boolean isRequired() { return false; }
        @Override protected String getFieldType() { return "List[File]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getRawFieldName() { return this.indexFieldName; }
        @Override protected String getFreezeFields() {
            if (originalIsMultiValued) {
                if (originalSuffix == null) {
                    return String.format(
                            ("%1$s ++= %2$s" +
                                    ".filter(orig => orig != null)" +
                                    ".map(orig => new File(orig.getPath + \"%3$s\"))%n"),
                            indexFieldName, originalFieldName, indexSuffix);
                } else {
                    return String.format(
                            ("%1$s ++= %2$s" +
                                    ".filter(orig => orig != null && orig.getName.endsWith(\"%4$s\"))" +
                                    ".flatMap(orig => Array(" +
                                    " new File(orig.getPath + \"%3$s\")," +
                                    " new File(orig.getPath.stripSuffix(\"%4$s\") + \"%3$s\") ))%n"),
                            indexFieldName, originalFieldName, indexSuffix, originalSuffix);
                }
            } else {
                if (originalSuffix == null) {
                    return String.format(
                            ("if (%2$s != null)%n  " +
                                    "%1$s :+= new File(%2$s.getPath + \"%3$s\")%n"),
                            indexFieldName, originalFieldName, indexSuffix);
                } else {
                    return String.format(
                            ("if (%2$s != null && %2$s.getName.endsWith(\"%4$s\"))%n  " +
                                    "%1$s ++= Array(" +
                                    " new File(%2$s.getPath + \"%3$s\")," +
                                    " new File(%2$s.getPath.stripSuffix(\"%4$s\") + \"%3$s\") )%n"),
                            indexFieldName, originalFieldName, indexSuffix, originalSuffix);
                }
            }
        }
    }

    // Tracks an automatically generated index
    private static abstract class OutputIndexArgumentField extends ArgumentField {
        protected final String indexFieldName;
        protected final String originalFieldName;
        public OutputIndexArgumentField(ArgumentDefinition originalArgumentDefinition) {
            this.indexFieldName = originalArgumentDefinition.fullName + "Index";
            this.originalFieldName = originalArgumentDefinition.fullName;
        }
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Output.class; }
        @Override public String getCommandLineAddition() { return ""; }
        @Override protected String getDoc() { return "Automatically generated index for " + this.originalFieldName; }
        @Override protected String getFullName() { return this.indexFieldName; }
        @Override protected boolean isRequired() { return false; }
        @Override protected String getFieldType() { return "File"; }
        @Override protected String getDefaultValue() { return "_"; }
        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getRawFieldName() { return this.indexFieldName; }

        @Override public boolean isGather() { return true; }
        @Override protected String getGatherAnnotation() {
            return String.format("@Gather(classOf[AutoIndexGatherFunction])%n");
        }
    }

    private static class VCFWriterIndexArgumentField extends OutputIndexArgumentField {
        public VCFWriterIndexArgumentField(ArgumentDefinition originalArgumentDefinition) {
            super(originalArgumentDefinition);
        }
        @Override protected String getFreezeFields() {
            return String.format(
                    ("if (%2$s != null)%n" +
                            "  if (!org.broadinstitute.sting.gatk.io.stubs.VCFWriterArgumentTypeDescriptor.isCompressed(%2$s.getPath))%n" +
                            "    %1$s = new File(%2$s.getPath + \"%3$s\")%n"),
                    indexFieldName, originalFieldName, Tribble.STANDARD_INDEX_EXTENSION);
        }
    }

    private static class SAMFileWriterIndexArgumentField extends OutputIndexArgumentField {
        public SAMFileWriterIndexArgumentField(ArgumentDefinition originalArgumentDefinition) {
            super(originalArgumentDefinition);
        }
        @Override protected String getFreezeFields() {
            return String.format(
                    ("if (%2$s != null)%n" +
                            "  if (!%3$s)%n" +
                            "    %1$s = new File(%2$s.getPath.stripSuffix(\".bam\") + \"%4$s\")%n"),
                    indexFieldName, originalFieldName, SAMFileWriterArgumentTypeDescriptor.DISABLE_INDEXING_FULLNAME, BAMIndex.BAMIndexSuffix);
        }
    }

    // Allows setting the format for floats and doubles
    private static class FormatterArgumentField extends ArgumentField {
        private final ArgumentField argumentField;
        public FormatterArgumentField(ArgumentField argumentField) {
            this.argumentField = argumentField;
        }
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Argument.class; }
        @Override public String getCommandLineAddition() { return ""; }
        @Override protected String getDoc() { return "Format string for " + this.argumentField.getFullName(); }
        @Override protected String getFullName() { return this.argumentField.getFullName() + "Format"; }
        @Override protected boolean isRequired() { return false; }
        @Override protected String getFieldType() { return "String"; }
        @Override protected String getDefaultValue() { return "\"%s\""; }
        @Override protected Class<?> getInnerType() { return String.class; }
        @Override protected String getRawFieldName() { return this.argumentField.getRawFieldName() + "Format"; }
    }
}

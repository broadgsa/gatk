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

package org.broadinstitute.gatk.queue.extensions.gatk;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.tribble.Tribble;
import org.broadinstitute.gatk.engine.arguments.GATKArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

import java.io.File;
import java.lang.annotation.Annotation;
import java.lang.reflect.Field;
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
                    "%5$sdef %3$s = this.%1$s%n" +
                    "%n" +
                    "/**%n" +
                    " * Short name of %1$s%n" +
                    " * @param value Short name of %1$s%n" +
                    " */%n" +
                    "%5$sdef %4$s(value: %2$s) { this.%1$s = value }%n",
                    getFieldName(),
                    getFieldType(),
                    getShortFieldGetter(),
                    getShortFieldSetter(),
                    getPrivacy());
    }

    protected static final String REQUIRED_TEMPLATE = " + required(\"%1$s\", %3$s, spaceSeparated=true, escape=true, format=%2$s)";
    protected static final String REPEAT_TEMPLATE = " + repeat(\"%1$s\", %3$s, spaceSeparated=true, escape=true, format=%2$s)";
    protected static final String OPTIONAL_TEMPLATE = " + optional(\"%1$s\", %3$s, spaceSeparated=true, escape=true, format=%2$s)";
    protected static final String FLAG_TEMPLATE = " + conditional(%3$s, \"%1$s\", escape=true, format=%2$s)";

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
        List<ArgumentField> argumentFields = new ArrayList<>();
        for (ArgumentSource argumentSource: parsingEngine.extractArgumentSources(classType))
            if (!argumentSource.isDeprecated()) {
                String gatherer = null;
                if (argumentSource.field.isAnnotationPresent(Gather.class)) {
                    Gather gather = argumentSource.field.getAnnotation(Gather.class);
                    if(! "".equals(gather.className()))
                        gatherer = gather.className();
                    else
                        gatherer = gather.value().getName();
                }
                for (ArgumentDefinition argumentDefinition: argumentSource.createArgumentDefinitions())
                    argumentFields.addAll(getArgumentFields(argumentDefinition, gatherer));
            }
        return argumentFields;
    }

    public static String getArgumentFullName(final Class<?> collection, final String fieldName) {
        try {
            final Field field = collection.getField(fieldName);
            final Argument arg = field.getAnnotation(Argument.class);
            if (arg != null)
                return arg.fullName();
            final Input inputAnnotation = field.getAnnotation(Input.class);
            if (inputAnnotation != null)
                return inputAnnotation.fullName();
            final Output outputAnnotation = field.getAnnotation(Output.class);
            if (outputAnnotation != null)
                return outputAnnotation.fullName();
        } catch (NoSuchFieldException e) {
            throw new IllegalStateException(String.format("Can't find field %s in ArgumentCollection %s", fieldName, collection.getSimpleName()), e);
        }
        throw new IllegalStateException(String.format("Field %s in class %s is not annotated as an argument", fieldName, collection.getName()));
    }

    private static final List<String> intervalFields = new ArrayList<>();
    private static final String inputFileArgument = getArgumentFullName(GATKArgumentCollection.class, "samFiles");

    static {
        intervalFields.add(getArgumentFullName(IntervalArgumentCollection.class, "intervals"));
        intervalFields.add(getArgumentFullName(IntervalArgumentCollection.class, "excludeIntervals"));
    }

    private static List<? extends ArgumentField> getArgumentFields(ArgumentDefinition argumentDefinition, String gatherer) {
        if (intervalFields.contains(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(
                    new IntervalFileArgumentField(argumentDefinition),
                    new IntervalStringArgumentField(argumentDefinition));

        } else if (NumThreadsArgumentField.NUM_THREADS_FIELD.equals(argumentDefinition.fullName)) {
            return Arrays.asList(new NumThreadsArgumentField(argumentDefinition));

        } else if (inputFileArgument.equals(argumentDefinition.fullName) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(new InputTaggedFileDefinitionField(argumentDefinition), new InputIndexesArgumentField(argumentDefinition, BAMIndex.BAMIndexSuffix, ".bam"));

        } else if ((RodBinding.class.equals(argumentDefinition.argumentType) || RodBinding.class.equals(argumentDefinition.componentType) || RodBindingCollection.class.equals(argumentDefinition.componentType)) && argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Arrays.asList(new InputTaggedFileDefinitionField(argumentDefinition), new InputIndexesArgumentField(argumentDefinition, Tribble.STANDARD_INDEX_EXTENSION));

        } else if (argumentDefinition.ioType == ArgumentIOType.INPUT) {
            return Collections.singletonList(new InputArgumentField(argumentDefinition));

        } else if (argumentDefinition.ioType == ArgumentIOType.OUTPUT) {

            List<ArgumentField> fields = new ArrayList<>();

            String gatherClass;

            // one can set the specific gatherer to use by adding @Gather before any output argument.
            // For example (used to be part of UG):
            //      @Gather(className = "org.broadinstitute.gatk.queue.extensions.gatk.CatVariantsGatherer")
            //      @Output(doc="File to which variants should be written",required=true)
            //      protected VariantContextWriter writer = null;
            if (gatherer != null)
                gatherClass = gatherer;
            else if (SAMFileWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                gatherClass = "BamGatherFunction";
            else if (VariantContextWriter.class.isAssignableFrom(argumentDefinition.argumentType))
                gatherClass = "CatVariantsGatherer"; // used to be "VcfGatherFunction";
            else
                gatherClass = "org.broadinstitute.gatk.queue.function.scattergather.SimpleTextGatherFunction";

            fields.add(new OutputArgumentField(argumentDefinition, gatherClass));

            if (SAMFileWriter.class.isAssignableFrom(argumentDefinition.argumentType)) {
                fields.add(new SAMFileWriterIndexArgumentField(argumentDefinition));
                fields.add(new SAMFileWriterMD5ArgumentField(argumentDefinition));
            }
            else if (VariantContextWriter.class.isAssignableFrom(argumentDefinition.argumentType)) {
                fields.add(new VCFWriterIndexArgumentField(argumentDefinition));
            }

            return fields;

        } else if (argumentDefinition.isFlag) {
            return Collections.singletonList(new FlagArgumentField(argumentDefinition));

        } else if (argumentDefinition.isMultiValued) {
            return Collections.singletonList(new MultiValuedArgumentField(argumentDefinition));

        } else if (!argumentDefinition.required && useOption(argumentDefinition.argumentType)) {
            boolean useFormat = useFormatter(argumentDefinition.argumentType);
            List<ArgumentField> fields = new ArrayList<>();
            ArgumentField field = new OptionedArgumentField(argumentDefinition, useFormat);
            fields.add(field);
            if (useFormat) fields.add(new FormatterArgumentField(field));
            return fields;

        } else {
            boolean useFormat = useFormatter(argumentDefinition.argumentType);
            List<ArgumentField> fields = new ArrayList<>();
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
        @Override protected String getFieldType() { return "Seq[String]"; }
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
        @Override protected String getFieldType() { return isMultiValued() ? "Seq[File]" : "File"; }
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
        @Override protected String getFieldType() { return String.format("Seq[%s]", getType(getInnerType())); }
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
    public static class NumThreadsArgumentField extends OptionedArgumentField {
        public static final String NUM_THREADS_FIELD = getArgumentFullName(GATKArgumentCollection.class, "numberOfDataThreads");
        public static final String NCT_FIELD = getArgumentFullName(GATKArgumentCollection.class, "numberOfCPUThreadsPerDataThread");

        public NumThreadsArgumentField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition, false);
        }

        @Override
        protected String getFreezeFields() {
            return String.format("if (%1$s.isDefined) nCoresRequest = %1$s%nif (%2$s.isDefined) nCoresRequest = Some(nCoresRequest.getOrElse(1) * %2$s.getOrElse(1))%n",
                    NUM_THREADS_FIELD, NCT_FIELD);
        }
    }

    // Tagged input_files or other rods.
    public static class InputTaggedFileDefinitionField extends ArgumentDefinitionField {
        public InputTaggedFileDefinitionField(ArgumentDefinition argumentDefinition) {
            super(argumentDefinition);
        }
        @Override protected Class<?> getInnerType() { return null; } // TaggedFile does not need to be imported.
        @Override protected String getFieldType() { return argumentDefinition.isMultiValued ? "Seq[File]" :  "File"; }
        @Override protected String getDefaultValue() { return argumentDefinition.isMultiValued ? "Nil" : "_"; }
        @Override protected String getCommandLineTemplate() {
            if (argumentDefinition.isMultiValued) {
                return " + repeat(\"%1$s\", %3$s, formatPrefix=TaggedFile.formatCommandLineParameter, spaceSeparated=true, escape=true, format=%2$s)";
            } else if (!argumentDefinition.required) {
                return " + optional(TaggedFile.formatCommandLineParameter(\"%1$s\", %3$s), %3$s, spaceSeparated=true, escape=true, format=%2$s)";
            } else {
                return " + required(TaggedFile.formatCommandLineParameter(\"%1$s\", %3$s), %3$s, spaceSeparated=true, escape=true, format=%2$s)";
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
        @Override protected String getFieldType() { return "Seq[File]"; }
        @Override protected String getDefaultValue() { return "Nil"; }
        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getRawFieldName() { return this.indexFieldName; }
        @Override protected String getPrivacy() { return "private "; }
        @Override protected String getFreezeFields() {
            if (originalIsMultiValued) {
                if (originalSuffix == null) {
                    return String.format(
                            ("%1$s ++= %2$s" +
                                    ".filter(orig => orig != null && (!orig.getName.endsWith(\".list\")))" +
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

    // Tracks an automatically generated index, md5, etc.
    private static abstract class AuxilliaryOutputArgumentField extends ArgumentField {
        protected final String originalFieldName;
        protected final String auxFieldName;
        protected final String auxFieldLabel;
        public AuxilliaryOutputArgumentField(ArgumentDefinition originalArgumentDefinition, String auxFieldLabel) {
            this.originalFieldName = originalArgumentDefinition.fullName;
            this.auxFieldName = originalArgumentDefinition.fullName + auxFieldLabel;
            this.auxFieldLabel = auxFieldLabel;
        }
        @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Output.class; }
        @Override public String getCommandLineAddition() { return ""; }
        @Override protected String getDoc() { return String.format("Automatically generated %s for %s", auxFieldLabel.toLowerCase(), this.originalFieldName); }
        @Override protected String getFullName() { return this.auxFieldName; }
        @Override protected boolean isRequired() { return false; }
        @Override protected String getFieldType() { return "File"; }
        @Override protected String getDefaultValue() { return "_"; }
        @Override protected Class<?> getInnerType() { return File.class; }
        @Override protected String getRawFieldName() { return this.auxFieldName; }
        @Override protected String getPrivacy() { return "private "; }

        @Override public boolean isGather() { return true; }
        @Override protected String getGatherAnnotation() {
            return String.format("@Gather(enabled=false)%n");
        }
    }

    private static class VCFWriterIndexArgumentField extends AuxilliaryOutputArgumentField {
        public VCFWriterIndexArgumentField(ArgumentDefinition originalArgumentDefinition) {
            super(originalArgumentDefinition, "Index");
        }
        @Override protected String getFreezeFields() {
            return String.format(
                    ("if (%2$s != null && !org.broadinstitute.gatk.utils.io.IOUtils.isSpecialFile(%2$s))%n" +
                            "  if (!org.broadinstitute.gatk.utils.commandline.ArgumentTypeDescriptor.isCompressed(%2$s.getPath))%n" +
                            "    %1$s = new File(%2$s.getPath + \"%3$s\")%n"),
                    auxFieldName, originalFieldName, Tribble.STANDARD_INDEX_EXTENSION);
        }
    }

    private static class SAMFileWriterIndexArgumentField extends AuxilliaryOutputArgumentField {
        public SAMFileWriterIndexArgumentField(ArgumentDefinition originalArgumentDefinition) {
            super(originalArgumentDefinition, "Index");
        }
        @Override protected String getFreezeFields() {
            return String.format(
                    ("if (%2$s != null && !org.broadinstitute.gatk.utils.io.IOUtils.isSpecialFile(%2$s))%n" +
                            "  if (!%3$s)%n" +
                            "    %1$s = new File(%2$s.getPath.stripSuffix(\".bam\") + \"%4$s\")%n"),
                    auxFieldName, originalFieldName, getArgumentFullName(GATKArgumentCollection.class, "disableBAMIndexing"), BAMIndex.BAMIndexSuffix);
        }
    }

    private static class SAMFileWriterMD5ArgumentField extends AuxilliaryOutputArgumentField {
        public SAMFileWriterMD5ArgumentField(ArgumentDefinition originalArgumentDefinition) {
            super(originalArgumentDefinition, "MD5");
        }
        @Override protected String getFreezeFields() {
            return String.format(
                    ("if (%2$s != null && !org.broadinstitute.gatk.utils.io.IOUtils.isSpecialFile(%2$s))%n" +
                            "  if (%3$s)%n" +
                            "    %1$s = new File(%2$s.getPath + \"%4$s\")%n"),
                    auxFieldName, originalFieldName, getArgumentFullName(GATKArgumentCollection.class, "enableBAMmd5"), ".md5");
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

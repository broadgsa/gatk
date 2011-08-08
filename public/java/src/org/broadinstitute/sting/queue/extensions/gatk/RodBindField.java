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

import org.broadinstitute.sting.commandline.Input;
import org.broadinstitute.sting.gatk.WalkerManager;
import org.broadinstitute.sting.gatk.refdata.tracks.RMDTrackBuilder;
import org.broadinstitute.sting.gatk.walkers.RMD;
import org.broadinstitute.sting.gatk.walkers.Walker;

import java.io.File;
import java.lang.annotation.Annotation;
import java.util.ArrayList;
import java.util.List;

/**
 * Allows user to specify the rod file but locks in the track name and the track type.
 */
public class RodBindField extends ArgumentField {
    public static final String ROD_BIND_FIELD = "rodBind";

    private final String trackName;
    private final String typeName;
    private final List<RodBindField> relatedFields;
    private final boolean isRequired;

    public RodBindField(String trackName, String typeName, List<RodBindField> relatedFields, boolean isRequired) {
        this.trackName = trackName;
        this.typeName = typeName;
        this.relatedFields = relatedFields;
        this.isRequired = isRequired;
    }

    @SuppressWarnings("unchecked")
    @Override protected Class<? extends Annotation> getAnnotationIOClass() { return Input.class; }
    @Override protected Class<?> getInnerType() { return File.class; }
    @Override protected String getFullName() { return escape(getRawFieldName()); }
    @Override protected String getFieldType() { return "File"; }
    @Override protected String getDefaultValue() { return "_"; }
    @Override protected String getRawFieldName() { return this.trackName + this.typeName; }
    @Override protected String getDoc() { return escape(this.typeName + " " + this.trackName); }
    @Override protected boolean isRequired() { return this.isRequired; }

    @Override public String getCommandLineAddition() {
        // TODO: Stop allowing the generic "rodBind" triplets to satisfy the requirement after @Requires are fixed.
        return String.format(" + optional(\" -B:%s,%s \", %s)",
        /*
        return String.format(this.useOption()
                ? " + optional(\" -B:%s,%s \", %s)"
                : " + \" -B:%s,%s \" + %s",
        */
                this.trackName, this.typeName, getFieldName());
    }

    private boolean useOption() {
        return !this.isRequired || (relatedFields.size() > 1);
    }

    @Override protected String getExclusiveOf() {
        StringBuilder exclusiveOf = new StringBuilder();
        // TODO: Stop allowing the generic "rodBind" triplets to satisfy the requirement after @Requires are fixed.
        if (this.isRequired)
            exclusiveOf.append(ROD_BIND_FIELD);
        for (RodBindField relatedField: relatedFields)
            if (relatedField != this) {
                if (exclusiveOf.length() > 0)
                    exclusiveOf.append(",");
                exclusiveOf.append(relatedField.getFieldName());
            }
        return exclusiveOf.toString();
    }
//
//    public static List<ArgumentField> getRodArguments(Class<? extends Walker> walkerClass, RMDTrackBuilder trackBuilder) {
//        List<ArgumentField> argumentFields = new ArrayList<ArgumentField>();
//
//        List<RMD> requires = WalkerManager.getRequiredMetaData(walkerClass);
//        List<RMD> allows = WalkerManager.getAllowsMetaData(walkerClass);
//
//        for (RMD required: requires) {
//            List<RodBindField> fields = new ArrayList<RodBindField>();
//            String trackName = required.name();
//            if ("*".equals(trackName)) {
//                // TODO: Add the field triplet for name=* after @Allows and @Requires are fixed on walkers
//                //fields.add(new RodBindArgumentField(argumentDefinition, true));
//            } else {
//                for (String typeName: trackBuilder.getFeatureManager().getTrackRecordTypeNames(required.type()))
//                    fields.add(new RodBindField(trackName, typeName, fields, true));
//            }
//            argumentFields.addAll(fields);
//        }
//
//        for (RMD allowed: allows) {
//            List<RodBindField> fields = new ArrayList<RodBindField>();
//            String trackName = allowed.name();
//            if ("*".equals(trackName)) {
//                // TODO: Add the field triplet for name=* after @Allows and @Requires are fixed on walkers
//                //fields.add(new RodBindArgumentField(argumentDefinition, false));
//            } else {
//                for (String typeName: trackBuilder.getFeatureManager().getTrackRecordTypeNames(allowed.type()))
//                    fields.add(new RodBindField(trackName, typeName, fields, true));
//            }
//            argumentFields.addAll(fields);
//        }
//
//        return argumentFields;
//    }
}

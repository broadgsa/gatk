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

package org.broadinstitute.gatk.utils.jna.drmaa.v1_0;

import com.sun.jna.Memory;
import com.sun.jna.NativeLong;
import com.sun.jna.Pointer;
import com.sun.jna.StringArray;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;
import org.ggf.drmaa.*;

import java.text.ParseException;
import java.util.*;

/**
 * JNA mapping from Java to C DRMAA binding.
 * See: Java and C Binding Documents on http://drmaa.org
 */
public class JnaSession implements Session {
    private static final PartialTimestampFormat PARTIAL_TIMESTAMP_FORMAT = new PartialTimestampFormat();
    private static final ThreadLocal<Memory> threadError = new ThreadLocal<Memory>() {
        @Override
        protected Memory initialValue() {
            return new Memory(LibDrmaa.DRMAA_ERROR_STRING_BUFFER);
        }
    };

    @Override
    public void init(String contact) throws DrmaaException {
        checkError(LibDrmaa.drmaa_init(contact, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
    }

    @Override
    public void exit() throws DrmaaException {
        checkError(LibDrmaa.drmaa_exit(getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
    }

    @Override
    public JobTemplate createJobTemplate() throws DrmaaException {
        PointerByReference jtRef = new PointerByReference();
        checkError(LibDrmaa.drmaa_allocate_job_template(jtRef, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        return new JnaJobTemplate(this, jtRef.getValue());
    }

    @Override
    public void deleteJobTemplate(JobTemplate jobTemplate) throws DrmaaException {
        JnaJobTemplate jnaJobTemplate = (JnaJobTemplate) jobTemplate;
        checkError(LibDrmaa.drmaa_delete_job_template(jnaJobTemplate.getPointer(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
    }

    @Override
    public String runJob(JobTemplate jobTemplate) throws DrmaaException {
        Memory jobId = new Memory(LibDrmaa.DRMAA_JOBNAME_BUFFER);
        JnaJobTemplate jnaJobTemplate = (JnaJobTemplate) jobTemplate;
        checkError(LibDrmaa.drmaa_run_job(jobId, LibDrmaa.DRMAA_JOBNAME_BUFFER_LEN, jnaJobTemplate.getPointer(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        return jobId.getString(0);
    }

    @Override
    public List runBulkJobs(JobTemplate jobTemplate, int start, int end, int incr) throws DrmaaException {
        PointerByReference jobIds = new PointerByReference();
        JnaJobTemplate jnaJobTemplate = (JnaJobTemplate) jobTemplate;
        checkError(LibDrmaa.drmaa_run_bulk_jobs(jobIds, jnaJobTemplate.getPointer(), start, end, incr, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        try {
            return getJobIds(jobIds);
        } finally {
            releaseJobIds(jobIds);
        }
    }

    @Override
    public void control(String jobId, int action) throws DrmaaException {
        checkError(LibDrmaa.drmaa_control(jobId, action, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
    }

    @SuppressWarnings("unchecked")
    @Override
    public void synchronize(List list, long timeout, boolean dispose) throws DrmaaException {
        StringArray jobIds = new StringArray((String[]) list.toArray(new String[list.size()]));
        checkError(LibDrmaa.drmaa_synchronize(jobIds, new NativeLong(timeout), dispose ? 1 : 0, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
    }

    @Override
    public JobInfo wait(String jobId, long timeout) throws DrmaaException {
        Memory jobIdOut = new Memory(LibDrmaa.DRMAA_JOBNAME_BUFFER);
        IntByReference stat = new IntByReference();
        PointerByReference rusage = new PointerByReference();
        IntByReference exited = new IntByReference();
        IntByReference exitStatus = new IntByReference();
        IntByReference signaled = new IntByReference();
        Memory signal = new Memory(LibDrmaa.DRMAA_SIGNAL_BUFFER);
        IntByReference coreDumped = new IntByReference();
        IntByReference aborted = new IntByReference();

        int errnum;

        errnum = LibDrmaa.drmaa_wait(jobId, jobIdOut, LibDrmaa.DRMAA_JOBNAME_BUFFER_LEN, stat, new NativeLong(timeout), rusage, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN);

        Map<String, String> rusageMap;
        if (errnum == LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_RUSAGE) {
            rusageMap = null;
        } else {
            try {
                rusageMap = collectionToMap(getAttrValues(rusage));
            } finally {
                releaseAttrValues(rusage);
            }
        }

        checkError(LibDrmaa.drmaa_wifexited(exited, stat.getValue(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));

        if (exited.getValue() != 0) {
            checkError(LibDrmaa.drmaa_wexitstatus(exitStatus, stat.getValue(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        }

        checkError(LibDrmaa.drmaa_wifsignaled(signaled, stat.getValue(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));

        if (signaled.getValue() != 0) {
            checkError(LibDrmaa.drmaa_wtermsig(signal, LibDrmaa.DRMAA_SIGNAL_BUFFER_LEN, stat.getValue(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
            checkError(LibDrmaa.drmaa_wcoredump(coreDumped, stat.getValue(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        }

        checkError(LibDrmaa.drmaa_wifaborted(aborted, stat.getValue(), getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));

        return new JnaJobInfo(jobIdOut.getString(0), rusageMap, exited.getValue() != 0, exitStatus.getValue(),
                signaled.getValue() != 0, signal.getString(0), coreDumped.getValue() != 0, aborted.getValue() != 0);
    }

    @Override
    public int getJobProgramStatus(String jobId) throws DrmaaException {
        IntByReference remotePs = new IntByReference();
        checkError(LibDrmaa.drmaa_job_ps(jobId, remotePs, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        return remotePs.getValue();
    }

    @Override
    public String getContact() {
        Memory contact = new Memory(LibDrmaa.DRMAA_CONTACT_BUFFER);
        try {
            checkError(LibDrmaa.drmaa_get_contact(contact, LibDrmaa.DRMAA_CONTACT_BUFFER_LEN, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        } catch (DrmaaException e) {
            // DRMAA spec says this method should throw DrmaaException.
            // Why doesn't interface implement this?
            throw new RuntimeException(e);
        }
        return contact.getString(0);
    }

    @Override
    public Version getVersion() {
        IntByReference major = new IntByReference();
        IntByReference minor = new IntByReference();
        try {
            checkError(LibDrmaa.drmaa_version(major, minor, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        } catch (DrmaaException e) {
            // DRMAA spec says this method should throw DrmaaException.
            // Why doesn't interface implement this?
            throw new RuntimeException(e);
        }
        return new Version(major.getValue(), minor.getValue());
    }

    @Override
    public String getDrmSystem() {
        Memory drmSystem = new Memory(LibDrmaa.DRMAA_DRM_SYSTEM_BUFFER);
        try {
            checkError(LibDrmaa.drmaa_get_DRM_system(drmSystem, LibDrmaa.DRMAA_DRM_SYSTEM_BUFFER_LEN, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        } catch (DrmaaException e) {
            // DRMAA spec says this method should throw DrmaaException.
            // Why doesn't interface implement this?
            throw new RuntimeException(e);
        }
        return drmSystem.getString(0);
    }

    @Override
    public String getDrmaaImplementation() {
        Memory drmaaImplementation = new Memory(LibDrmaa.DRMAA_DRMAA_IMPLEMENTATION_BUFFER);
        try {
            checkError(LibDrmaa.drmaa_get_DRMAA_implementation(drmaaImplementation, LibDrmaa.DRMAA_DRMAA_IMPLEMENTATION_BUFFER_LEN, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        } catch (DrmaaException e) {
            // DRMAA spec says this method should throw DrmaaException.
            // Why doesn't interface implement this?
            throw new RuntimeException(e);
        }
        return drmaaImplementation.getString(0);
    }

    public static void setAttribute(Pointer jt, String name, String value) throws DrmaaException {
        if (getAttrNames().contains(name)) {
            checkError(LibDrmaa.drmaa_set_attribute(jt, name, value, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        }
        else {
            throw new InvalidAttributeValueException("Attribute " + name + " is not supported by this implementation of DRMAA");
        }
    }

    public static String getAttribute(Pointer jt, String name) throws DrmaaException {
        if (getAttrNames().contains(name)) {
            Memory attrBuffer = new Memory(LibDrmaa.DRMAA_ATTR_BUFFER);
            checkError(LibDrmaa.drmaa_get_attribute(jt, name, attrBuffer, LibDrmaa.DRMAA_ATTR_BUFFER_LEN, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
            return attrBuffer.getString(0);
        }
        else {
            throw new InvalidAttributeValueException("Attribute " + name + " is not supported by this implementation of DRMAA");
        }
    }

    public static void setVectorAttribute(Pointer jt, String name, Collection<String> values) throws DrmaaException {
        StringArray valuesArray = new StringArray(values.toArray(new String[values.size()]));
        checkError(LibDrmaa.drmaa_set_vector_attribute(jt, name, valuesArray, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
    }

    public static List<String> getVectorAttribute(Pointer jt, String name) throws DrmaaException {
        PointerByReference values = new PointerByReference();
        checkError(LibDrmaa.drmaa_get_vector_attribute(jt, name, values, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        try {
            return getAttrValues(values);
        } finally {
            releaseAttrValues(values);
        }
    }

    public static void setPartialTime(Pointer jt, String name, PartialTimestamp partialTimestamp) throws DrmaaException {
        setAttribute(jt, name, PARTIAL_TIMESTAMP_FORMAT.format(partialTimestamp));
    }

    public static PartialTimestamp getPartialTime(Pointer jt, String name) throws DrmaaException {
        String time = getAttribute(jt, name);
        if (time == null)
            return null;
        try {
            return PARTIAL_TIMESTAMP_FORMAT.parse(time);
        } catch (ParseException e) {
            throw new InternalException(name + " property is unparsable");
        }
    }

    public static Set<String> getAttrNames() throws DrmaaException {
        PointerByReference values = new PointerByReference();
        checkError(LibDrmaa.drmaa_get_attribute_names(values, getError(), LibDrmaa.DRMAA_ERROR_STRING_BUFFER_LEN));
        try {
            return new LinkedHashSet<String>(getAttrNames(values));
        } finally {
            releaseAttrNames(values);
        }
    }

    public static Collection<String> mapToCollection(Map<String, String> map) {
        Collection<String> collection = new LinkedHashSet<String>();
        for (Map.Entry<String, String> entry: map.entrySet())
            collection.add(entry.getKey() + "=" + entry.getValue());
        return collection;
    }

    public static Map<String, String> collectionToMap(Collection<String> list) {
        Map<String, String> map = new LinkedHashMap<String, String>();
        for (String entry: list) {
            if (entry == null)
                continue;
            int equals = entry.indexOf('=');
            if (equals < 0)
                continue;
            map.put(entry.substring(0, equals), entry.substring(equals + 1));
        }
        return map;
    }

    public static String formatLimit(long secs) {
        long seconds = (secs % 60);
        long minutes = (secs / 60) % 60;
        long hours = (secs / 3600);
        return String.format("%d:%02d:%02d", hours, minutes, seconds);
    }

    public static long parseLimit(String limit) {
        long seconds = 0;
        if (limit != null) {
            for (String token: limit.split(":")) {
                seconds *= 60;
                seconds += Long.parseLong(token);
            }
        }
        return seconds;
    }

    private static List<String> getAttrNames(PointerByReference names) throws DrmaaException {
        List<String> namesList = new ArrayList<String>();
        IntByReference size = new IntByReference();
        int errnum;

        errnum = LibDrmaa.drmaa_get_num_attr_names(names.getValue(), size);
        checkError(errnum, "unable to get attribute names");
        int num = size.getValue();

        Memory value = new Memory(LibDrmaa.DRMAA_ATTR_BUFFER);
        for (int i = 1; i <= num; i++) {
            errnum = LibDrmaa.drmaa_get_next_attr_name(names.getValue(), value, LibDrmaa.DRMAA_ATTR_BUFFER_LEN);
            checkError(errnum, "unable to get attribute name " + i);
            if (errnum == LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_MORE_ELEMENTS)
                break;
            namesList.add(value.getString(0));
        }

        return namesList;
    }

    private static List<String> getAttrValues(PointerByReference values) throws DrmaaException {
        List<String> valuesList = new ArrayList<String>();
        IntByReference size = new IntByReference();
        int errnum;

        errnum = LibDrmaa.drmaa_get_num_attr_values(values.getValue(), size);
        checkError(errnum, "unable to get attribute values");
        int num = size.getValue();

        Memory value = new Memory(LibDrmaa.DRMAA_ATTR_BUFFER);
        for (int i = 1; i <= num; i++) {
            errnum = LibDrmaa.drmaa_get_next_attr_value(values.getValue(), value, LibDrmaa.DRMAA_ATTR_BUFFER_LEN);
            checkError(errnum, "unable to get attribute value " + i);
            if (errnum == LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_MORE_ELEMENTS)
                break;
            valuesList.add(value.getString(0));
        }

        return valuesList;
    }

    private static List<String> getJobIds(PointerByReference jobIds) throws DrmaaException {
        List<String> jobIdsList = new ArrayList<String>();
        IntByReference size = new IntByReference();
        int errnum;

        errnum = LibDrmaa.drmaa_get_num_job_ids(jobIds.getValue(), size);
        checkError(errnum, "unable to get jobIds");
        int num = size.getValue();

        Memory value = new Memory(LibDrmaa.DRMAA_JOBNAME_BUFFER);
        for (int i = 1; i <= num; i++) {
            errnum = LibDrmaa.drmaa_get_next_job_id(jobIds.getValue(), value, LibDrmaa.DRMAA_JOBNAME_BUFFER_LEN);
            checkError(errnum, "unable to get jobId " + i);
            if (errnum == LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_MORE_ELEMENTS)
                break;
            jobIdsList.add(value.getString(0));
        }

        return jobIdsList;
    }

    private static void releaseAttrNames(PointerByReference names) throws DrmaaException {
        LibDrmaa.drmaa_release_attr_names(names.getValue());
    }

    private static void releaseAttrValues(PointerByReference values) throws DrmaaException {
        LibDrmaa.drmaa_release_attr_values(values.getValue());
    }

    private static void releaseJobIds(PointerByReference jobIds) throws DrmaaException {
        LibDrmaa.drmaa_release_job_ids(jobIds.getValue());
    }

    private static Memory getError() {
        return threadError.get();
    }

    private static void checkError(int errnum) throws DrmaaException {
        if (errnum != LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS)
            checkError(errnum, getError().getString(0));
    }

    private static void checkError(int errnum, String error) throws DrmaaException {
        switch (errnum) {
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUCCESS:
                break;
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_INTERNAL_ERROR:
                throw new InternalException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE:
                throw new DrmCommunicationException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_AUTH_FAILURE:
                throw new AuthorizationException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_INVALID_ARGUMENT:
                throw new IllegalArgumentException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_ACTIVE_SESSION:
                throw new NoActiveSessionException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_MEMORY:
                throw new OutOfMemoryError(error);

                /* -------------- init and exit specific --------------- */
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_INVALID_CONTACT_STRING:
                throw new InvalidContactStringException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_DEFAULT_CONTACT_STRING_ERROR:
                throw new DefaultContactStringException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_DEFAULT_CONTACT_STRING_SELECTED:
                throw new NoDefaultContactStringException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_DRMS_INIT_FAILED:
                throw new DrmsInitException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_ALREADY_ACTIVE_SESSION:
                throw new AlreadyActiveSessionException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_DRMS_EXIT_ERROR:
                throw new DrmsExitException(error);

                /* ---------------- job attributes specific -------------- */
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_INVALID_ATTRIBUTE_FORMAT:
                throw new InvalidAttributeFormatException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE:
                throw new InvalidAttributeValueException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_CONFLICTING_ATTRIBUTE_VALUES:
                throw new ConflictingAttributeValuesException(error);

                /* --------------------- job submission specific -------------- */
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_TRY_LATER:
                throw new TryLaterException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_DENIED_BY_DRM:
                throw new DeniedByDrmException(error);

                /* ------------------------------- job control specific ---------------- */
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_INVALID_JOB:
                throw new InvalidJobException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_RESUME_INCONSISTENT_STATE:
                throw new ResumeInconsistentStateException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_SUSPEND_INCONSISTENT_STATE:
                throw new SuspendInconsistentStateException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_HOLD_INCONSISTENT_STATE:
                throw new HoldInconsistentStateException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_RELEASE_INCONSISTENT_STATE:
                throw new ReleaseInconsistentStateException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_EXIT_TIMEOUT:
                throw new ExitTimeoutException(error);
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_RUSAGE:
                break;
            case LibDrmaa.DRMAA_ERRNO.DRMAA_ERRNO_NO_MORE_ELEMENTS:
                break;
            default:
                throw new IllegalArgumentException(String.format("Unknown error code %d: %s", errnum, error));
        }
    }
}

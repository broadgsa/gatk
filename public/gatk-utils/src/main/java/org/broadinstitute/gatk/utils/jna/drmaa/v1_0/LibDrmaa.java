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

import com.sun.jna.*;
import com.sun.jna.ptr.IntByReference;
import com.sun.jna.ptr.PointerByReference;

@SuppressWarnings("unused")
public class LibDrmaa {
    static {
        Native.register("drmaa");
    }

/* see www.drmaa.org for more details on the DRMAA specification */
/****** DRMAA/-DRMAA_Interface *************************************************
*  NAME
*     DRMAA_Interface -- DRMAA interface
*
*  FUNCTION
*     The enlisted functions specify the C/C++ binding of the DRMAA interface
*     specification.
*
*  SEE ALSO
*     DRMAA/drmaa_get_next_attr_name()
*     DRMAA/drmaa_get_next_attr_value()
*     DRMAA/drmaa_get_next_job_id()
*     DRMAA/drmaa_release_attr_names()
*     DRMAA/drmaa_release_attr_values()
*     DRMAA/drmaa_release_job_ids()
*     DRMAA/drmaa_init()
*     DRMAA/drmaa_exit()
*     DRMAA/drmaa_allocate_job_template()
*     DRMAA/drmaa_delete_job_template()
*     DRMAA/drmaa_set_attribute()
*     DRMAA/drmaa_get_attribute()
*     DRMAA/drmaa_set_vector_attribute()
*     DRMAA/drmaa_get_vector_attribute()
*     DRMAA/drmaa_get_attribute_names()
*     DRMAA/drmaa_get_vector_attribute_names()
*     DRMAA/drmaa_run_job()
*     DRMAA/drmaa_run_bulk_jobs()
*     DRMAA/drmaa_control()
*     DRMAA/drmaa_synchronize()
*     DRMAA/drmaa_wait()
*     DRMAA/drmaa_wifexited()
*     DRMAA/drmaa_wexitstatus()
*     DRMAA/drmaa_wifsignaled()
*     DRMAA/drmaa_wtermsig()
*     DRMAA/drmaa_wcoredump()
*     DRMAA/drmaa_wifaborted()
*     DRMAA/drmaa_job_ps()
*     DRMAA/drmaa_strerror()
*     DRMAA/drmaa_get_contact()
*     DRMAA/drmaa_version()
*     DRMAA/drmaa_get_DRM_system()
*******************************************************************************/

/* ------------------- Constants ------------------- */
/*
 * some not yet agreed buffer length constants
 * these are recommended minimum values
 */

/* drmaa_get_attribute() */
public static final long DRMAA_ATTR_BUFFER = 1024;
public static final NativeLong DRMAA_ATTR_BUFFER_LEN = new NativeLong(DRMAA_ATTR_BUFFER - 1);

/* drmaa_get_contact() */
public static final long DRMAA_CONTACT_BUFFER = 1024;
public static final NativeLong DRMAA_CONTACT_BUFFER_LEN = new NativeLong(DRMAA_CONTACT_BUFFER - 1);

/* drmaa_get_DRM_system() */
public static final long DRMAA_DRM_SYSTEM_BUFFER = 1024;
public static final NativeLong DRMAA_DRM_SYSTEM_BUFFER_LEN = new NativeLong(DRMAA_DRM_SYSTEM_BUFFER - 1);

/* drmaa_get_DRM_system() */
public static final long DRMAA_DRMAA_IMPLEMENTATION_BUFFER = 1024;
public static final NativeLong DRMAA_DRMAA_IMPLEMENTATION_BUFFER_LEN = new NativeLong(DRMAA_DRMAA_IMPLEMENTATION_BUFFER - 1);

/*
 * Agreed buffer length constants
 * these are recommended minimum values
 */
public static final long DRMAA_ERROR_STRING_BUFFER = 1024;
public static final long DRMAA_JOBNAME_BUFFER = 1024;
public static final long DRMAA_SIGNAL_BUFFER = 32;

public static final NativeLong DRMAA_ERROR_STRING_BUFFER_LEN = new NativeLong(DRMAA_ERROR_STRING_BUFFER - 1);
public static final NativeLong DRMAA_JOBNAME_BUFFER_LEN = new NativeLong(DRMAA_JOBNAME_BUFFER - 1);
public static final NativeLong DRMAA_SIGNAL_BUFFER_LEN = new NativeLong(DRMAA_SIGNAL_BUFFER - 1);

/*
 * Agreed constants
 */
public static final NativeLong DRMAA_TIMEOUT_WAIT_FOREVER = new NativeLong(-1);
public static final NativeLong DRMAA_TIMEOUT_NO_WAIT = new NativeLong(0);

public static final String DRMAA_JOB_IDS_SESSION_ANY = "DRMAA_JOB_IDS_SESSION_ANY";
public static final String DRMAA_JOB_IDS_SESSION_ALL = "DRMAA_JOB_IDS_SESSION_ALL";

public static final String DRMAA_SUBMISSION_STATE_ACTIVE = "drmaa_active";
public static final String DRMAA_SUBMISSION_STATE_HOLD = "drmaa_hold";

/*
 * Agreed placeholder names
 */
public static final String DRMAA_PLACEHOLDER_INCR = "$drmaa_incr_ph$";
public static final String DRMAA_PLACEHOLDER_HD = "$drmaa_hd_ph$";
public static final String DRMAA_PLACEHOLDER_WD = "$drmaa_wd_ph$";

/*
 * Agreed names of job template attributes
 */
public static final String DRMAA_REMOTE_COMMAND = "drmaa_remote_command";
public static final String DRMAA_JS_STATE = "drmaa_js_state";
public static final String DRMAA_WD = "drmaa_wd";
public static final String DRMAA_JOB_CATEGORY = "drmaa_job_category";
public static final String DRMAA_NATIVE_SPECIFICATION = "drmaa_native_specification";
public static final String DRMAA_BLOCK_EMAIL = "drmaa_block_email";
public static final String DRMAA_START_TIME = "drmaa_start_time";
public static final String DRMAA_JOB_NAME = "drmaa_job_name";
public static final String DRMAA_INPUT_PATH = "drmaa_input_path";
public static final String DRMAA_OUTPUT_PATH = "drmaa_output_path";
public static final String DRMAA_ERROR_PATH = "drmaa_error_path";
public static final String DRMAA_JOIN_FILES = "drmaa_join_files";
public static final String DRMAA_TRANSFER_FILES = "drmaa_transfer_files";
public static final String DRMAA_DEADLINE_TIME = "drmaa_deadline_time";
public static final String DRMAA_WCT_HLIMIT = "drmaa_wct_hlimit";
public static final String DRMAA_WCT_SLIMIT = "drmaa_wct_slimit";
public static final String DRMAA_DURATION_HLIMIT = "drmaa_duration_hlimit";
public static final String DRMAA_DURATION_SLIMIT = "drmaa_duration_slimit";

/* names of job template vector attributes */
public static final String DRMAA_V_ARGV = "drmaa_v_argv";
public static final String DRMAA_V_ENV = "drmaa_v_env";
public static final String DRMAA_V_EMAIL = "drmaa_v_email";

/*
 * DRMAA errno values
 *
 * do not touch these values are agreed !!!
 */
public static interface DRMAA_ERRNO {
   /* -------------- these are relevant to all sections ---------------- */
   public static final int DRMAA_ERRNO_SUCCESS = 0; /* Routine returned normally with success. */
   public static final int DRMAA_ERRNO_INTERNAL_ERROR = 1; /* Unexpected or internal DRMAA error like memory allocation, system call failure, etc. */
   public static final int DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE = 2; /* Could not contact DRM system for this request. */
   public static final int DRMAA_ERRNO_AUTH_FAILURE = 3; /* The specified request is not processed successfully due to authorization failure. */
   public static final int DRMAA_ERRNO_INVALID_ARGUMENT = 4; /* The input value for an argument is invalid. */
   public static final int DRMAA_ERRNO_NO_ACTIVE_SESSION = 5; /* Exit routine failed because there is no active session */
   public static final int DRMAA_ERRNO_NO_MEMORY = 6; /* failed allocating memory */

   /* -------------- init and exit specific --------------- */
   public static final int DRMAA_ERRNO_INVALID_CONTACT_STRING = 7; /* Initialization failed due to invalid contact string. */
   public static final int DRMAA_ERRNO_DEFAULT_CONTACT_STRING_ERROR = 8; /* DRMAA could not use the default contact string to connect to DRM system. */
   public static final int DRMAA_ERRNO_NO_DEFAULT_CONTACT_STRING_SELECTED = 9; /* No default contact string was provided or selected. DRMAA requires that the default contact string is selected when there is more than one default contact string due to multiple DRMAA implementation contained in the binary module. */
   public static final int DRMAA_ERRNO_DRMS_INIT_FAILED = 10; /* Initialization failed due to failure to init DRM system. */
   public static final int DRMAA_ERRNO_ALREADY_ACTIVE_SESSION = 11; /* Initialization failed due to existing DRMAA session. */
   public static final int DRMAA_ERRNO_DRMS_EXIT_ERROR = 12; /* DRM system disengagement failed. */

   /* ---------------- job attributes specific -------------- */
   public static final int DRMAA_ERRNO_INVALID_ATTRIBUTE_FORMAT = 13; /* The format for the job attribute value is invalid. */
   public static final int DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE = 14; /* The value for the job attribute is invalid. */
   public static final int DRMAA_ERRNO_CONFLICTING_ATTRIBUTE_VALUES = 15; /* The value of this attribute is conflicting with a previously set attributes. */

   /* --------------------- job submission specific -------------- */
   public static final int DRMAA_ERRNO_TRY_LATER = 16; /* Could not pass job now to DRM system. A retry may succeed however (saturation). */
   public static final int DRMAA_ERRNO_DENIED_BY_DRM = 17; /* The DRM system rejected the job. The job will never be accepted due to DRM configuration or job template settings. */

   /* ------------------------------- job control specific ---------------- */
   public static final int DRMAA_ERRNO_INVALID_JOB = 18; /* The job specified by the 'jobid' does not exist. */
   public static final int DRMAA_ERRNO_RESUME_INCONSISTENT_STATE = 19; /* The job has not been suspended. The RESUME request will not be processed. */
   public static final int DRMAA_ERRNO_SUSPEND_INCONSISTENT_STATE = 20; /* The job has not been running, and it cannot be suspended. */
   public static final int DRMAA_ERRNO_HOLD_INCONSISTENT_STATE = 21; /* The job cannot be moved to a HOLD state. */
   public static final int DRMAA_ERRNO_RELEASE_INCONSISTENT_STATE = 22; /* The job is not in a HOLD state. */
   public static final int DRMAA_ERRNO_EXIT_TIMEOUT = 23; /* We have encountered a time-out condition for drmaa_synchronize or drmaa_wait. */
   public static final int DRMAA_ERRNO_NO_RUSAGE = 24; /* This error code is returned by drmaa_wait() when a job has finished but no rusage and stat data could be provided. */
   public static final int DRMAA_ERRNO_NO_MORE_ELEMENTS = 25; /* There are no more elements in the opaque string vector. */

   public static final int DRMAA_NO_ERRNO = 26;
}

/*
 * Agreed DRMAA job states as returned by drmaa_job_ps()
 */
public static interface DRMAA_PS {
 public static final int DRMAA_PS_UNDETERMINED = 0x00; /* process status cannot be determined */
 public static final int DRMAA_PS_QUEUED_ACTIVE = 0x10; /* job is queued and active */
 public static final int DRMAA_PS_SYSTEM_ON_HOLD = 0x11; /* job is queued and in system hold */
 public static final int DRMAA_PS_USER_ON_HOLD = 0x12; /* job is queued and in user hold */
 public static final int DRMAA_PS_USER_SYSTEM_ON_HOLD = 0x13; /* job is queued and in user and system hold */
 public static final int DRMAA_PS_RUNNING = 0x20; /* job is running */
 public static final int DRMAA_PS_SYSTEM_SUSPENDED = 0x21; /* job is system suspended */
 public static final int DRMAA_PS_USER_SUSPENDED = 0x22; /* job is user suspended */
 public static final int DRMAA_PS_USER_SYSTEM_SUSPENDED = 0x23; /* job is user and system suspended */
 public static final int DRMAA_PS_DONE = 0x30; /* job finished normally */
 public static final int DRMAA_PS_FAILED = 0x40;  /* job finished, but failed */
}

/*
 * Agreed DRMAA actions for drmaa_control()
 */
public static interface DRMAA_CONTROL {
 public static final int DRMAA_CONTROL_SUSPEND = 0;
 public static final int DRMAA_CONTROL_RESUME = 1;
 public static final int DRMAA_CONTROL_HOLD = 2;
 public static final int DRMAA_CONTROL_RELEASE = 3;
 public static final int DRMAA_CONTROL_TERMINATE = 4;
}

/* ------------------- Data types ------------------- */
/*
 * Agreed opaque DRMAA job template
 * struct drmaa_job_template_s is in japiP.h
 */
//typedef struct drmaa_job_template_s drmaa_job_template_t;

/* ---------- C/C++ language binding specific interfaces -------- */

//typedef struct drmaa_attr_names_s drmaa_attr_names_t;
//typedef struct drmaa_attr_values_s drmaa_attr_values_t;
//typedef struct drmaa_job_ids_s  drmaa_job_ids_t;

/*
 * get next string attribute from iterator
 *
 * returns DRMAA_ERRNO_SUCCESS or DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE
 * if no such exists
 */

public static native int drmaa_get_next_attr_name(/* drmaa_attr_names_t* */ Pointer values, Pointer value,
                             NativeLong value_len);
public static native int drmaa_get_next_attr_value(/* drmaa_attr_names_t* */ Pointer values, Pointer value,
                              NativeLong value_len);
public static native int drmaa_get_next_job_id(/* drmaa_job_ids_t* */ Pointer values, Pointer value,
                          NativeLong value_len);

/*
 * get element count of opaque string vector
 *
 * Gives the number of elements in the opaque string vector.  Useful for
 * copying the contents into an array.
 */
public static native int drmaa_get_num_attr_names(/* drmaa_attr_names_t* */ Pointer values, IntByReference size);
public static native int drmaa_get_num_attr_values(/* drmaa_attr_values_t* */ Pointer values, IntByReference size);
public static native int drmaa_get_num_job_ids(/* drmaa_job_ids_t* */ Pointer values, IntByReference size);

/*
 * release opaque string vector
 *
 * Opaque string vectors can be used without any constraint
 * until the release function has been called.
 */
public static native void drmaa_release_attr_names(/* drmaa_attr_names_t* */ Pointer values);
public static native void drmaa_release_attr_values(/* drmaa_attr_values_t* */ Pointer values);
public static native void drmaa_release_job_ids(/* drmaa_job_ids_t* */ Pointer values);

/* ------------------- init/exit routines ------------------- */
/*
 * Initialize DRMAA API library and create a new DRMAA Session. 'Contact'
 * is an implementation dependent string which MAY be used to specify
 * which DRM system to use. This routine MUST be called before any
 * other DRMAA calls, except for drmaa_version().
 * If 'contact' is NULL, the default DRM system SHALL be used provided there is
 * only one DRMAA implementation in the provided binary module.  When these is
 * more than one DRMAA implementation in the binary module, drmaa_init() SHALL
 * return the DRMAA_ERRNO_NO_DEFAULT_CONTACT_STRING_SELECTED error. drmaa_init()
 * SHOULD be called by only one of the threads. The main thread is RECOMMENDED.
 * A call by another thread SHALL return DRMAA_ERRNO_ALREADY_ACTIVE_SESSION.
 * When 'contact' is a a semi-colon separated list of name=value strings, the
 * strings will be parsed and interpreted.  The current list of accepted names
 * is:
 *    session -- the id of the session to which to reconnect
#if 0
 *    sge_root -- the SGE_ROOT to use
 *    sge_cell -- the SGE_CELL to use
#endif
 *
 * drmaa_init() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_INVALID_CONTACT_STRING,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_ALREADY_ACTIVE_SESSION,
 *    DRMAA_ERRNO_NO_DEFAULT_CONTACT_STRING_SELECTED, or
 *    DRMAA_ERRNO_DEFAULT_CONTACT_STRING_ERROR.
 */
public static native int drmaa_init(String contact, Pointer error_diagnosis, NativeLong error_diag_len);


/*
 * Disengage from DRMAA library and allow the DRMAA library to perform
 * any necessary internal clean up.
 * This routine SHALL end the current DRMAA Session, but SHALL NOT effect any
 * jobs (e.g., queued and running jobs SHALL remain queued and running).
 * drmaa_exit() SHOULD be called by only one of the threads. Other thread calls
 * to drmaa_exit() MAY fail since there is no active session.
 *
 * drmaa_exit() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_DRMS_EXIT_ERROR or
 *    DRMAA_ERRNO_NO_ACTIVE_SESSION.
 */
public static native int drmaa_exit(Pointer error_diagnosis, NativeLong error_diag_len);

/* ------------------- job template routines ------------------- */

/*
 * Allocate a new job template.
 *
 * drmaa_allocate_job_template() SHALL return DRMAA_ERRNO_SUCCESS on success,
 * otherwise:
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE,
 *    DRMAA_ERRNO_INTERNAL_ERROR or
 *    DRMAA_ERRNO_NO_MEMORY.
 */
public static native int drmaa_allocate_job_template(/* drmaa_job_template_t** */ PointerByReference jt, Pointer error_diagnosis, NativeLong error_diag_len);

/*
 * Deallocate a job template. This routine has no effect on jobs.
 *
 * drmaa_delete_job_template() SHALL return DRMAA_ERRNO_SUCCESS on success,
 * otherwise:
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE or
 *    DRMAA_ERRNO_INTERNAL_ERROR.
 */
public static native int drmaa_delete_job_template(/* drmaa_job_template_t* */ Pointer jt, Pointer error_diagnosis,
                              NativeLong error_diag_len);


/*
 * Adds ('name', 'value') pair to list of attributes in job template 'jt'.
 * Only non-vector attributes SHALL be passed.
 *
 * drmaa_set_attribute() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_INVALID_ATTRIBUTE_FORMAT,
 *    DRMAA_ERRNO_INVALID_ARGUMENT,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE or
 *    DRMAA_ERRNO_CONFLICTING_ATTRIBUTE_VALUES.
 */
public static native int drmaa_set_attribute(/* drmaa_job_template_t* */ Pointer jt, String name,
                        String value, Pointer error_diagnosis,
                        NativeLong error_diag_len);


/*
 * If 'name' is an existing non-vector attribute name in the job
 * template 'jt', then the value of 'name' SHALL be returned; otherwise,
 * NULL is returned.
 *
 * drmaa_get_attribute() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE.
 */
public static native int drmaa_get_attribute(/* drmaa_job_template_t* */ Pointer jt, String name, Pointer value,
                        NativeLong value_len, Pointer error_diagnosis,
                        NativeLong error_diag_len);

/* Adds ('name', 'values') pair to list of vector attributes in job template
 * 'jt'. Only vector attributes SHALL be passed.
 * A 'value' string vector containing n elements must be n+1 elements long, with
 * the nth value, i.e. value[n], being set to NULL as a delimitor.
 *
 * drmaa_set_vector_attribute() SHALL return DRMAA_ERRNO_SUCCESS on success,
 * otherwise:
 *    DRMAA_ERRNO_INVALID_ATTRIBUTE_FORMAT,
 *    DRMAA_ERRNO_INVALID_ARGUMENT,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_CONFLICTING_ATTRIBUTE_VALUES.
 */
public static native int drmaa_set_vector_attribute(/* drmaa_job_template_t* */ Pointer jt, String name,
                               Pointer value, Pointer error_diagnosis,
                               NativeLong error_diag_len);


/*
 * If 'name' is an existing vector attribute name in the job template 'jt',
 * then the values of 'name' are returned; otherwise, NULL is returned.
 *
 * drmaa_get_vector_attribute() SHALL return DRMAA_ERRNO_SUCCESS on success,
 * otherwise:
 *    DRMAA_ERRNO_INVALID_ATTRIBUTE_VALUE.
 */
public static native int drmaa_get_vector_attribute(/* drmaa_job_template_t* */ Pointer jt, String name,
                               /* drmaa_attr_values_t ** */ PointerByReference values,
                               Pointer error_diagnosis, NativeLong error_diag_len);


/*
 * SHALL return the set of supported attribute names whose associated
 * value type is String. This set SHALL include supported DRMAA reserved
 * attribute names and native attribute names.
 *
 * drmaa_get_attribute_names() SHALL return DRMAA_ERRNO_SUCCESS on success,
 * otherwise:
 *    DRMAA_ERRNO_NO_MEMORY.
 */
public static native int drmaa_get_attribute_names(/* drmaa_attr_names_t ** */ PointerByReference values,
                              Pointer error_diagnosis, NativeLong error_diag_len);

/*
 * SHALL return the set of supported attribute names whose associated
 * value type is String Vector.  This set SHALL include supported DRMAA reserved
 * attribute names and native attribute names.
 *
 * drmaa_get_vector_attribute_names() SHALL return DRMAA_ERRNO_SUCCESS on
 * success, otherwise:
 *    DRMAA_ERRNO_NO_MEMORY.
 */
public static native int drmaa_get_vector_attribute_names(/* drmaa_attr_names_t ** */ PointerByReference values,
                                     Pointer error_diagnosis,
                                     NativeLong error_diag_len);

/* ------------------- job submission routines ------------------- */

/*
 * Submit a job with attributes defined in the job template 'jt'.
 * The job identifier 'job_id' is a printable, NULL terminated string,
 * identical to that returned by the underlying DRM system.
 *
 * drmaa_run_job() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_TRY_LATER,
 *    DRMAA_ERRNO_DENIED_BY_DRM,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE or
 *    DRMAA_ERRNO_AUTH_FAILURE.
 */
public static native int drmaa_run_job(Pointer job_id, NativeLong job_id_len,
                  /* drmaa_job_template_t * */ Pointer jt, Pointer error_diagnosis,
                  NativeLong error_diag_len);

/*
 * Submit a set of parametric jobs, dependent on the implied loop index, each
 * with attributes defined in the job template 'jt'.
 * The job identifiers 'job_ids' SHALL all be printable,
 * NULL terminated strings, identical to those returned by the underlying
 * DRM system. Nonnegative loop bounds SHALL NOT use file names
 * that start with minus sign like command line options.
 * DRMAA defines a special index placeholder, drmaa_incr_ph, (which has the
 * value "$incr_pl$") that is used to construct parametric job templates.
 * For example:
 * //C++ string syntax used
 * drmaa_set_attribute(pjt, "stderr", drmaa_incr_ph + ".err" );
 *
 * drmaa_run_bulk_jobs() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_TRY_LATER,
 *    DRMAA_ERRNO_DENIED_BY_DRM,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE or
 *    DRMAA_ERRNO_AUTH_FAILURE.
 */
public static native int drmaa_run_bulk_jobs(/* drmaa_job_ids_t ** */ PointerByReference jobids,
                        /* drmaa_job_template_t * */ Pointer jt, int start, int end,
                        int incr, Pointer error_diagnosis, NativeLong error_diag_len);

/* ------------------- job control routines ------------------- */

/*
 * Start, stop, restart, or kill the job identified by 'job_id'.
 * If 'job_id' is DRMAA_JOB_IDS_SESSION_ALL then this routine
 * acts on all jobs *submitted* during this DRMAA session.
 * The legal values for 'action' and their meanings SHALL be:
 * DRMAA_CONTROL_SUSPEND:     stop the job,
 * DRMAA_CONTROL_RESUME:      (re)start the job,
 * DRMAA_CONTROL_HOLD:        put the job on-hold,
 * DRMAA_CONTROL_RELEASE:     release the hold on the job, and
 * DRMAA_CONTROL_TERMINATE:   kill the job.
 *
 * This routine SHALL return once the action has been acknowledged by
 * the DRM system, but does not necessarily wait until the action
 * has been completed.
 *
 * drmaa_control() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE,
 *    DRMAA_ERRNO_AUTH_FAILURE,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_RESUME_INCONSISTENT_STATE,
 *    DRMAA_ERRNO_SUSPEND_INCONSISTENT_STATE,
 *    DRMAA_ERRNO_HOLD_INCONSISTENT_STATE,
 *    DRMAA_ERRNO_RELEASE_INCONSISTENT_STATE or
 *    DRMAA_ERRNO_INVALID_JOB.
 */
public static native int drmaa_control(String jobid, int action, Pointer error_diagnosis,
                  NativeLong error_diag_len);


/*
 * Wait until all jobs specified by 'job_ids' have finished
 * execution. If 'job_ids' is DRMAA_JOB_IDS_SESSION_ALL then this routine
 * waits for all jobs *submitted* during this DRMAA session. The timeout value
 * is used to specify the number of seconds to wait for the job to fail finish
 * before returning if a result is not immediately available.  The value
 * DRMAA_TIMEOUT_WAIT_FOREVER can be used to specify that routine should wait
 * indefinitely for a result. The value DRMAA_TIMEOUT_NO_WAIT can be used to
 * specify that the routine should return immediately if no result is available.
 * If the call exits before timeout, all the jobs have
 * been waited on or there was an interrupt.
 * If the invocation exits on timeout, the return code is
 * DRMAA_ERRNO_EXIT_TIMEOUT. The caller SHOULD check system time before and
 * after this call in order to check how much time has passed.
 *
 * The dispose parameter specifies how to treat reaping information:
 * True=1      "fake reap", i.e. dispose of the rusage data
 * False=0     do not reap
 *
 * A 'job_ids' string vector containing n elements must be n+1 elements long,
 * with the nth value, i.e. job_ids[n], being set to NULL as a delimitor.
 *
 * drmaa_synchronize() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE,
 *    DRMAA_ERRNO_AUTH_FAILURE,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_EXIT_TIMEOUT or
 *    DRMAA_ERRNO_INVALID_JOB.
 */
public static native int drmaa_synchronize(Pointer job_ids, NativeLong timeout, int dispose,
                      Pointer error_diagnosis, NativeLong error_diag_len);


/*
 * This routine SHALL wait for a job with job_id to fail or finish execution. If
 * the special string, DRMAA_JOB_IDS_SESSION_ANY is provided as the job_id,
 * this routine SHALL wait for any job from the session. This routine is modeled
 * on the wait3 POSIX routine. The timeout value is used to specify the number
 * of seconds to wait for the job to fail finish before returning if a result is
 * not immediately available.  The value DRMAA_TIMEOUT_WAIT_FOREVER can be
 * used to specify that routine should wait indefinitely for a result. The value
 * DRMAA_TIMEOUT_NO_WAIT may be specified that the routine should return
 * immediately if no result is available.
 * If the call exits before timeout ,the job has been waited on
 * successfully or there was an interrupt.
 * If the invocation exits on timeout, the return code is
 * DRMAA_ERRNO_EXIT_TIMEOUT. The caller SHOULD check system time before and
 * after this call in order to check how much time has passed.
 * The routine reaps jobs on a successful call, so any subsequent calls
 * to drmaa_wait SHOULD fail returning an error DRMAA_ERRNO_INVALID_JOB meaning
 * that the job has been already reaped. This error is the same as if the job
 * was unknown. Failing due to an elapsed timeout has an effect that it is
 * possible to issue drmaa_wait multiple times for the same job_id.  When
 * successful, the rusage information SHALL be provided as an array of strings,
 * where each string complies with the format <name>=<value>. The string portion
 * <value> contains the amount of resources consumed by the job and is opaque.
 * The 'stat' drmaa_wait parameter is used in the drmaa_w* functions for
 * providing more detailed information about job termination if available. An
 * analogous set of macros is defined in POSIX for analyzing the wait3(2) OUT
 * parameter 'stat'.
 *
 * drmaa_wait() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE,
 *    DRMAA_ERRNO_AUTH_FAILURE,
 *    DRMAA_ERRNO_NO_RUSAGE,
 *    DRMAA_ERRNO_NO_MEMORY,
 *    DRMAA_ERRNO_EXIT_TIMEOUT or
 *    DRMAA_ERRNO_INVALID_JOB.
 */
public static native int drmaa_wait(String job_id, Pointer job_id_out, NativeLong job_id_out_len,
               IntByReference stat, NativeLong timeout, /* drmaa_attr_values_t ** */ PointerByReference rusage,
               Pointer error_diagnosis, NativeLong error_diag_len);

/*
 * Evaluates into 'exited' a non-zero value if stat was returned for a
 * job that terminated normally. A zero value can also indicate that
 * altough the job has terminated normally an exit status is not available
 * or that it is not known whether the job terminated normally. In both
 * cases drmaa_wexitstatus() SHALL NOT provide exit status information.
 * A non-zero 'exited' value indicates more detailed diagnosis can be provided
 * by means of drmaa_wifsignaled(), drmaa_wtermsig() and drmaa_wcoredump().
 */
public static native int drmaa_wifexited(IntByReference exited, int stat, Pointer error_diagnosis,
                    NativeLong error_diag_len);

/*
 * If the OUT parameter 'exited' of drmaa_wifexited() is non-zero,
 * this function evaluates into 'exit_code' the exit code that the
 * job passed to _exit() (see exit(2)) or exit(3C), or the value that
 * the child process returned from main.
 */
public static native int drmaa_wexitstatus(IntByReference exit_status, int stat, Pointer error_diagnosis,
                      NativeLong error_diag_len);

/*
 * Evaluates into 'signaled' a non-zero value if status was returned
 * for a job that terminated due to the receipt of a signal. A zero value
 * can also indicate that altough the job has terminated due to the receipt
 * of a signal the signal is not available or that it is not known whether
 * the job terminated due to the receipt of a signal. In both cases
 * drmaa_wtermsig() SHALL NOT provide signal information.
 */
public static native int drmaa_wifsignaled(IntByReference signaled, int stat, Pointer error_diagnosis,
                      NativeLong error_diag_len);

/*
 * If the OUT parameter 'signaled' of drmaa_wifsignaled(stat) is
 * non-zero, this function evaluates into signal a string representation of the
 * signal that caused the termination of the job. For signals declared by POSIX,
 * the symbolic names SHALL be returned (e.g., SIGABRT, SIGALRM).
 * For signals not declared by POSIX, any other string MAY be returned.
 */
public static native int drmaa_wtermsig(Pointer signal, NativeLong signal_len, int stat,
                   Pointer error_diagnosis, NativeLong error_diag_len);

/*
 * If the OUT parameter 'signaled' of drmaa_wifsignaled(stat) is
 * non-zero, this function evaluates into 'core_dumped' a non-zero value
 * if a core image of the terminated job was created.
 */
public static native int drmaa_wcoredump(IntByReference core_dumped, int stat, Pointer error_diagnosis,
                    NativeLong error_diag_len);

/*
 * Evaluates into 'aborted' a non-zero value if 'stat'
 * was returned for a job that ended before entering the running state.
 */
public static native int drmaa_wifaborted(IntByReference aborted, int stat, Pointer error_diagnosis,
                     NativeLong error_diag_len);



/*
 * Get the program status of the job identified by 'job_id'.
 * The possible values returned in 'remote_ps' and their meanings SHALL be:
 *
 * DRMAA_PS_UNDETERMINED          = 0x00: process status cannot be determined
 * DRMAA_PS_QUEUED_ACTIVE         = 0x10: job is queued and active
 * DRMAA_PS_SYSTEM_ON_HOLD        = 0x11: job is queued and in system hold
 * DRMAA_PS_USER_ON_HOLD          = 0x12: job is queued and in user hold
 * DRMAA_PS_USER_SYSTEM_ON_HOLD   = 0x13: job is queued and in user and system
 *                                        hold
 * DRMAA_PS_RUNNING               = 0x20: job is running
 * DRMAA_PS_SYSTEM_SUSPENDED      = 0x21: job is system suspended
 * DRMAA_PS_USER_SUSPENDED        = 0x22: job is user suspended
 * DRMAA_PS_USER_SYSTEM_SUSPENDED = 0x23: job is user and system suspended
 * DRMAA_PS_DONE                  = 0x30: job finished normally
 * DRMAA_PS_FAILED                = 0x40: job finished, but failed
 *
 * DRMAA SHOULD always get the status of job_id from DRM system, unless the
 * previous status has been DRMAA_PS_FAILED or DRMAA_PS_DONE and the status has
 * been successfully cached. Terminated jobs get DRMAA_PS_FAILED status.
 *
 * drmaa_synchronize() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_DRM_COMMUNICATION_FAILURE,
 *    DRMAA_ERRNO_AUTH_FAILURE,
 *    DRMAA_ERRNO_NO_MEMORY or
 *    DRMAA_ERRNO_INVALID_JOB.
 */
public static native int drmaa_job_ps(String job_id, IntByReference remote_ps, Pointer error_diagnosis,
                 NativeLong error_diag_len);

/* ------------------- auxiliary routines ------------------- */

/*
 * SHALL return the error message text associated with the errno number. The
 * routine SHALL return null string if called with invalid ERRNO number.
 */
public static native String drmaa_strerror(int drmaa_errno);

/*
 * If called before drmaa_init(), it SHALL return a comma delimited default
 * DRMAA implementation contacts string, one per each DRM system provided
 * implementation. If called after drmaa_init(), it SHALL return the selected
 * contact string. The output string is Implementation dependent.
 * drmaa_get_contact() SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_INTERNAL_ERROR.
 */
public static native int drmaa_get_contact(Pointer contact, NativeLong contact_len,
         Pointer error_diagnosis, NativeLong error_diag_len);

/*
 * OUT major - major version number (non-negative integer)
 * OUT minor - minor version number (non-negative integer)
 * SHALL return the major and minor version numbers of the DRMAA library;
 * for DRMAA 1.0, 'major' is 1 and 'minor' is 0.
 */
public static native int drmaa_version(IntByReference major, IntByReference minor,
         Pointer error_diagnosis, NativeLong error_diag_len);


/*
 * If called before drmaa_init(), it SHALL return a comma delimited DRM systems
 * string, one per each DRM system provided implementation. If called after
 * drmaa_init(), it SHALL return the selected DRM system. The output string is
 * implementation dependent.
 *
 * drmaa_get_DRM_system() SHALL return DRMAA_ERRNO_SUCCESS on success,
 * otherwise:
 *    DRMAA_ERRNO_INTERNAL_ERROR.
 */
public static native int drmaa_get_DRM_system(Pointer drm_system, NativeLong drm_system_len,
         Pointer error_diagnosis, NativeLong error_diag_len);


/*
 * If called before drmaa_init(), it SHALL return a comma delimited DRMAA
 * implementations string, one per each DRM system provided implementation. If
 * called after drmaa_init(), it SHALL return the selected DRMAA implementation.
 * The output (string) is implementation dependent. drmaa_get_DRM_implementation
 * routine SHALL return DRMAA_ERRNO_SUCCESS on success, otherwise:
 *    DRMAA_ERRNO_INTERNAL_ERROR.
 */
public static native int drmaa_get_DRMAA_implementation(Pointer drmaa_impl, NativeLong drmaa_impl_len,
         Pointer error_diagnosis, NativeLong error_diag_len);
}

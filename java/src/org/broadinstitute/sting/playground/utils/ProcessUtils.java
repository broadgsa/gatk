package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.xReadLines;
import org.broadinstitute.sting.utils.StingException;
import org.apache.log4j.Logger;

/**
 * A set of utilities for managing external processes.
 */
public class ProcessUtils {
	private static Logger logger = Logger.getLogger(ProcessUtils.class);

	/**
	 * Runs a command line and returns the result code.
	 * @param command Command line to execute.
	 * @return Result code of the command.
	 */
	public static int runCommandAndWait(String command) {
		try {
			logger.debug("Running command: " + command);

			Process p = Runtime.getRuntime().exec(command);
			int result = p.waitFor();

			if (logger.isDebugEnabled()) {
				for (String line : new xReadLines(p.getInputStream())) {
					logger.debug("command: " + line);
				}
				for (String line : new xReadLines(p.getErrorStream())) {
					logger.error("command: " + line);
				}
			}

			logger.debug("Command exited with result: " + result);

			return result;
		} catch (Exception e) {
			throw new StingException("Error running command:" + command, e);
		}
	}
}

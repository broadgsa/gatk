/*
 * Copyright (c) 2010 The Broad Institute
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

package org.broadinstitute.sting.playground.utils;

import org.broadinstitute.sting.utils.GATKException;
import org.broadinstitute.sting.utils.text.XReadLines;
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
				for (String line : new XReadLines(p.getInputStream())) {
					logger.debug("command: " + line);
				}
				for (String line : new XReadLines(p.getErrorStream())) {
					logger.error("command: " + line);
				}
			}

			logger.debug("Command exited with result: " + result);

			return result;
		} catch (Exception e) {
			throw new GATKException("Error running command:" + command, e);
		}
	}
}

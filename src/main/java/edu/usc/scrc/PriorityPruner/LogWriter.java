/**
Copyright (c) 2014 Christopher K. Edlund, Malin Anker, Fredrick R. Schumacher, W. James Gauderman, David V. Conti,
University of Southern California,
Los Angeles, CA  90033, USA.

The MIT License (MIT)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */

package edu.usc.scrc.PriorityPruner;

import org.apache.log4j.Logger;

/**
 * This class handles the logging of PriorityPruner by using the Apache Commons
 * Log4j library. The class is currently implementing the Singleton-pattern.
 */
public class LogWriter {

	private static LogWriter singletonObject = null;
	private static Logger logger = Logger.getLogger(LogWriter.class);

	/**
	 * Private constructor for LogWriter, to implement the Singleton-pattern.
	 */
	private LogWriter() {
	}

	/**
	 * To implement the Singleton-pattern, the following getInstance-method
	 * manages the public access to this class.
	 * 
	 * @return instance of this class
	 */
	public static LogWriter getInstance() {

		if (singletonObject == null) {
			singletonObject = new LogWriter();
			return singletonObject;
		} else {
			return singletonObject;
		}
	}

	// public getters and setter for private fields of this class

	public static Logger getLogger() {
		return logger;
	}

	public static void setLogger(Logger logger) {
		LogWriter.logger = logger;
	}
}
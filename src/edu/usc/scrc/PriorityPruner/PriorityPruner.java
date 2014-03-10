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

import java.io.IOException;
import java.util.*;

import org.apache.log4j.FileAppender;
import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PatternLayout;

/**
 * This is the main class of PriorityPruner. It initiates crucial functions of
 * this program such as parsing of the command line options and the creation of
 * object responsible for the pruning of provided SNPs. At the moment, the
 * tplink-format is the only supported input format for genotype and family
 * data. A SNP input file is also required for parsing of the tped-format.
 */
public class PriorityPruner {
	public static void main(String[] args) {
		try {
			CommandLineOptions options = null;
			try {
				// parses the command line options
				options = CommandLineOptions.getInstance(args);

				// in case a problem is encountered during parsing, a default
				// log file is created with this information. This is because
				// the -out option, defining name and path of output files,
				// hasen't been parsed if an exception is caught.
			} catch (PriorityPrunerException e) {
				// defines logging properties
				PatternLayout layout = new PatternLayout(
						"org.apache.log4j.PatternLayout");
				layout.setConversionPattern("%m%n");
				FileAppender logFileAppender = new FileAppender(layout,
						"PriorityPruner.log", false);
				Logger.getRootLogger().addAppender(logFileAppender);
				// TODO: update download info
				LogWriter
						.getLogger()
						.info("Welcome to PriorityPruner version 1.0 \n\nFor latest version please visit: http://sourceforge.net/projects/prioritypruner/\n(C) 2014 University of Southern California, The MIT License (MIT)\n-------------------------------------------------------------------------------\n");
				LogWriter.getLogger().warn(e.getMessage());
				System.exit(0);
			}

			// defines logging properties
			PatternLayout layout = new PatternLayout(
					"org.apache.log4j.PatternLayout");
			layout.setConversionPattern("%m%n");
			FileAppender logFileAppender = new FileAppender(layout,
					options.getOutputFile() + ".log", false);
			Logger.getRootLogger().addAppender(logFileAppender);
			// TODO: update download info
			LogWriter
					.getLogger()
					.info("Welcome to PriorityPruner version 1.0 \n\nFor latest version please visit: http://sourceforge.net/projects/prioritypruner/\n(C) 2014 University of Southern California, The MIT License (MIT)\n-------------------------------------------------------------------------------\n");
			LogWriter.getLogger().info(options.getOptionsInEffect());
			long start = System.currentTimeMillis();
			LogWriter.getLogger()
					.info("Analysis started: " + new Date() + "\n");
			LogWriter.getLogger().info(
					"Writing this text to log file [ "
							+ options.getOutputFile() + ".log ]");
			// sets logging level to debug if the -verbose option is specified
			// in command line
			if (options.getVerbose()) {
				Logger.getRootLogger().setLevel(Level.DEBUG);
			}

			// checks that required input files are provided, else program will
			// terminate. At the moment only the tplink-format is supported.
			// When other formats are added a more generic check like
			// if(Genotypes.getSnpGenotypes()!=null &&
			// Genotypes.getSnpGenotypes()!=null) could be used instead.
			if (options.getTped() != null && options.getTfam() != null
					&& options.getSnpFilePath() != null) {
				// initiates pruning
				new Pruner();
				LogWriter.getLogger()
						.info("\nAnalysis finished: " + new Date());
				long end = System.currentTimeMillis();
				// prints duration time
				printDuration(end - start);
			} else {
				throw new PriorityPrunerException(
						"Missing genotype or family data, please rerun program with correct input files using following options: \n-tfile : specifies the tped and tfam files, if they have the same prefix (entered without suffix) \n-tped : specifies the tped file (entered with suffix: \".tped\") \n-tfam : specifies the tfam file (entered with suffix: \".tfam\") \n-snpfile : specifies the SNP input file (entered with suffix: \".txt\")");
			}
			// catches all PriorityPrunerExceptions that get thrown during
			// execution of this program and exits it
		} catch (PriorityPrunerException e) {
			LogWriter.getLogger().warn(e.getMessage());
			System.exit(0);
		} catch (IOException e) {
			LogWriter
					.getLogger()
					.warn("Could not create file: "
							+ e.getMessage()
							+ "\nPlease check that correct file path is provided.");
			System.exit(0);
		}
	}

	/**
	 * This method prints the duration time in correct format.
	 * 
	 * @param diff
	 *            duration in milliseconds
	 */
	private static void printDuration(long diff) {

		int days = 0;
		int hours = 0;
		int minutes = 0;
		int seconds = 0;
		int millis = 0;

		// if diff is to large to be stored as a long and wraps around
		if (diff < 0) {
			LogWriter.getLogger().info("Duration time not available.");
			return;
		}
		if (diff >= 86400000) {
			days = (int) diff / 86400000;
			diff = diff % 86400000;
		}
		if (diff >= 3600000) {
			hours = (int) diff / 3600000;
			diff = diff % 3600000;
		}
		if (diff >= 60000) {
			minutes = (int) diff / 60000;
			diff = diff % 60000;
		}
		if (diff >= 1000) {
			seconds = (int) diff / 1000;
			diff = diff % 1000;
		}
		millis = (int) diff;

		if (days != 0) {
			LogWriter.getLogger().info(
					"Duration: " + days + " day(s), " + hours + " hour(s), "
							+ minutes + " minute(s), " + seconds
							+ " second(s), " + millis + " millisecond(s) ");
		} else if (hours != 0) {
			LogWriter.getLogger().info(
					"Duration: " + hours + " hour(s), " + minutes
							+ " minute(s), " + seconds + " second(s), "
							+ millis + " millisecond(s) ");
		} else if (minutes != 0) {
			LogWriter.getLogger().info(
					"Duration: " + minutes + " minute(s), " + seconds
							+ " second(s), " + millis + " millisecond(s) ");
		} else if (seconds != 0) {
			LogWriter.getLogger().info(
					"Duration: " + seconds + " second(s), " + millis
							+ " millisecond(s) ");
		} else {
			LogWriter.getLogger().info(
					"Duration: " + millis + " millisecond(s) ");
		}
	}
}
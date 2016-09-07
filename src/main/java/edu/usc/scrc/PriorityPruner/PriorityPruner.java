/**
Copyright (c) 2016 Christopher K. Edlund, Malin Anker, Fredrick R. Schumacher, W. James Gauderman, David V. Conti,
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
		
		// LD output file
		LinkageDisequilibriumFile ldFile = null;
		int returnCode = 0;
		FileAppender logFileAppender = null;
		
		try {
			CommandLineOptions options = null;
			try {
				// parses the command line options
				options = new CommandLineOptions(args);

				// in case a problem is encountered during parsing, a default
				// log file is created with this information. This is because
				// the -out option, defining name and path of output files,
				// hasen't been parsed if an exception is caught.
			} catch (PriorityPrunerException e) {
				// defines logging properties
				PatternLayout layout = new PatternLayout("org.apache.log4j.PatternLayout");
				layout.setConversionPattern("%m%n");
				FileAppender defaultLogFileAppender = new FileAppender(layout,"PriorityPruner.log", false);
				Logger.getRootLogger().addAppender(defaultLogFileAppender);
				// TODO: update download info
				LogWriter.getLogger()
						.info("Welcome to PriorityPruner version 0.1.4 \n\nFor latest version please visit: http://prioritypruner.sourceforge.net\n(C) 2016 Christopher K. Edlund et al., The MIT License (MIT)\n-------------------------------------------------------------------------------\n");
				LogWriter.getLogger().warn("ERROR: " + e.getMessage());	
				defaultLogFileAppender.close();
				System.exit(1);
			}

			// defines logging properties
			PatternLayout layout = new PatternLayout(
					"org.apache.log4j.PatternLayout");
			layout.setConversionPattern("%m%n");
			logFileAppender = new FileAppender(layout,
					options.getOutputPrefix() + ".log", false);
			Logger.getRootLogger().addAppender(logFileAppender);
			// TODO: update download info
			LogWriter
					.getLogger()
					.info("Welcome to PriorityPruner version 0.1.4 \n\nFor latest version please visit: http://prioritypruner.sourceforge.net\n(C) 2016 Christopher K. Edlund et al., The MIT License (MIT)\n-------------------------------------------------------------------------------\n");
			
			long start = System.currentTimeMillis();
			if (options.getOutputPrefix()!= null){
			LogWriter.getLogger().info(
					"Writing this text to log file [ "
							+ options.getOutputPrefix() + ".log ]");
			}
			
			LogWriter.getLogger()
					.info("Analysis started: " + new Date() + "\n");

			LogWriter.getLogger().info(options.getOptionsInEffect());
			
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
					&& options.getSnpTablePath() != null) {
				
				// parse the list of SNPs to prune; all pruning results are 
				// stored in this object
				SnpListFile snpListFile = new SnpListFile(options.getSnpTablePath(), options.getNumMetrics(), options);
				
				// create an LD file if specified by user
				if (options.isOutputLDTable()){
					ldFile = new LinkageDisequilibriumFile(options);
				}
				
				// parse keep/remove samples list in case --keep or --remove is specified by user
				PlinkSampleListFile keepRemoveSamples = null;
				if (options.getKeep() != null){
					keepRemoveSamples = new PlinkSampleListFile(options.getKeep());
				}else if (options.getRemove() != null){
					keepRemoveSamples = new PlinkSampleListFile(options.getRemove());
				}
				
				// parse genotypes (at this time only Transposed PLINK is supported)
				Genotypes genotypes = new TPlink(options.getTped(), options.getTfam(), snpListFile, keepRemoveSamples, options);
				
				// verify all SNPs from snpListFile are in genotypes
				checkSnpsAreInGenotypeFile(snpListFile);
						
				// prune the list of SNPs
				Pruner pruner = new Pruner(genotypes, snpListFile, ldFile, options);
				
				// write the results file
				new ResultsFile(snpListFile, options);
				
				LogWriter.getLogger().info("\nAnalysis finished: " + new Date());
				long end = System.currentTimeMillis();
				// prints duration time
				printDuration(end - start);
			} else {
				throw new PriorityPrunerException(
						"Missing genotype dataset or SNP input table. You must define input files with the following options: \n-tfile : specifies the tped and tfam files, if they have the same prefix (entered without suffix) \n-tped : specifies the tped file (entered with suffix: \".tped\") \n-tfam : specifies the tfam file (entered with suffix: \".tfam\") \n-snp_table : specifies the SNP input table");
			}
			// catches all PriorityPrunerExceptions that get thrown during
			// execution of this program and exits it
		} catch (PriorityPrunerException e) {
			LogWriter.getLogger().warn("\r\nERROR: " + e.getMessage());
			returnCode = 1;
		} catch (IOException e) {
			LogWriter
					.getLogger()
					.warn("Could not create file: "
							+ e.getMessage()
							+ "\nPlease check that correct file path is provided.");
			returnCode = 1;
		} finally{
			// close all open files
			if (ldFile != null){
				try{
					ldFile.close();
				}catch(Exception e){
					returnCode = 1;
					e.printStackTrace();
				}
			}
			if (logFileAppender != null){
				logFileAppender.close();
				Logger.getRootLogger().removeAppender(logFileAppender);
			}
			
			//exit
			System.exit(returnCode);
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
	
	/**
	 * Checks that all SNPs defined in the SNP Input Table have corresponding data in the
	 * genotype dataset.
	 * 
	 * @throws PriorityPrunerException
	 *             if SNP not found in the genotype dataset
	 */
	private static void checkSnpsAreInGenotypeFile(SnpListFile snpListFile) throws PriorityPrunerException{
		// checks that corresponding SNP is provided in genotypes file
		for (SnpInfo snp : snpListFile.getSnps()) {
			if (!snp.getInTped()) {
				throw new PriorityPrunerException(snp.getSnpName()
						+ " at chromsome " + snp.getChr()
						+ ":" + snp.getPos()
						+ " with alleles: " + snp.getAllele1() + " "
						+ snp.getAllele2() + ", not found in genotype dataset.");
			}
		}
	}

	
}
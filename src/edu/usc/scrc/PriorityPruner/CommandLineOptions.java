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



//import javax.swing.text.html.HTMLDocument.Iterator;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.cli.AlreadySelectedException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.MissingArgumentException;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.UnrecognizedOptionException;

//import org.apache.commons.cli.*;

/**
 * This class defines and parses the command line options given as input to the
 * program. It uses the Apache Commons CLI 1.2 library for both specifying and
 * parsing the options.
 */
public class CommandLineOptions {

	// parameters storing command line options. In cases where it apply, they
	// are set to default values.
	private static CommandLineOptions singletonObject = null;
	private long maxDistance = 500000;
	private double minMaf = 0;
	//private double minHwe = 0;
	private double minSnpCallRate = 0.1;
	private double minDesignScore = 0;
	//private double absoluteMinDesignScore = 0;
	private ArrayList<Metric> metricWeights = new ArrayList<Metric>();
	private ArrayList<SurrogateThreshold> sortedSurrogateThresholds = new ArrayList<SurrogateThreshold>();
	private ArrayList<R2Threshold> sortedR2Thresholds = new ArrayList<R2Threshold>();
	private String chr = null;
	private String outputPrefix = null;
	//private String forceIncludeFilePath = null;
	private String optionsInEffect = null;
	private boolean sortByForceIncludeAndPValue = true;
	private String snpTablePath = null;
	private int numMetrics = 0;
	private String tped = null;
	private String tfam = null;
	private String tfile = null;
	private Option[] parsedOptions = null;
	private boolean verbose = false;
	private boolean outputLDTable = false;
	private String remove = null;
	private String keep = null;
	private double keepPercentage = -1;
	// flag indicating if surrogates should be added for force included SNPs
	// this value goes together with option "nos", [no s(urrogates)]
	private boolean addSurrogatesForForceIncludedSnps = true;
	// flag indicating if an additional surrogate should be added in case index
	// SNP doesn't pass the design score minimum
	//private boolean useSurrogateForNonPassingIndexSnp = false;

	/**
	 * Private constructor that initiates and parses command line options.
	 * 
	 * @param args
	 *            command line arguments
	 * @throws PriorityPrunerException
	 *             if errors are encountered during parsing
	 */
	private CommandLineOptions(String[] args) throws PriorityPrunerException {
		Options helpOptions = initHelpOptions();
		Options options = initOptions();
		// if help option is encountered
		if (helpParse(helpOptions, options, args)) {
			System.exit(0);
		} else {
			parse(options, args);
		}
	}

	/**
	 * To implement the Singleton-pattern, the following two getInstance-methods
	 * manages the public access to this class. This method calls the private
	 * constructor to create a CommandLineOption-object if it's not already
	 * created.
	 * 
	 * @param args
	 *            command line arguments
	 * @return instance of this class
	 * @throws PriorityPrunerException
	 *             if errors occur during parsing
	 */
	public static CommandLineOptions getInstance(String[] args)
			throws PriorityPrunerException {

		if (singletonObject == null) {
			singletonObject = new CommandLineOptions(args);
			return singletonObject;
		} else {
			return singletonObject;
		}
	}

	/**
	 * To implement the Singleton-pattern, this and the previous
	 * getInstance-method manages the public access to this class. This method
	 * can only return an instance of CommandLineOptions if previous
	 * getInstance-method has been called first.
	 * 
	 * @return instance of this class
	 * @throws PriorityPrunerException
	 *             if an instance isn't already initialized
	 */
	public static CommandLineOptions getInstance()
			throws PriorityPrunerException {

		if (singletonObject == null) {
			// TODO: What kind of message should go here? This should hopefully
			// never happen?
			throw new PriorityPrunerException(
					"An instance of CommandLineOptions is not yet initialized.");
		} else {
			return singletonObject;
		}
	}

	/**
	 * Creates an Option-object with both a full and a shorthand version of its
	 * name.
	 * 
	 * @param numArgs
	 *            number of arguments required for every argument entry
	 * @param typeOfArg
	 *            type of argument expected
	 * @param longOpt
	 *            full name of this option
	 * @param desc
	 *            description of this option
	 * @param isReq
	 *            specifies if this option is required or not
	 * @param shortOpt
	 *            shorthand name for this option
	 * @return a static Option-object
	 */
	private static Option createOptionTwoNames(int numArgs, String typeOfArg,
			String longOpt, String desc, boolean isReq, String shortOpt) {

		OptionBuilder.hasArgs(numArgs);
		OptionBuilder.withArgName(typeOfArg);
		OptionBuilder.withLongOpt(longOpt);
		OptionBuilder.withDescription(desc);
		OptionBuilder.isRequired(isReq);
		return OptionBuilder.create(shortOpt);
	}

	/**
	 * Creates an Option-object with specified parameters, without a shorthand
	 * version of its name.
	 * 
	 * @param numArgs
	 *            number of arguments required for an argument entry
	 * @param typeOfArg
	 *            type of argument expected
	 * @param desc
	 *            description of this option
	 * @param isReq
	 *            specifies if this option is required or not
	 * @param shortOpt
	 *            short name for this option
	 * @return a static Option-object
	 */
	private static Option createOptionOneName(int numArgs, String typeOfArg,
			String desc, boolean isReq, String longOpt) {

		OptionBuilder.hasArgs(numArgs);
		OptionBuilder.withArgName(typeOfArg);
		OptionBuilder.withDescription(desc);
		OptionBuilder.isRequired(isReq);
		return OptionBuilder.create(longOpt);
	}

	/**
	 * Initiates help option and adds it to a collection.
	 * 
	 * @return an Options-object, containing the help option
	 */
	private Options initHelpOptions() {
		Options helpOptions = new Options();
		Option help = createOptionTwoNames(0, "none", "help", "Print help",
				false, "h");
		helpOptions.addOption(help);
		return helpOptions;
	}

	/**
	 * Initiates options with input arguments from user, and adds them to a
	 * collection.
	 * 
	 * @return an Options-object, which is a collection of all Option-objects
	 */
	private Options initOptions() {
		// when defining required options/changing descriptions, notice the
		// temporary solution for printing out missing required options (when a
		// MissingOptionException is caught)
		
		Options options = new Options();
		
		//max-distance
		Option maxDistance = createOptionOneName(
				1,
				"integer",
				"Max distance in base pairs between SNPs to calculate pair-wise LD. Default is 500,000.",
				false, "max_distance");
		
		//min-maf
		Option minMaf = createOptionOneName(1, "float",
				"Minimum minor allele frequency (MAF) for a SNP to be selected. Does not apply to force-selected SNPs. Default is 0.", false, "min_maf");
		
		//min-hwe
//		Option minHwPvalue = createOptionOneName(1, "float",
//				"Minimum Hardy-Weinberg p-value for a SNP to be selected. Does not apply to force-selected SNPs. Default is 0.", false, "min_hwe");
		
		//min-snp-callrate
		Option minSnpCallRate = createOptionOneName(1, "float",
				"Minimum SNP call rate for a SNP to be selected. Does not apply to force-selected SNPs. Default is 0.", false, "min_snp_callrate");
		
		//min-design-score
		Option minDesignScore = createOptionOneName(1, "float",
				"Minimum design score for a SNP to be selected. Does not apply to force-selected SNPs. Default is 0.", false, "min_design_score");
		
//		Option absMinDesignScore = createOptionTwoNames(1,
//				"0.0 - no upper limit", "absmindesignscore",
//				"Set absolute minimum design score value allowed", false,
//				"amds");
		
		//metric
		Option metric = createOptionOneName(2,
				"name, weight",
				"Specify name and weight of a metric.", false, "metric");
		
		//st
		Option surrogateThreshold = createOptionOneName(
				2,
				"p-value, num_surrogates",
				"Add surrogate-threshold, specify both p-value and the number of surrogates",
				false, "st");
		
		//r2t
		Option r2Threshold = createOptionOneName(
				2,
				"p-value, r2",
				"Add p-value dependent r^2-threshold, specify both p-value and r^2-threshold.",
				false, "r2t");
		
		//r2
		Option fixedR2 = createOptionOneName(1, "r2",
				"Add a fixed r^2-threshold (valid for all p-values)", false,
				"r2");
		
		//chromosome
		Option chromosome = createOptionOneName(
				1,
				"chromosome",
				"Specify chromosome, if none chosen - all chromosomes will be included",
				false, "chr");
		
		//out
		Option outputFile = createOptionOneName(1, "output_prefix",
				"Define output file prefix.", false, "out");
		
		//force_include
		Option forceIncludeFilePath = createOptionOneName(
				1,
				"filename",
				"Specify file path for a file containing SNPs that are to be force included",
				false, "force_include");
		
		//do_not_pick_force_included_first
		Option doNotPickForceIncludedFirst = createOptionOneName(
				0,
				"none",
				"Force-included SNPs are processed in the same order as all other SNPs (i.e., by increasing p-value)",
				false, "do_not_pick_force_included_first");
		
		//snp_table 
		Option snpTable = createOptionOneName(1, "input_filename",
				"Specify file path for SNP input table", false, "snp_table");
		
		//tped
		Option tped = createOptionOneName(1, "input_filename",
				"Specify file path for tped file", false, "tped");
		
		//tfam
		Option tfam = createOptionOneName(1, "input_filename",
				"Specify file path for tfam file", false, "tfam");
		
		//tfile
		Option tfile = createOptionOneName(1, "input_file_prefix",
				"Specify file prefix of tped and tfam files", false,
				"tfile");
		
		//add_surrogate_for_non_passing_index_snps
		Option additionalSurrogates = createOptionOneName(
				0,
				"none",
				"Specify that an additional surrogate will be added instead of choosing a new index SNP, in the case current index SNP does not pass the design score minimum",
				false, "add_surrogate_for_non_passing_index_snps");
		
		//verbose
		Option verbose = createOptionOneName(
				0,
				"none",
				"Specify that detailed information should be printed out and stored in the log file",
				false, "verbose");
		
		//outputLDTable
		Option outputLDTable = createOptionOneName(
				0,
				"none",
				"Directs the software to output a table containing all LD calculations made by the software",
				false, "ld");
		
		
		//remove
		Option remove = createOptionOneName(
				1,
				"filename",
				"Specify file path for a file containing family ID and individual ID of individuals that should be removed from current sample",
				false, "remove");
		
		//keep
		Option keep = createOptionOneName(
				1,
				"filename",
				"Specify file path for a file containing family ID and individual ID of individuals that should be kept, rest will be removed from current sample",
				false, "keep");
		
		//keep_random
		Option keepRandom = createOptionOneName(
				1,
				"float",
				"Specify a percentage value (as a decimal) to get a random selection of individuals of that size",
				false, "keep_random");
		
		//no_surrogates_for_force_included_snps
		Option noSurrogatesForForceIncluded = createOptionOneName(
				0,
				"none",
				"Specify that surrogates should not be added for force included SNPs",
				false, "no_surrogates_for_force_included_snps");
		
		//help
		Option help = createOptionTwoNames(0, "none", "help", "Print help",
				false, "h");

		// adds a special set of sample options that are mutually exclusive, to
		// an option group
		OptionGroup keepRemoveGroup = new OptionGroup();
		keepRemoveGroup.addOption(remove);
		keepRemoveGroup.addOption(keep);
		keepRemoveGroup.addOption(keepRandom);

		// adds the two different r^2-options to an option group, to make them
		// mutually exclusive
		OptionGroup r2Group = new OptionGroup();
		r2Group.addOption(r2Threshold);
		r2Group.addOption(fixedR2);
		//r2Group.setRequired(true);

		// adds all options to an option collection
		

		options.addOption(maxDistance);
		options.addOption(minMaf);
		//options.addOption(minHwPvalue);
		options.addOption(minSnpCallRate);
		options.addOption(minDesignScore);
		//options.addOption(absMinDesignScore);
		options.addOption(metric);
		options.addOption(surrogateThreshold);
		options.addOption(r2Threshold);
		options.addOption(chromosome);
		options.addOption(outputFile);
		options.addOption(forceIncludeFilePath);
		options.addOption(doNotPickForceIncludedFirst);
		options.addOption(snpTable);
		options.addOption(tped);
		options.addOption(tfam);
		options.addOption(tfile);
		options.addOption(additionalSurrogates);
		options.addOption(verbose);
		options.addOption(outputLDTable);
		options.addOption(noSurrogatesForForceIncluded);
		options.addOption(r2Threshold);
		options.addOption(fixedR2);
		options.addOption(help);
		options.addOptionGroup(keepRemoveGroup);
		options.addOptionGroup(r2Group);

		return options;
	}

	/**
	 * Parses command line arguments for the help-option and prints help message
	 * in case this option is encountered.
	 * 
	 * @param helpOptions
	 *            collection containing the help option
	 * @param options
	 *            collection of Option-objects (all other options)
	 * @param args
	 *            command line arguments
	 * @throws PriorityPrunerException
	 *             if problems are encountered during parsing
	 */
	@SuppressWarnings("unchecked")
	private boolean helpParse(Options helpOptions, Options options,
			String[] args) throws PriorityPrunerException {
		try {
			CommandLine helpCommandLine = new GnuParser().parse(helpOptions,
					args, true);
			//String empty = "";
			if (helpCommandLine.hasOption("h")) {
				printHelp();
				return true;
			} // checks list of arguments manually to make sure help is printed
				// correctly
			else {
				List<String> list = new ArrayList<String>();
				list = helpCommandLine.getArgList();
				if (list.contains("-h") || list.contains("--help")
						|| list.contains("-help")) {
					printHelp();
					return true;
				}
			}
		} catch (ParseException e) {
			throw new PriorityPrunerException(
					"Invalid command line options. Type \"-h\" option to view a list of available program options.");
		}
		return false;
	}

	private void printHelp(){
		
		System.out.println("PriorityPruner v0.1.0\r\n");
		System.out.println("(c) 2014 Christopher K. Edlund, Malin Anker, Fredrick R. Schumacher, W. James Gauderman, David V. Conti\r\n");
		System.out.println("For a list of options, please go to:");
		System.out.println("            http://sourceforge.net/projects/prioritypruner/\r\n");
		System.out.println("Questions, comments, and suggestions should be emailed to:");
		System.out.println("            PriorityPruner-help@lists.sourceforge.net");
		
	}
	/**
	 * Parses command line arguments, evaluates them, and stores values if
	 * they're correct.
	 * 
	 * @param options
	 *            collection of Option-objects
	 * @param args
	 *            command line arguments
	 * @throws PriorityPrunerException
	 *             if error occurs during parsing
	 */
	private void parse(Options options, String[] args)
			throws PriorityPrunerException {
		try {
			CommandLineParser cmdLineParser = new GnuParser();
			CommandLine commandLine = cmdLineParser.parse(options, args);
			setParsedOptions(commandLine.getOptions());
			// gets set to true if -tfile option is entered through the command
			// line
			//boolean tFlag = false;
			// counts --tped & --tfam options entered through the command line
			//int tpedCount = 0;
			//int tfamCount = 0;

			// loops through all options to save the ones in effect as a String
			// for printing
			String tmp = "Options in effect:";
			for (Option opt : commandLine.getOptions()) {
				tmp += ("\r\n	-" + opt.getOpt());
				if (opt.getValues() != null) {
					for (int i = 0; i < opt.getValues().length; i++) {
						tmp += (" " + opt.getValues()[i]);
					}
				}
			}


			// saves options in effect for printing
			this.setOptionsInEffect(tmp + "\r\n");

			// parse max_distance
			if (commandLine.hasOption("max_distance")) {
				this.maxDistance(getLongArgument("max_distance",
						commandLine.getOptionValue("max_distance"), 0,
						Long.MAX_VALUE));
				checkInput(1, "max_distance", commandLine);
			}

			// parse min_maf
			if (commandLine.hasOption("min_maf")) {
				this.setMinMaf(getDoubleArgument("min_maf",
						commandLine.getOptionValue("min_maf"), 0, 0.5));
				checkInput(1, "min_maf", commandLine);
			}

			// parse min_hwe
//			if (commandLine.hasOption("min_hwe")) {
//				this.setMinHwe(getDoubleArgument("min_hwe",
//						commandLine.getOptionValue("min_hwe"), 0, 1));
//				checkInput(1, "min_hwe", commandLine);
//			}

			// parses min_snp_callrate
			if (commandLine.hasOption("min_snp_callrate")) {
				this.setMinSnpCallRate(getDoubleArgument("min_snp_callrate",
						commandLine.getOptionValue("min_snp_callrate"), 0, 1));
				checkInput(1, "min_snp_callrate", commandLine);
			}

			// parses min_design_score
			if (commandLine.hasOption("min_design_score")) {
				this.setMinDesignScore(getDoubleArgument("min_design_score",
						commandLine.getOptionValue("min_design_score"), 0, Double.MAX_VALUE));
				checkInput(1, "min_design_score", commandLine);
			}

//			// parses option that sets the absolute minimum design score
//			// value, allowed arguments are doubles between 0.0 and
//			// Double.MAX_VALUE
//			if (commandLine.hasOption("amds")) {
//				this.setAbsoluteMinDesignScore(getDoubleArgument(
//						"absmindesignscore",
//						commandLine.getOptionValue("amds"), 0, Double.MAX_VALUE));
//				checkInput(1, "amds", commandLine);
//			}

			// parse metric
			if (commandLine.hasOption("metric")) {
				String[] str = commandLine.getOptionValues("metric");
				if (str.length % 2 == 1) {
					throw new PriorityPrunerException(
							"Only one argument specified for option \"metric\", a column name (text string) and a weight (decimal number) are required.");
				} else {
					for (int i = 0; i < str.length; i++) {
						if (str[i].equals("design_score")) {
							throw new PriorityPrunerException(
									"Invalid metric: \""
											+ str[i]
											+ "\". To filter on design score, simply add correct column in the SNP input file, no metric is needed.");
						}
						addMetric(str[i], this.getDoubleArgument("metric",
								str[i + 1], 0, Double.MAX_VALUE));
						i++;
					}
				}
			}

			// parse st
			if (commandLine.hasOption("st")) {
				String[] str = commandLine.getOptionValues("st");
				if (str.length % 2 == 1) {
					throw new PriorityPrunerException(
							"Only one argument specified for option: \"st\", a p-value (decimal number) and a number of surrogates (integer) are are required.");
				} else {
					for (int i = 0; i < str.length; i++) {
						addSurrogateThreshold(this.getDoubleArgument("st",
								str[i], 0, 1), this.getIntegerArgument("st",
								str[i + 1], 0, Integer.MAX_VALUE));
						i++;
					}
				}
			}

			// parse r2t
			if (commandLine.hasOption("r2t")) {
				String[] str = commandLine.getOptionValues("r2t");
				if (str.length % 2 == 1) {
					throw new PriorityPrunerException(
							"Only one argument specified for option: \"r2t\", a p-value and a threshold (decimal numbers between 0 and 1) are required.");
				} else {
					for (int i = 0; i < str.length; i++) {
						addR2Threshold(
								this.getDoubleArgument("r2t", str[i], 0, 1),
								this.getDoubleArgument("r2t", str[i + 1], 0, 1));
						i++;
					}
				}
			}

			// parse r2
			if (commandLine.hasOption("r2")) {
				String value = commandLine.getOptionValue("r2");
				this.addR2Threshold(1,
						this.getDoubleArgument("r2", value, 0, 1));
				checkInput(1, "r2", commandLine);
			}

			// parse chr
			if (commandLine.hasOption("chr")) {
				String value = commandLine.getOptionValue("chr");
				// recoding chromosome X-representations to "23"
				if (value.toUpperCase().equals("X")
						|| value.toUpperCase().equals("CHRX")) {
					//value = "23";
				}
				// if chromosome Y or mitochondrial DNA is encountered, an
				// exception is thrown
				else if (value.toUpperCase().equals("M")
						|| value.toUpperCase().equals("MT")
						|| value.toUpperCase().equals("CHRM")
						|| value.toUpperCase().equals("CHRMT")
						|| value.toUpperCase().equals("Y")
						|| value.toUpperCase().equals("CHRY")
						|| value.toUpperCase().equals("24")) {
					throw new PriorityPrunerException(
							"Chromosome \""
									+ value
									+ "\" specified in the command line, is not supported. Please update input files and rerun program.");
				}
				this.setChr(value);
				checkInput(1, "chr", commandLine);
			}

			// parse out
			if (commandLine.hasOption("out")) {
				String value = commandLine.getOptionValue("out");
				if (value.endsWith("/")) {
					value = new StringBuilder(value).append("prioritypruner")
							.toString();
				}
				this.setOutputPrefix(value);
				checkInput(1, "out", commandLine);
			}

//			// parse forceInclude
//			if (commandLine.hasOption("force_include")) {
//				String value = commandLine.getOptionValue("force_include");
//				this.setForceIncludeFilePath(value);
//				checkInput(1, "force_include", commandLine);
//			}

			// parse do_not_pick_force_included_snps_first
			if (commandLine.hasOption("do_not_pick_force_included_snps_first")) {
				this.setSortByForceIncludeAndPValue(false);
			}

			// parse snp_table
			if (commandLine.hasOption("snp_table")) {
				String[] str = commandLine.getOptionValues("snp_table");
				this.setSnpTablePath(str[0]);
				// counts number of metrics added to command line
				int metrics = 0;
				for (int i = 0; i < getParsedOptions().length; i++) {
					if (getParsedOptions()[i].getOpt().equals("metric")) {
						metrics++;
					}
				}
				this.setNumMetrics(metrics);
				checkInput(1, "snp_table", commandLine);
			}

			// parse tped
			if (commandLine.hasOption("tped")) {
				String value = commandLine.getOptionValue("tped");
				checkInput(1, "tped", commandLine);
				this.setTped(value);
			}

			// parse tfam
			if (commandLine.hasOption("tfam")) {
				String value = commandLine.getOptionValue("tfam");
				checkInput(1, "tfam", commandLine);
				this.setTfam(value);
			}

			// parse tfile
			if (commandLine.hasOption("tfile")) {
				String value = commandLine.getOptionValue("tfile");
				checkInput(1, "tfile", commandLine);
				this.setTped(value + ".tped");
				this.setTfam(value + ".tfam");
				this.setTfile(value);
			}

			// parse use_surrogate_for_non_passing_index_snps
//			if (commandLine.hasOption("use_surrogate_for_non_passing_index_snps")) {
//				this.setUseSurrogateForNonPassingIndexSnp(true);
//			}

			// parses verbose
			if (commandLine.hasOption("verbose")) {
				this.setVerbose(true);
			}
			
			// parse outputLDTable
			if (commandLine.hasOption("ld")){
				this.setOutputLDTable(true);
			}

			// parse remove
			if (commandLine.hasOption("remove")) {
				String value = commandLine.getOptionValue("remove");
				checkInput(1, "remove", commandLine);
				this.setRemove(value);
			}

			// parse keep
			if (commandLine.hasOption("keep")) {
				String value = commandLine.getOptionValue("keep");
				checkInput(1, "keep", commandLine);
				this.setKeep(value);
			}

			// parse keep_random
			if (commandLine.hasOption("keep_random")) {
				String value = commandLine.getOptionValue("keep_random");
				this.setKeepPercentage(this.getDoubleArgument("keep_random",
						value, 0, 1));
				checkInput(1, "keep_random", commandLine);
			}
			
			
			// parse no_surrogates_for_force_included_snps
			if (commandLine.hasOption("no_surrogates_for_force_included_snps")) {
				this.setAddSurrogatesForForceIncludedSnps(false);
			}

			// check that we have all required arguments
			checkRequiredArguments(commandLine);
			
			// checks if any unrecognized arguments been entered
			checkLeftArguments(commandLine);

			// if several options from the same options group been entered,
			// an AlreadySelectedException gets thrown. A custom message is
			// generated since we wanted another design than the one provided in
			// the library gave
		} catch (AlreadySelectedException e) {
			String message = "";
			for (int i = 0; i < e.getOptionGroup().getNames().toArray().length; i++) {
				message += "\"" + e.getOptionGroup().getNames().toArray()[i]
						+ "\" ";
			}
			throw new PriorityPrunerException(
					"The options: "
							+ message
							+ " may not both be defined. Type --help for help.");

			// if an undefined option is entered an UnrecognizedOptionException
			// gets thrown
		} catch (UnrecognizedOptionException e) {
			throw new PriorityPrunerException(e.getOption()
							+ " is not a valid option. Type --help for help.");

			// if an option that is supposed to have arguments is missing,
			// a MissingArgumentException gets thrown
		} catch (MissingArgumentException e) {
			throw new PriorityPrunerException("Missing argument for option \""
					+ e.getOption().getOpt() + "\". Expected: "
					+ e.getOption().getArgName()
					+ ". Type --help for help.");

			// if a required option is missing, a MissingOptionException gets
			// thrown
		} catch (MissingOptionException e) {
			// if any other problem occurs while parsing, a general
			// ParseException gets thrown
		} catch (ParseException parseException) {
			throw new PriorityPrunerException(
					"Invalid command line options. Type --help for help.");
		}
	}

	
	private void checkRequiredArguments(CommandLine commandLine) throws PriorityPrunerException{
		
		String messages = "";
		
		// check for genotype input
		if (this.tfam == null || this.tped == null){
			messages += "A genotype dataset must be specified with the --tfile or --tped/--tfam options.\r\n";
		}
		
		// check that both tfile and 
		if (commandLine.hasOption("tfile") && (commandLine.hasOption("tped") || commandLine.hasOption("tfam"))){
			messages += "Only the --tfile option or the --tped/--tfam options may be specified.\r\n";
		}
		
		// check for r2 
		if (this.sortedR2Thresholds.isEmpty()){
			messages += "At least one r-squared threshold must be specified with the --r2 or --r2t options.\r\n";
		}
		
		// check for out
		if (this.outputPrefix == null){
			this.outputPrefix = "prioritypruner";
		}
	
		// check snpTable
		if (this.snpTablePath == null){
			messages += "A SNP input table file must be specified with the --snp_table option.\r\n";
		}
		
		if (messages.length() >0){
			throw new PriorityPrunerException(messages);
		}
	}
	
	
	/**
	 * Parses a String-argument to a long value and makes sure it's valid.
	 * 
	 * @param option
	 *            name of the Option to which the argument belongs
	 * @param argument
	 *            String that will be parsed
	 * @param min
	 *            minimum value allowed for this parameter
	 * @param max
	 *            maximum value allowed for this parameter
	 * @return a long value
	 * @throws NumberFormatException
	 *             if invalid String-argument were provided
	 * @throws PriorityPrunerException
	 *             if invalid String-argument were provided
	 */
	private long getLongArgument(String option, String argument, long min,
			long max) throws PriorityPrunerException {
		try {
			Long longValue = Long.parseLong(argument);
			if (longValue >= min && longValue <= max) {
				return longValue;
			} else {
				// exception gets thrown either if parsing is impossible or
				// values are out of range
				throw new NumberFormatException();
			}
		} catch (NumberFormatException e) {
			throw new PriorityPrunerException(
							 argument
							+ " is not a valid input for option " + option + ", please specify an integer value between "
							+ min + " and " + max
							+ ". For help type --help.");
		}
	}

	/**
	 * Parses a String-argument to an integer value and makes sure it's valid.
	 * 
	 * @param option
	 *            name of the Option to that the argument belongs
	 * @param argument
	 *            String that will be parsed
	 * @param min
	 *            minimum value allowed for this parameter
	 * @param max
	 *            maximum value allowed for this parameter
	 * @return an integer value
	 * @throws NumberFormatException
	 *             if invalid String-argument were provided
	 * @throws PriorityPrunerException
	 *             if invalid String-argument were provided
	 */
	private int getIntegerArgument(String option, String argument, int min,
			int max) throws PriorityPrunerException {
		try {
			Integer intValue = Integer.parseInt(argument);
			if (intValue >= min && intValue <= max) {
				return intValue;
			} else {
				// exception get thrown either if parsing is impossible or
				// values are out of range
				throw new NumberFormatException();
			}
		} catch (NumberFormatException e) {
			throw new PriorityPrunerException(
					"\""
							+ argument
							+ "\" is not a valid input for option \""
							+ option
							+ "\", please specify an integer value between "
							+ min + " and " + max
							+ ".\r\nFor more information type \"-h\".");
		}
	}

	/**
	 * Parses a String-argument to a double value and makes sure it's valid.
	 * 
	 * @param option
	 *            name of the Option to that the argument belongs
	 * @param argument
	 *            String that will be parsed
	 * @param min
	 *            minimum value allowed for this parameter
	 * @param max
	 *            maximum value allowed for this parameter
	 * @return a double value
	 * @throws NumberFormatException
	 *             if invalid String-argument were provided
	 * @throws PriorityPrunerException
	 *             if invalid String-argument were provided
	 */
	private double getDoubleArgument(String option, String argument,
			double min, double max) throws PriorityPrunerException {
		try {
			Double doubleValue = Double.parseDouble(argument);
			if (doubleValue >= min && doubleValue <= max) {
				return doubleValue;
			} else {
				// exception get thrown either if parsing is impossible or
				// values are out of range
				throw new NumberFormatException();
			}
		} catch (NumberFormatException e) {
			throw new PriorityPrunerException(
					"\""
							+ argument
							+ "\" is not a valid input for option \r\n"
							+ option
							+ "\r\n, please specify an decimal value between "
							+ min + " and " + max
							+ ".\r\nFor more information type \"-h\".");
		}
	}

	/**
	 * Checks that options with one argument are entered correct number of times
	 * (1 time for options currently using this method) by counting
	 * corresponding arguments. If more than the specified number of arguments
	 * are provided in the command line an exception will get thrown.
	 * 
	 * @param numArgs
	 *            allowed number of times the option should be entered
	 * @param option
	 *            name of the option
	 * @param commandLine
	 *            CommandLine-object, which holds the command line arguments
	 * @throws PriorityPrunerException
	 *             if too many arguments were provided to the command line
	 */
	private void checkInput(int numArgs, String option, CommandLine commandLine)
			throws PriorityPrunerException {
		String[] args = commandLine.getOptionValues(option);
		if (args.length > numArgs) {
			throw new PriorityPrunerException(
					"Option '"
							+ option
							+ "' specified more than once. Type --help for help.");
		}
	}

	/**
	 * Checks if unexpected arguments were encountered during parsing, and if
	 * so, prints them out when throwing PriorityPrunerException.
	 * 
	 * @param commandLine
	 *            CommandLine-object, which holds the command line arguments
	 * @throws PriorityPrunerException
	 *             if unexpected arguments were found
	 */
	private void checkLeftArguments(CommandLine commandLine)
			throws PriorityPrunerException {

		String[] left = commandLine.getArgs();
		if (left.length > 0) {
			String leftArgs = "";
			for (int i = 0; i < left.length; i++) {
				leftArgs = leftArgs + left[i] + " ";
			}
			throw new PriorityPrunerException(
					"Unknown argument(s): "
							+ leftArgs
							+ ". Type --help for help.");
		}
	}

	/**
	 * Adds surrogate threshold to a sorted ArrayList called
	 * sortedSurrogateThresholds.
	 * 
	 * @param pValue
	 *            p-value for current threshold, represented as a double between
	 *            0.0 and 1.0
	 * @param numSurrogates
	 *            number of surrogates, 0 - no upper limit
	 */
	private void addSurrogateThreshold(double pValue, int numSurrogates) {
		this.sortedSurrogateThresholds.add(new SurrogateThreshold(pValue,
				numSurrogates));
		Collections.sort(sortedSurrogateThresholds);
	}

	/**
	 * Adds r^2 threshold to a sorted ArrayList called sortedR2Thresholds.
	 * 
	 * @param pValue
	 *            p-value for current threshold, represented as a double between
	 *            0.0 and 1.0
	 * @param r2Threshold
	 *            r^2 threshold represented as a double between 0.0 and 1.0
	 */
	private void addR2Threshold(double pValue, double r2Threshold) {
		sortedR2Thresholds.add(new R2Threshold(pValue, r2Threshold));
		Collections.sort(sortedR2Thresholds);
	}

	/**
	 * Creates and adds metric to an ArrayList called metricWeights. If a
	 * duplicated metric is encountered an exception gets thrown.
	 * 
	 * @param name
	 *            name of metric
	 * @param weight
	 *            weight of metric
	 * @throws PriorityPrunerException
	 *             if same metric is added more than once
	 */
	private void addMetric(String name, double weight)
			throws PriorityPrunerException {
		for (Metric m : metricWeights) {
			if (name.equals(m.getName())) {
				throw new PriorityPrunerException(
						"Metric \""
								+ name
								+ "\" added several times. Please update arguments in the command line and rerun program.");
			}
		}
		metricWeights.add(new Metric(name, weight));
	}

	// public getters and setters for private fields of this class

	public ArrayList<SurrogateThreshold> getSortedSurrogateThresholds() {
		return sortedSurrogateThresholds;
	}

	public ArrayList<R2Threshold> getSortedR2Thresholds() {
		return sortedR2Thresholds;
	}

	public long getMaxDistance() {
		return maxDistance;
	}

	public void maxDistance(long distance) {
		this.maxDistance = distance;
	}

	public double getMinMaf() {
		return minMaf;
	}

	public void setMinMaf(double minMaf) {
		this.minMaf = minMaf;
	}

//	public double getMinHwe() {
//		return minHwe;
//	}
//
//	public void setMinHwe(double minHwe) {
//		this.minHwe = minHwe;
//	}

	public double getMinSnpCallRate() {
		return minSnpCallRate;
	}

	public void setMinSnpCallRate(double minimumGenotypePercentage) {
		this.minSnpCallRate = minimumGenotypePercentage;
	}

	public double getMinDesignScore() {
		return minDesignScore;
	}

	public void setMinDesignScore(double minDesignScore) {
		this.minDesignScore = minDesignScore;
	}

//	public double getAbsoluteMinDesignScore() {
//		return absoluteMinDesignScore;
//	}
//
//	public void setAbsoluteMinDesignScore(double absoluteMinDesignScore) {
//		this.absoluteMinDesignScore = absoluteMinDesignScore;
//	}

	public ArrayList<Metric> getMetrics() {
		return metricWeights;
	}

	public String getChr() {
		return chr;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public String getOutputPrefix() {
		return outputPrefix;
	}

	public void setOutputPrefix(String outputFile) {
		this.outputPrefix = outputFile;
	}

//	public String getForceIncludeFilePath() {
//		return forceIncludeFilePath;
//	}
//
//	public void setForceIncludeFilePath(String forceIncludeFilePath) {
//		this.forceIncludeFilePath = forceIncludeFilePath;
//	}

	public Boolean getSortByForceIncludeAndPValue() {
		return sortByForceIncludeAndPValue;
	}

	public void setSortByForceIncludeAndPValue(Boolean value) {
		this.sortByForceIncludeAndPValue = value;
	}

	public String getSnpTablePath() {
		return snpTablePath;
	}

	public void setSnpTablePath(String snpFilePath) {
		this.snpTablePath = snpFilePath;
	}

	public int getNumMetrics() {
		return numMetrics;
	}

	public void setNumMetrics(int numMetrics) {
		this.numMetrics = numMetrics;
	}

	public String getTped() {
		return tped;
	}

	public void setTped(String tped) {
		this.tped = tped;
	}

	public String getTfam() {
		return tfam;
	}

	public void setTfam(String tfam) {
		this.tfam = tfam;
	}

	public String getTfile() {
		return tfile;
	}

	public void setTfile(String tfile) {
		this.tfile = tfile;
	}

	public Option[] getParsedOptions() {
		return parsedOptions;
	}

	public void setParsedOptions(Option[] parsedOptions) {
		this.parsedOptions = parsedOptions;
	}

//	public boolean isUseSurrogateForNonPassingIndexSnp() {
//		return useSurrogateForNonPassingIndexSnp;
//	}
//
//	public void setUseSurrogateForNonPassingIndexSnp(boolean value) {
//		this.useSurrogateForNonPassingIndexSnp = value;
//	}

	public boolean getVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	public boolean isOutputLDTable() {
		return this.outputLDTable;
	}

	public void setOutputLDTable(boolean outputLDTable) {
		this.outputLDTable = outputLDTable;
	}

	
	public String getRemove() {
		return remove;
	}

	public void setRemove(String remove) {
		this.remove = remove;
	}

	public String getKeep() {
		return keep;
	}

	public void setKeep(String keep) {
		this.keep = keep;
	}

	public double getKeepPercentage() {
		return keepPercentage;
	}

	public void setKeepPercentage(double keepPercentage) {
		this.keepPercentage = keepPercentage;
	}

	public String getOptionsInEffect() {
		return optionsInEffect;
	}

	public void setOptionsInEffect(String optionsInEffect) {
		this.optionsInEffect = optionsInEffect;
	}

	public boolean getAddSurrogatesForForceIncludedSnps() {
		return addSurrogatesForForceIncludedSnps;
	}

	public void setAddSurrogatesForForceIncludedSnps(boolean value) {
		this.addSurrogatesForForceIncludedSnps = value;
	}
}
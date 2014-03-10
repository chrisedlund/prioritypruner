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

import java.io.*;
import java.util.*;

/**
 * This class handles both the parsing of the SNP input file, as well as the
 * comparison between it and the TPlink data.
 */
public class SnpListFile {

	private String filePath;
	// number of metrics defined in command line
	private int numMetrics;
	private ArrayList<SnpInfo> snps = new ArrayList<SnpInfo>();
	private ArrayList<SnpInfo> snpsSortedByChrPos = new ArrayList<SnpInfo>();
	private ArrayList<String> chromsomes = new ArrayList<String>();
	// hashtable with SnpInfo-objects, which allows easy look up
	private Hashtable<String, ArrayList<SnpInfo>> snpNameToSnpInfoListHash = new Hashtable<String, ArrayList<SnpInfo>>();
	private BufferedReader reader;
	// allows access to options entered in the command line
	private CommandLineOptions options = CommandLineOptions.getInstance();
	// stores information from SNP input file about which SNPs that should be
	// force included
	private HashSet<String> forceIncludeHashSet = new HashSet<String>();
	// splitting regex for input files that allows single tabs or spaces
	private String delim = "[\\s|\\t]";
	private int designScorePlacement = -1;

	/**
	 * Constructor for SnpListFile. Initiates parsing of the SNP input file, if
	 * applicable - it initiates parsing of force include info, sets the force
	 * include flag in each SNP, and then initiates sorting of the SNPs.
	 * 
	 * @param filePath
	 *            file path for SNP input file
	 * @param numMetrics
	 *            number of metrics defined in command line
	 * @throws PriorityPrunerException
	 *             if problem are encountered during parsing
	 */
	public SnpListFile(String filePath, int numMetrics)
			throws PriorityPrunerException {
		this.filePath = filePath;
		this.numMetrics = numMetrics;
		parseFile();
		if (options.getForceIncludeFilePath() != null) {
			parseForceInclude();
			setForceInclude();
		}
		sortSnps();
	}

	/**
	 * Sorts the parsed SNPs according to the parameters defined by the user in
	 * command line.
	 */
	private void sortSnps() {
		// sort by force include then p-value, if a specific chromosome is
		// selected
		if (options.getSortByForceIncludeAndPValue()
				&& options.getChr() != null) {
			Collections.sort(snps, new ForceIncludePValueSorter());
			// sort by p-value then name, if a specific chromosome is selected
		} else if (!options.getSortByForceIncludeAndPValue()
				&& options.getChr() != null) {
			Collections.sort(snps);
		}
		// sort by force include, then by chromosome and p-value, if all
		// chromosomes are selected
		if (options.getSortByForceIncludeAndPValue()
				&& options.getChr() == null) {
			Collections.sort(snps, new FIChrPValueSorter());
			// sort by chromosome and then p-value, if all chromosomes are
			// selected
		} else if (!options.getSortByForceIncludeAndPValue()
				&& options.getChr() == null) {
			Collections.sort(snps, new ChrPValueSorter());
		}

		for (SnpInfo snpInfo : snps) {
			snpsSortedByChrPos.add(snpInfo);
		}
		// sort by chromosome then position
		Collections.sort(snpsSortedByChrPos, new PosSorter());
		for (int i = 0; i < snpsSortedByChrPos.size(); i++) {
			snpsSortedByChrPos.get(i).setSortedByPosIndex(i);
		}
	}

	/**
	 * Parses the SNP input file, checks that values are valid and if they are,
	 * creates a SnpInfo-object with these as parameters. The combination of SNP
	 * name, chromosome, position and alleles in a SnpInfo-object should be
	 * unique. To ensure this, a check for duplicates of this type is provided
	 * before the new SnpInfo-object is stored. In case a duplicate is
	 * encountered, a PriorityPrunerException will get thrown.
	 * 
	 * @throws PriorityPrunerException
	 *             if problems are encountered during parsing
	 */
	private void parseFile() throws PriorityPrunerException {
		try {
			reader = new BufferedReader(new FileReader(filePath));
			int lineNum = 1;
			String[] header = reader.readLine().split(delim);

			// checks that necessary columns are provided
			if (!header[0].toLowerCase().equals("snpname")
					|| !header[1].toLowerCase().equals("chr")
					|| !header[2].toLowerCase().equals("pos")
					|| !header[3].toLowerCase().equals("a1")
					|| !header[4].toLowerCase().equals("a2")
					|| !header[5].toLowerCase().equals("p")
					|| !header[6].toLowerCase().equals("num_assays")
					|| !header[7].toLowerCase().equals("force_include")) {
				throw new PriorityPrunerException(
						"Invalid column headers in SNP input file. First 8 columns should equal: \"snpname\", \"chr\", \"pos\", \"a1\", \"a2\", \"p\", \"num_assays\", \"force_include\".\nMake sure columns are separated by only one tab or space character.");
			}

			ArrayList<MetricNamePos> metricNames = new ArrayList<MetricNamePos>();
			// loops through all Metric-objects (the metrics specified in
			// the command line) and saves their name and position in the SNP
			// input file, to metricNames.
			// If gone through whole header in SNP input file without finding
			// current metric, an exception will get thrown
			outerLoop: for (int j = 0; j < options.getMetrics().size(); j++) {
				for (int i = 8; i < header.length; i++) {
					if (options.getMetrics().get(j).getName()
							.equals(header[i])) {
						metricNames.add(new MetricNamePos(header[i], i));
						continue outerLoop;
					}
				}
				throw new PriorityPrunerException(
						"Specified metric \""
								+ options.getMetrics().get(j).getName()
								+ "\" is missing in file \""
								+ filePath
								+ "\". \nPlease add it to the file or respecify your options in command line.");
			}

			// loops through header to check for duplicated metrics (column
			// names)
			for (int i = 0; i < header.length; i++) {
				for (int j = 0; j < header.length; j++) {
					if (header[i].equals(header[j]) && i != j) {
						throw new PriorityPrunerException(
								"Duplicated column name \""
										+ header[i]
										+ "\" in file \""
										+ filePath
										+ "\". \nPlease remove/rename incorrect column(s) and rerun program.");
					}
				}
			}

			// loops through header to be able to warn if metric provided in
			// SNP input file isn't in command line
			outerLoop: for (int i = 8; i < header.length; i++) {
				if (header[i].equals("design_score")) {
					designScorePlacement = i;
					continue outerLoop;
				}
				for (int j = 0; j < options.getMetrics().size(); j++) {
					if (header[i]
							.equals(options.getMetrics().get(j).getName())) {
						continue outerLoop;
					}
				}
				LogWriter
						.getLogger()
						.warn("\nMetric \""
								+ header[i]
								+ "\" is specified in \""
								+ filePath
								+ "\" but not in command line. \nTo use this metric, please specify it in command line. For more information about available program options type \"-h\".");
			}
			
			if(designScorePlacement==-1&&(options.getAbsoluteMinDesignScore()!=0 || options.getMinDesignScore()!=0)){
				LogWriter.getLogger().warn("\nThreshold value(s) for minimum design score and/or absolute minimum design score entered in command line \nwithout any design score data provided in the SNP input file. No design score filter could be applied.");
			}
			
			// reads and parses every line in the SNP input file, and checks and
			// stores this information
			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);
				lineNum += 1;
				// makes sure that at least the expected number of columns are
				// present
				if (splitString.length < (8 + numMetrics)) {
					throw new PriorityPrunerException(
							"Missing essential columns in \"" + filePath
									+ "\". Expected " + (numMetrics + 8)
									+ " columns, but found "
									+ splitString.length + " on line "
									+ lineNum
									+ ".\nPlease update SNP input file.");
				} else {
					// makes sure that values aren't missing for chosen metrics
					for (MetricNamePos metric : metricNames) {
						try {
							if (splitString[metric.getPos()].equals("")) {
								throw new PriorityPrunerException(
										"Missing values in the \""
												+ metric.getName()
												+ "\" column on line "
												+ lineNum
												+ " in \""
												+ filePath
												+ "\". \nPlease update your file or respecify the chosen metrics.");
							}
						} catch (ArrayIndexOutOfBoundsException e) {
							throw new PriorityPrunerException(
									"Missing values in the \""
											+ metric.getName()
											+ "\" column on line "
											+ lineNum
											+ " in \""
											+ filePath
											+ "\". \nPlease update your file or respecify the chosen metrics.");
						}
					}
				}
				// makes sure that information in mandatory columns isn't
				// missing
				for (int i = 0; i < 8; i++) {
					if (splitString[i].equals("")) {
						throw new PriorityPrunerException(
								"Missing mandatory information on line "
										+ lineNum
										+ " in \""
										+ filePath
										+ "\""
										+ ". \nPlease update your input file and rerun program.");
					}
				}

				String snpName = new String(splitString[0]);
				String chr = new String(splitString[1]);

				// in case chromosome X is encountered it gets recoded to "23"
				if (chr.toUpperCase().equals("X")
						|| chr.toUpperCase().equals("CHRX")
						|| chr.toUpperCase().equals("CHR_X")
						|| chr.toUpperCase().equals("X_NONPAR")
						|| chr.toUpperCase().equals("23")) {
					chr = "23";
				}
				// stores parsed chromosomes
				if (!chromsomes.contains(chr.toUpperCase())) {
					chromsomes.add(chr.toUpperCase());
				}
				// continues parsing if either the chromosome parameter provided
				// in command line matches the one in SNP input file, or if no
				// chromosome option was provided in command line (i.e., all
				// chromosomes in SNP input file will be parsed)
				if (CommandLineOptions.getInstance().getChr() == null
						|| chr.toUpperCase().equals(
								CommandLineOptions.getInstance().getChr()
										.toUpperCase())) {
					// parses base pair position
					int pos;
					try {
						pos = Integer.parseInt(new String(splitString[2]));
					} catch (NumberFormatException e) {
						throw new PriorityPrunerException(
								"Invalid value specified for 'pos' at line "
										+ lineNum + " in file \"" + filePath
										+ "\". Integer value expected.");
					}
					// parses alleles
					String allele1 = new String(splitString[3].toUpperCase());
					String allele2 = new String(splitString[4].toUpperCase());

					// parses p-value
					double pValue;
					try {
						pValue = Double.parseDouble(new String(splitString[5]));
					} catch (NumberFormatException e) {
						throw new PriorityPrunerException(
								"Invalid value specified for 'p' at line "
										+ lineNum + " in file \"" + filePath
										+ "\". Decimal value expected.");
					}
					// parses design score
					int numAssays;
					try {
						numAssays = Integer
								.parseInt(new String(splitString[6]));
						if (numAssays != 1 && numAssays != 2) {
							throw new NumberFormatException();
						}
					} catch (NumberFormatException e) {
						throw new PriorityPrunerException(
								"Invalid value specified for 'num_assays' at line "
										+ lineNum + " in file \"" + filePath
										+ "\". '1' or '2' expected.");
					}
					// parses force include flag
					boolean forceInclude;
					if (splitString[7].equals("0")) {
						forceInclude = false;
					} else if (splitString[7].equals("1")) {
						forceInclude = true;
					} else {
						throw new PriorityPrunerException(
								"Invalid value specified for 'force_include' at line "
										+ lineNum + " in file \"" + filePath
										+ "\". Binary value expected.");
					}

					// parses metric weights
					double[] metrics = new double[numMetrics];
					try {
						for (int i = 0; i < metricNames.size(); i++) {
							metrics[i] = Double.parseDouble(new String(
									splitString[metricNames.get(i).getPos()]));
						}

					} catch (NumberFormatException e) {
						throw new PriorityPrunerException(
								"Invalid value specified for customized metric at line "
										+ lineNum + " in file \"" + filePath
										+ "\". Decimal value expected.");
					}

					// creates new SnpInfo-object with the parsed info from one
					// line in SNP input file
					SnpInfo snp = new SnpInfo(snpName, chr, pos, allele1,
							allele2, pValue, numAssays, forceInclude, metrics);

					// checks if a design score column is provided in the SNP
					// input file, if it is - add this value to the
					// SnpInfo-object
					if (designScorePlacement != -1) {
						double designScore;
						try {
							designScore = Double.parseDouble(new String(
									splitString[designScorePlacement]));
						} catch (NumberFormatException e) {
							throw new PriorityPrunerException(
									"Invalid value specified for 'design_score' at line "
											+ lineNum + " in file \""
											+ filePath
											+ "\". Decimal value expected.");
						}
						snp.setDesignScore(designScore);
					}

					// checks whether the hashtable contains the SNP name,
					// if it doesn't, it adds it
					if (!snpNameToSnpInfoListHash.containsKey(snpName)) {
						snpNameToSnpInfoListHash.put(snpName,
								new ArrayList<SnpInfo>());
					}
					ArrayList<SnpInfo> snpInfoList = snpNameToSnpInfoListHash
							.get(snpName);
					// loops through the ArrayList corresponding to current SNP
					// name in snpNameToSnpInfoListHash, to determine if this is
					// a duplicate
					for (SnpInfo currentSnp : snpInfoList) {
						if ((currentSnp.getSnpName().equals(snpName)
								&& currentSnp.getChr().equals(chr) && currentSnp
								.getPos() == pos)
								&& ((currentSnp.getAllele1().equals(allele1) && currentSnp
										.getAllele2().equals(allele2)) || (currentSnp
										.getAllele1().equals(allele2) && currentSnp
										.getAllele2().equals(allele1)))) {
							reader.close();
							throw new PriorityPrunerException(
									"Duplicated SNP \""
											+ snpName
											+ "\" at line "
											+ lineNum
											+ " in file \""
											+ filePath
											+ "\". \nThe combination of SNP name, chromosome, base pair position, and alleles should be unique. \nPlease remove or update the information about this SNP in input files and rerun program.");
						}
					}
					// if not a duplicate, add the SNP to following data
					// structures:
					snpInfoList.add(snp);
					snpNameToSnpInfoListHash.put(snpName, snpInfoList);
					snps.add(snp);
				}
			}
			reader.close();
			if (CommandLineOptions.getInstance().getChr() != null
					&& !chromsomes.contains(CommandLineOptions.getInstance()
							.getChr().toUpperCase())) {
				throw new PriorityPrunerException(
						"No SNP at chromosome \""
								+ CommandLineOptions.getInstance().getChr()
								+ "\" could be found in \""
								+ filePath
								+ "\". \nPlease respecify the argument in the command line or update input files.");
			}
		} catch (FileNotFoundException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage()
					+ "\nPlease check that correct file path is provided.");
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage()
					+ "\nPlease check that correct file path is provided.");
		}
	}

	/**
	 * Parses the force include file, in case that option is chosen in command
	 * line.
	 * 
	 * @throws PriorityPrunerException
	 *             if problems are encountered during parsing
	 */
	private void parseForceInclude() throws PriorityPrunerException {
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(
					options.getForceIncludeFilePath()));
			int index = 1;
			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);
				// expecting one SNP per line, comments are allowed and won't be
				// parsed
				if (splitString.length < 1) {
					reader.close();
					throw new PriorityPrunerException(
							"Missing SNP name at line: "
									+ (index)
									+ " in file \""
									+ filePath
									+ "\". \nPlease update file and rerun program.");
				}
				String snpName = new String(splitString[0]);

				forceIncludeHashSet.add(snpName);
				index++;
			}
			reader.close();
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		}
	}

	/**
	 * Sets the "force include"-flag in each SnpInfo-object according to the
	 * force include file.
	 */
	private void setForceInclude() {
		for (SnpInfo snp : snps) {
			if (forceIncludeHashSet.contains(snp.getSnpName())) {
				snp.setForceInclude(true);
			}
		}
	}

	/**
	 * This method makes sure information about a certain SNP matches in both
	 * the tped and SNP input file. It gets called by the TPlink-class during
	 * parsing, and if a SnpInfo-object with correct properties is found, it's
	 * returned to TPlink.
	 * 
	 * @param snpName
	 *            name of this SNP as provided in tped
	 * @param chr
	 *            chromosome of this SNP as provided in tped
	 * @param pos
	 *            base pair position provided in tped
	 * @param allele1
	 *            allele 1 provided in tped
	 * @param allele2
	 *            allele 2 provided in tped
	 * @return a matching SnpInfo-object
	 */
	public SnpInfo getSnpInfo(String snpName, String chr, int pos,
			String allele1, String allele2) {

		if (snpNameToSnpInfoListHash.containsKey(snpName)) {
			ArrayList<SnpInfo> snpInfoList = snpNameToSnpInfoListHash
					.get(snpName);

			// loops through all SnpInfo-objects and first checks that
			// chromosome and position matches
			for (SnpInfo snpInfo : snpInfoList) {
				if ((snpInfo.getChr().toUpperCase().equals(chr.toUpperCase()) && snpInfo
						.getPos() == pos)) {
					// if allele 1 is provided in tped, but allele 2 is missing,
					// check if allele 1 matches any of the two alleles in SNP
					// input file. If so, return the SnpInfo-object
					if (!allele1.equals("0") && allele2.equals("0")) {
						if ((snpInfo.getAllele1().equals(allele1) || snpInfo
								.getAllele2().equals(allele1))) {
							return snpInfo;
						}
						// if allele 2 is provided in tped, but allele 1 is
						// missing, check if allele 2 matches any of the two
						// alleles in SNP input file. If so, return the
						// SnpInfo-object
					} else if (!allele2.equals("0") && allele1.equals("0")) {
						if (snpInfo.getAllele1().equals(allele2)
								|| snpInfo.getAllele2().equals(allele2)) {
							return snpInfo;
						}
						// if either none or both of the alleles are missing
						// check that the information in tped matches SNP input
						// file. If so, return the SnpInfo-object
					} else if ((snpInfo.getAllele1().equals(allele1) && snpInfo
							.getAllele2().equals(allele2))
							|| (snpInfo.getAllele1().equals(allele2) && snpInfo
									.getAllele2().equals(allele1))) {
						return snpInfo;
					}
				}
			}
		}
		// if SNP name isn't in the hash or if it didn't match all criteria
		return null;
	}

	// public access to two of the private fields of this class

	public ArrayList<SnpInfo> getSnps() {
		return snps;
	}

	public ArrayList<SnpInfo> getSnpsSortedByChrPos() {
		return snpsSortedByChrPos;
	}
}
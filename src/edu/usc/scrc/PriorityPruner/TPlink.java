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
import java.util.HashSet;

/**
 * This class handles the parsing of the transposed plink-formats: tped & tfam.
 * The tped-format is used to store information regarding SNPs and its genotypes
 * and the tfam-format is used for family data. This class is a subclass of
 * Genotypes, which stores all the parsed information. This class supports
 * partial parsing of genotypes trough the options: remove, keep and keep_random
 * that help select subsets of individuals used during tped parsing. These
 * options are mutually exclusive by definition and are specified as arguments
 * to the command line. Both remove and keep requires that a file path to an
 * indlist is specified.
 */
public class TPlink extends Genotypes {

	private HashSet<String> keepRemoveHashSet = new HashSet<String>();
	private HashSet<String> individualHashSet = new HashSet<String>();
	private boolean remove;
	private boolean keep;
	private boolean keep_random;
	private CommandLineOptions options = CommandLineOptions.getInstance();

	// splitting regex for input files that allows single tabs or spaces
	private String delim = "[\\s|\\t]";

	/**
	 * Constructor used to initialize parsing according to specified options.
	 * 
	 * @param filePathTPed
	 *            relative or absolute file path for the tped file
	 * @param filePathTFam
	 *            relative or absolute file path for the tfam file
	 * @param snpListFile
	 *            the SNP input file, providing information about which SNPS in
	 *            the tped file to parse
	 * @throws PriorityPrunerException
	 *             if files aren't found or if problems are encountered during
	 *             parsing
	 */
	public TPlink(String filePathTPed, String filePathTFam,
			SnpListFile snpListFile) throws PriorityPrunerException {
		super();

		// checks if any of the remove, keep or keep_random options has been
		// chosen
		if (options.getKeep() != null || options.getRemove() != null
				|| options.getKeepPercentage() != -1) {
			parseKeepRemove();
		}

		// initiates tfam parsing
		parseTfam(filePathTFam);

		// sets the keep-flag in Individual-objects according
		// to keep-/remove-file
		setKeepRemove();

		// initiates tped parsing
		parseTped(filePathTPed, snpListFile);
	}

	/**
	 * Sets the boolean value to true for the option in effect (remove, keep or
	 * keep_random). If the option is remove or keep, the method also parses the
	 * associated file.
	 * 
	 * @throws PriorityPrunerException
	 *             if the associated file isn't found
	 */
	private void parseKeepRemove() throws PriorityPrunerException {

		String filePath = "";
		BufferedReader reader;

		// checks which option that has been entered
		if (options.getRemove() != null) {
			filePath = options.getRemove();
			remove = true;
		} else if (options.getKeep() != null) {
			filePath = options.getKeep();
			keep = true;
		} else {
			keep_random = true;
			return;
		}
		// parses the file with the individuals to remove/keep
		try {
			reader = new BufferedReader(new FileReader(filePath));
			int index = 1;
			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);

				if (splitString.length < 2) {
					reader.close();
					throw new PriorityPrunerException(
							"Invalid number of columns specified in file \""
									+ filePath
									+ "\" at line: "
									+ (index)
									+ "\nExpected: 2 (Family ID, Individual ID) Found: "
									+ splitString.length + ".");
				}
				String famID = new String(splitString[0]);
				String indID = new String(splitString[1]);

				if (famID.equals("") || indID.equals("")) {
					reader.close();
					throw new PriorityPrunerException(
							"Invalid formatting in file \""
									+ filePath
									+ "\" at line: "
									+ (index)
									+ "\nPlease check that values are separated by single space or tab characters only.");
				}
				// stores family ID + individual ID together as a unique key, to
				// enable easy and accurate retrieval of these individuals
				keepRemoveHashSet.add(famID + " " + indID);
				index++;
			}
			reader.close();
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		}
	}

	/**
	 * Sets the keep-flag in Individual-objects according to the keep/remove
	 * input file. It also stores the gender and index of each individual (in
	 * variables accessible from Genotypes).
	 */
	private void setKeepRemove() {
		int index = 0;

		// loops through all individuals to set their keep-flags, and save
		// gender and index
		for (Individual individual : individuals) {
			if (keep) {
				if (keepRemoveHashSet.contains(individual.getFamilyID() + " "
						+ individual.getIndividualID())) {
					individual.setKeep(true);
					subjectSexes.add(individual.getGender());
					founderIndices.add(index);
					index++;
				} else {
					individual.setKeep(false);
				}
			} else if (remove) {
				if (keepRemoveHashSet.contains(individual.getFamilyID() + " "
						+ individual.getIndividualID())) {
					individual.setKeep(false);
				} else {
					individual.setKeep(true);
					subjectSexes.add(individual.getGender());
					founderIndices.add(index);
					index++;
				}
			} else if (keep_random) {
				individual.setKeep(false);
				index++;
			} else {
				individual.setKeep(true);
				subjectSexes.add(individual.getGender());
				founderIndices.add(index);
				index++;
			}
		}

		// chooses random individuals for the option keep_random
		if (keep_random) {
			int numKeep = Math
					.round((float) (options.getKeepPercentage() * index));
			int founderIndex = 0;

			for (int i = 0; i < numKeep; i++) {
				int random = (int) (Math.random() * index);
				while (individuals.get(random).getKeep() == true) {
					random = (int) (Math.random() * index);
				}
				individuals.get(random).setKeep(true);
				subjectSexes.add(individuals.get(random).getGender());
				founderIndices.add(founderIndex);
				founderIndex++;
			}
		}
	}

	/**
	 * Parses the tped file. Since it's important that the tped file is
	 * correctly formatted, several checks for that are provided in this method.
	 * 
	 * @param filePath
	 *            relative or absolute file path for tped file
	 * @param snpListFile
	 *            relative or absolute file path for SNP input file
	 * @throws PriorityPrunerException
	 *             if any problems are encountered during parsing
	 */
	private void parseTped(String filePath, SnpListFile snpListFile)
			throws PriorityPrunerException {
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(filePath));
			int index = 0;

			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);
				// checks that no double tabs or spaces been entered in tped
				// file
				for (String str : splitString) {
					if (str.equals("")) {
						reader.close();
						throw new PriorityPrunerException(
								"Invalid formatting in file \""
										+ filePath
										+ "\" at line: "
										+ (index + 1)
										+ "\nPlease check that values are separated by single space or tab characters only.");
					}
				}

				// checks that number of columns in tped file is not less than
				// the required four
				if (splitString.length < 4) {
					reader.close();
					throw new PriorityPrunerException(
							"Missing required columns (\"Chromsome\", \"SNP Name\", \"Distance\", \"Position\") in file \""
									+ filePath
									+ "\" (at line: "
									+ (index + 1)
									+ ").");
				}

				// checks that the number of genotypes in tped file matches
				// number of individuals
				if (((splitString.length - 4) / 2 != individuals.size())) {
					reader.close();
					throw new PriorityPrunerException(
							"The number of individuals and genotype data are not matching: \nNumber of individuals provided in tfam file: "
									+ individuals.size()
									+ "\nNumber of genotypes provided in \""
									+ filePath
									+ "\" : "
									+ (splitString.length - 4)
									/ 2
									+ " (at line: " + (index + 1) + ").");
				}

				// stores chromosome X as "23"
				String chr = new String(splitString[0]);
				if (chr.toUpperCase().equals("X")
						|| chr.toUpperCase().equals("CHRX")
						|| chr.toUpperCase().equals("CHR_X")
						|| chr.toUpperCase().equals("X_NONPAR")
						|| chr.toUpperCase().equals("23")) {
					chr = "23";
				}

				// checks if all chromosomes should be parsed or if a specific
				// chromosome is specified in command line, and if it in that
				// case matches the chromosome on the line in the tped where
				// we are. If it does - continue parsing this line, else - go to
				// next line
				if (CommandLineOptions.getInstance().getChr() == null
						|| chr.toUpperCase().equals(
								CommandLineOptions.getInstance().getChr()
										.toUpperCase())) {
					String snpName = new String(splitString[1]);
					int pos;
					try {
						pos = Integer.parseInt(new String(splitString[3]));
					} catch (NumberFormatException e) {
						reader.close();
						throw new PriorityPrunerException(
								"Invalid value: \""
										+ splitString[3]
										+ "\", specified for base pair position in file \""
										+ filePath + "\" at line: "
										+ (index + 1)
										+ ". \nInteger value expected.");
					}
					String allele1 = "0";
					String allele2 = "0";
					String[] genotypes = new String[2 * founderIndices.size()];
					int individualIndex = 0;
					int genotypesIndex = 0;

					// goes through genotypes for this SNP (from individuals
					// that are set to be kept) and checks that no more than two
					// alleles are provided. Correct genotypes are also stored
					// in a String-array which will be saved together with other
					// info in a SnpGenotypes-object.
					for (int k = 4; k < splitString.length; k += 2) {
						if (individuals.get(individualIndex).getKeep()) {

							for (int j = 0; j < 2; j++) {
								if (allele1.equals("0")) {
									allele1 = new String(splitString[k + j]);
								} else if (allele2.equals("0")
										&& !allele1.equals(splitString[k + j])) {
									allele2 = new String(splitString[k + j]);
								}
								if (!allele1.equals(splitString[k + j])
										&& !allele2.equals(splitString[k + j])
										&& !splitString[k + j].equals("0")) {
									reader.close();
									throw new PriorityPrunerException(
											"To many types of alleles provided in file \""
													+ filePath
													+ "\" at line: "
													+ (index + 1)
													+ "\nFound: \""
													+ allele1
													+ "\", \""
													+ allele2
													+ "\", \""
													+ splitString[k + j]
													+ "\"\nA maximum of 2 alleles are allowed.");
								}
								genotypes[genotypesIndex] = new String(
										splitString[k + j]);
								genotypesIndex++;
							}
						}
						individualIndex++;
					}

					// gets SnpInfo-object from SnpListFile, if there is a
					// matching SNP in there
					SnpInfo snpInfo = snpListFile.getSnpInfo(snpName, chr, pos,
							allele1, allele2);

					// if the SnpInfo-object we got from SnpListFile is
					// not null, we incorporate it in a SnpGenotypes-object
					if (snpInfo != null) {
						SnpGenotypes snpGenotypesLocal = new SnpGenotypes(
								snpName, snpInfo, allele1, allele2, genotypes);
						snpInfo.setSnpGenotypes(snpGenotypesLocal);
						snpInfo.setInTped(true);
						snpGenotypes.add(snpGenotypesLocal);
					}
				}
				index++;
				splitString = null;
			}
			reader.close();
		} catch (FileNotFoundException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		}
	}

	/**
	 * Parses the tfam file. Since it's important that the tfam file is
	 * correctly formatted, several checks for that are provided in this method.
	 * 
	 * @param filePath
	 *            relative or absolute file path for tfam file
	 * @throws PriorityPrunerException
	 *             if problem is encountered during parsing
	 */
	private void parseTfam(String filePath) throws PriorityPrunerException {
		BufferedReader reader;

		try {
			reader = new BufferedReader(new FileReader(filePath));
			int index = 0;

			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);
				// checks that no double tabs or spaces been entered in tped
				// file
				for (String str : splitString) {
					if (str.equals("")) {
						reader.close();
						throw new PriorityPrunerException(
								"Invalid formatting in file \""
										+ filePath
										+ "\" at line: "
										+ (index + 1)
										+ "\nPlease check that values are separated by single space or tab characters only.");
					}
				}
				// checks that correct number of columns are provided
				if (splitString.length != 6) {
					reader.close();
					throw new PriorityPrunerException(
							"Invalid number of columns specified in file \""
									+ filePath
									+ "\" at line: "
									+ (index + 1)
									+ ".\nExpected: 6 (Family ID, Individual ID, Paternal ID, Maternal ID, Sex [1=male; 2=female], Phenotype). \nFound: "
									+ splitString.length + ".");
				}
				// checks that gender information is correct, it has to be
				// specified
				if (!splitString[4].equals("1") && !splitString[4].equals("2")) {
					reader.close();
					throw new PriorityPrunerException(
							"Invalid value specified for 'gender': "
									+ splitString[4]
									+ ", at line: "
									+ (index + 1)
									+ " in file \""
									+ filePath
									+ "\". \nIf male: '1' expected, if female: '2' expected.");
				}
				// checks that only founders are provided
				if (!splitString[2].equals("0") || !splitString[3].equals("0")) {
					reader.close();
					throw new PriorityPrunerException(
							"Only founders allowed, please remove/update information about individual at line: "
									+ (index + 1)
									+ " in file \""
									+ filePath
									+ "\".");
				}
				// checks that no duplicates get entered
				if (!individualHashSet
						.contains((splitString[0] + " " + splitString[1]))) {
					String famID = new String(splitString[0]);
					String indID = new String(splitString[1]);

					Individual individual = new Individual(famID, indID,
							new String(splitString[2]), new String(
									splitString[3]), new String(splitString[4]));
					individualHashSet.add(famID + " " + indID);
					individuals.add(individual);
					index++;
				} else {
					reader.close();
					throw new PriorityPrunerException(
							"Duplicate individual data found in file \""
									+ filePath
									+ "\" at line: "
									+ (index + 1)
									+ "\nPlease remove data, update tped file and rerun program.");
				}
			}
			reader.close();
		} catch (FileNotFoundException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		}
	}
}
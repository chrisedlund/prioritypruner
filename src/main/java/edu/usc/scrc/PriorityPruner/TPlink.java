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


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

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
	 * @param keepRemoveSamples
	 *            The PlinkSampleList containing the list samples to keep or remove
	 *            based on the --keep or --remove options. If null, then no list
	 *            is defined.
	 * @throws PriorityPrunerException
	 *             if files aren't found or if problems are encountered during
	 *             parsing
	 */
	public TPlink(String filePathTPed, String filePathTFam,
			SnpListFile snpListFile, PlinkSampleListFile keepRemoveSamples, CommandLineOptions options) throws PriorityPrunerException {
		super(keepRemoveSamples, options);

		LogWriter.getLogger().info("Reading pedigree information from [ " + filePathTFam + " ]");
		

		// parse the tfam file
		parseTfam(filePathTFam);
		
		// checks if any of the remove, keep or keep_random options has been
		// chosen
//		if (options.getKeep() != null || options.getRemove() != null
//				|| options.getKeepPercentage() > 0) {
//			parseKeepRemove();
//		}
		

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
//	private void parseKeepRemove() throws PriorityPrunerException {
//
//		String filePath = "";
//		BufferedReader reader;
//
//		// checks which option that has been entered
//		if (options.getRemove() != null) {
//			filePath = options.getRemove();
//			remove = true;
//		} else if (options.getKeep() != null) {
//			filePath = options.getKeep();
//			keep = true;
//		} else {
//			keep_random = true;
//			return;
//		}
//		// parses the file with the individuals to remove/keep
//		try {
//			reader = new BufferedReader(new FileReader(filePath));
//			int index = 1;
//			while (reader.ready()) {
//				String[] splitString = reader.readLine().split(delim);
//
//				if (splitString.length != 2) {
//					reader.close();
//					throw new PriorityPrunerException(
//							"Invalid number of columns specified in "
//									+ filePath
//									+ " at line "
//									+ (index)
//									+ ". Expected 2 columns but found: "
//									+ splitString.length + ".");
//				}
//				String famID = new String(splitString[0]);
//				String indID = new String(splitString[1]);
//
//				if (famID.equals("") || indID.equals("")) {
//					reader.close();
//					throw new PriorityPrunerException(
//							"Invalid formatting in file \""
//									+ filePath
//									+ "\" at line: "
//									+ (index)
//									+ "\nPlease check that values are separated by single space or tab characters only.");
//				}
//				// stores family ID + individual ID together as a unique key, to
//				// enable easy and accurate retrieval of these individuals
//				keepRemoveHashSet.add(famID + " " + indID);
//				index++;
//			}
//			reader.close();
//		} catch (IOException e) {
//			throw new PriorityPrunerException("Could not open file: "
//					+ e.getMessage());
//		}
//	}

	/**
	 * This is a transposed plink file, so we've already parsed the sample list (TFAM file)
	 *  before parsing the bulk of the data. For each individual, set its 
	 * keep flag based on the mutually exclusive --keep, --remove, or --keep_random flags. 
	 * Set the keptFounders member variable.
	 */
	private void setKeepRemove() {
		//int index = 0;
		int numKept = 0;
		int numRemoved = 0;
		
		
		if (options.getKeep() != null ){ // the user has defined a --keep file
			// loop through each individual and set its keep flag to true
			//  if it exists in the keepRemoveSamples member variable, otherwise set to false
			for (Individual individual: this.individuals){
				if (this.keepRemoveSamples.contains(individual.getFamilyID(), individual.getIndividualID())){
					individual.setKeep(true);
					numKept++;
				}else{
					individual.setKeep(false);
				}
			}
		}else if (options.getRemove() != null){ // the user has defined a --remove file
			// loop through each individual and set its keep flag to false
			//  if it exists in the keepRemoveSamples member variable, otherwise set to true
			for (Individual individual: this.individuals){
				if (this.keepRemoveSamples.contains(individual.getFamilyID(), individual.getIndividualID())){
					individual.setKeep(false);
					numRemoved++;
				}else{
					individual.setKeep(true);
				}
			}
		}else if (options.getKeepPercentage() > 0){ // the user has specified the --keep_random option
			// create a temporary array list with all individuals
			ArrayList<Individual> tempList = new ArrayList<Individual>(this.individuals.size());
			for (Individual individual: this.individuals){
				tempList.add(individual);
			}
			
			// randomly shuffle the temporary list, use a seed if defined by the user
			Random random;
			if (options.getSeed() != null){
				random = new Random(options.getSeed());
				LogWriter.getLogger().info("Using " + options.getSeed() + " as seed for randomly selecting individuals");
			}else{
				random = new Random();
			}
			Collections.shuffle(tempList, random);
				
			// determine the number of individuals we need to keep 
			int numKeep = Math.round((float) (options.getKeepPercentage() * this.individuals.size()));
			
			// for the first numKeep individuals, set keep to true; for the rest set keep to false
			numKept = numKeep;
			for (int i = 0; i < tempList.size(); i++){
				Individual individual = tempList.get(i);
				individual.setKeep(i < numKeep);
			}
		}
		
//		// loops through all individuals to set their keep-flags, and save
//		// gender and index
//		for (Individual individual : individuals) {
//			if (keep) {
//				if (keepRemoveHashSet.contains(individual.getFamilyID() + " "
//						+ individual.getIndividualID())) {
//					individual.setKeep(true);
//					numKept++;
//					index++;
//				} else {
//					individual.setKeep(false);
//				}
//			} else if (remove) {
//				if (keepRemoveHashSet.contains(individual.getFamilyID() + " "
//						+ individual.getIndividualID())) {
//					individual.setKeep(false);
//					numRemoved++;
//				} else {
//					individual.setKeep(true);
//					index++;
//				}
//			} else if (keep_random) {
//				individual.setKeep(false);
//				index++;
//			} else {
//				individual.setKeep(true);
//				index++;
//			}
//		}
//
//		// chooses random individuals for the option keep_random
//		if (keep_random) {
//			
//			final class RandomInd implements Comparable<RandomInd>{
//				public Individual individual;
//				public Double random;
//				public RandomInd(Individual individual, double random){this.individual = individual; this.random = random;}
//				@Override
//				public int compareTo(RandomInd arg0) {return this.random.compareTo(arg0.random);}
//			}
//			
//			ArrayList<RandomInd> randomIndList = new ArrayList<RandomInd>();
//			for (Individual individual: individuals){
//				randomIndList.add(new RandomInd(individual, new Double(Math.random())));
//			}
//			Collections.sort(randomIndList);
//			
//			int numKeep = Math.round((float) (options.getKeepPercentage() * index));
//
//			for (int i = 0; i < numKeep; i++) {
//				randomIndList.get(i).individual.setKeep(true);
//				numKept++;
//			}
//		}
		
		for (Individual ind: individuals){
			if (ind.getKeep()){
				keptFounders.add(ind);
			}
		}
		
		if (options.getKeep() != null){
			LogWriter.getLogger().info("Reading individuals to keep [ " + options.getKeep() + " ] ... " + numKept + " read");
		}else if (options.getRemove() != null){
			LogWriter.getLogger().info("Reading individuals to remove [ " + options.getRemove() + " ] ... " + numRemoved + " read");
		}else if (options.getKeepPercentage() > 0){
			LogWriter.getLogger().info("Selecting " + numKept + " random individuals to keep");
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
		BufferedReader reader = null;
		//HashMap<String,Integer> uniqueSnpNameHash = new HashMap<String,Integer>();
		
		int notFoundInSnpInputTable = 0;
		try {
			reader = new BufferedReader(new FileReader(filePath));
			int line = 1;

			LogWriter.getLogger().info("Reading genotypes from [ " + filePath + " ]");
			
			if (this.options.getChr() != null){
				LogWriter.getLogger().info("Extracting SNPs from chromosome " + this.options.getChr());
			}
			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);
				// checks that no double tabs or spaces been entered in tped file
				for (String str : splitString) {
					if (str.length() == 0) {
						throw new PriorityPrunerException(
								"Problem with line " + line + " in [ " + filePath + " ]\r\n"
										+ "Ensure values are separated by a single space or tab character.");
					}
				}
				
				// checks that the number of columns is what we expect
				int expectedColumns = individuals.size() * 2 + 4;
				if (splitString.length != expectedColumns) {
					throw new PriorityPrunerException(
							"Problem with line " + line + " in [ " + filePath + " ]\r\n"
									+ "Expecting 4 + 2 * " + individuals.size() + " = " + expectedColumns
									+ " columns, but found " + splitString.length);
				}

				// stores chromosome X as "23"
				String chr = new String(splitString[0]);
//				if (chr.toUpperCase().equals("X")
//						|| chr.toUpperCase().equals("CHRX")
//						|| chr.toUpperCase().equals("23")) {
//					chr = "23";
//				}

				// checks if all chromosomes should be parsed or if a specific
				// chromosome is specified in command line, and if it in that
				// case matches the chromosome on the line in the tped where
				// we are. If it does - continue parsing this line, else - go to
				// next line
				if (this.options.getChr() == null
						|| chr.toUpperCase().equals(
								this.options.getChr()
										.toUpperCase())) {
					
					String snpName = new String(splitString[1]);
					
//					//make sure there are no duplicate snps
//					if (uniqueSnpNameHash.containsKey(snpName)){
//						throw new PriorityPrunerException("Duplicate SNP found in tped file: " + snpName);
//					}else{
//						uniqueSnpNameHash.put(snpName, 0);
//					}
					
					int pos;
					try {
						pos = Integer.parseInt(new String(splitString[3]));
						if (pos < 1){
							throw new NumberFormatException();
						}
					} catch (NumberFormatException e) {
						throw new PriorityPrunerException(
								"Problem with line " + line + " in [ " + filePath + " ]\r\n"
								+ "Invalid value: \""
										+ splitString[3]
										+ "\", specified for position in column 4.");
					}
					String allele1 = "0";
					String allele2 = "0";
					String[] genotypes = new String[2 * keptFounders.size()];
					int individualIndex = 0;
					int genotypesIndex = 0;

					// goes through genotypes for this SNP (from individuals
					// that are set to be kept) and checks that no more than two
					// alleles are provided. Correct genotypes are also stored
					// in a String-array which will be saved together with other
					// info in a SnpGenotypes-object.
					for (int k = 4; k < splitString.length; k += 2) {
						Individual individual = individuals.get(individualIndex);
						if (individual.getKeep()) {

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
											"Locus " + snpName + " has >2 alleles:\r\n"
													+ "individual " + individual.getFamilyID() + " " + individual.getIndividualID() 
													+ " has genotype [ " + splitString[k] + " " + splitString[k + 1] + " ]\r\n"
													+ "but we've already seen [ " + allele1 + " ] and [ " + allele2 + " ]");
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
						if (snpInfo.getSnpGenotypes() != null){
							throw new PriorityPrunerException(
									"Duplicated SNP \""
											+ snpName
											+ "\" at line "
											+ line
											+ " in TPED file. " 
											+ "The combination of snpname, chr, pos, allele1/allele2 must be unique.");
						}
						snpInfo.setSnpGenotypes(snpGenotypesLocal);
						snpInfo.setInTped(true);
						snpGenotypes.add(snpGenotypesLocal);
					}else{
						notFoundInSnpInputTable++;
					}
				}
				line++;
				splitString = null;
			}
			
			
			LogWriter.getLogger().info("Excluding " + notFoundInSnpInputTable + " SNPs missing from [ " + snpListFile.getFilePath() + " ]");
			
			
			LogWriter.getLogger().info(snpGenotypes.size() + " (of " + (line - 1) + ") SNPs to be included from [ " + filePath + " ]");
			
		} catch (FileNotFoundException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} catch (PriorityPrunerException e){
			throw e;
		} finally{
			try {
				if (reader != null){
					reader.close();
				}
			} catch (IOException e) {
			}
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
		BufferedReader reader = null;

		try {
			reader = new BufferedReader(new FileReader(filePath));
			int line = 1;
			int maleCount = 0;
			int femaleCount = 0;

			while (reader.ready()) {
				String[] splitString = reader.readLine().split(delim);
				// checks that no double tabs or spaces been entered in tped
				// file
				for (String str : splitString) {
					if (str.equals("")) {
						throw new PriorityPrunerException(
								"Problem with line " + line + " in [ " + filePath + " ]\r\n"
										+ "Ensure values are separated by a single space or tab character.");
					}
				}
				// checks that correct number of columns are provided
				if (splitString.length != 6) {
					throw new PriorityPrunerException(
							"Problem with line " + line + " in [ " + filePath + " ]\r\n"
									+ "Expecting 6 columns, but found " + splitString.length);
				}
				// checks that gender information is correct, it has to be
				// specified
				if (!splitString[4].equals("1") && !splitString[4].equals("2")) {
					throw new PriorityPrunerException(
							"Problem with line " + line + " in [ " + filePath + " ]\r\n"
									+ "Individual " + splitString[0] + " " + splitString[1] + " has invalid sex code " + splitString[4] 
								    + ". Must be either 1 for male or 2 for female.");
				}
				
				// checks that only founders are provided
				if (!splitString[2].equals("0") || !splitString[3].equals("0")) {
					throw new PriorityPrunerException(
							"Problem with line " + line + " in [ " + filePath + " ]\r\n"
									+ "Individual " + splitString[0] + " " + splitString[1] + " is a non-founder but only founders are allowed. ");
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
					if (splitString[4].equals("1")){
						maleCount++;
					}else if (splitString[4].equals("2")){
						femaleCount++;
					}
					line++;
				} else {
					reader.close();
					throw new PriorityPrunerException(
							"Duplicate individual found: [ " + splitString[0] + " " + splitString[1] + " ]");
				}
			}
			LogWriter.getLogger().info(individuals.size() + " individuals read from from [ " + filePath + " ]");
			LogWriter.getLogger().info(maleCount + " males, " + femaleCount + " females, and 0 of unspecified sex");
			
		} catch (FileNotFoundException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} catch (PriorityPrunerException e){
			throw e;
		} finally{
			try {
				if (reader != null){
					reader.close();
				}
			} catch (IOException e) {
			}
		}
	}
}
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
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;


/**
 * This class handles both the parsing of the SNP input file, as well as the
 * comparison between it and the TPlink data.
 */
public class SnpListFile {

	// path of the file
	private String filePath;
	
	// number of metrics defined in command line
	private int numMetrics;
	
	// ArrayList of SnpInfo objects contained in the file, in the same order as the file
	private ArrayList<SnpInfo> snps = new ArrayList<SnpInfo>();
	
	// ArrayList of SnpInfo objects contained in the file, sorted by chr, pos
	private ArrayList<SnpInfo> snpsSortedByChrPos = new ArrayList<SnpInfo>();
	
	// HashSet of chromosomes contained in the file
	private HashSet<String> chromsomeHash = new HashSet<String>();
	
	// Hashtable with SnpInfo-objects, which allows easy look up
	private Hashtable<String, ArrayList<SnpInfo>> snpNameToSnpInfoListHash = 
			new Hashtable<String, ArrayList<SnpInfo>>();
	
	// file reader
	//private BufferedReader reader;
	
	// allows access to options entered in the command line
	private CommandLineOptions options = null;
	
	// stores information from SNP input file about which SNPs that should be
	// force included
	//private HashSet<String> forceIncludeHashSet = new HashSet<String>();
	
	// splitting regex for input files that allows single tabs or spaces
	private String delim = "[\\s|\\t]";

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
	public SnpListFile(String filePath, int numMetrics, CommandLineOptions options)
			throws PriorityPrunerException {
		this.filePath = filePath;
		this.numMetrics = numMetrics;
		this.options = options;
		parseFile();
//		if (options.getForceIncludeFilePath() != null) {
//			parseForceInclude();
//			setForceInclude();
//		}
		sortSnps();
	}

	// public access to private fields 
	public ArrayList<SnpInfo> getSnps() {
		return snps;
	}

	public ArrayList<SnpInfo> getSnpsSortedByChrPos() {
		return snpsSortedByChrPos;
	}


	public String getFilePath(){
		return this.filePath;
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
		
		BufferedReader reader = null;
		int lineNum = 1;
		String line;
		
		try {
			reader = new BufferedReader(new FileReader(filePath));
			String[] header = reader.readLine().split(delim);

			// checks that necessary columns are provided
			if (header.length < 8 || !header[0].equals("name")
					|| !header[1].equals("chr")
					|| !header[2].equals("pos")
					|| !header[3].equals("a1")
					|| !header[4].equals("a2")
					|| !header[5].equals("p")
					|| !header[6].equals("forceSelect")
					|| !header[7].equals("designScore")) {
				throw new PriorityPrunerException(
						"The first 8 columns in the SNP Input Table should " +
						" equal: \"name\", \"chr\", \"pos\", \"a1\", \"a2\", \"p\", " + 
						"\"forceSelect\", \"designScore\", separated by tabs " + 
						" or spaces. Note that column names are case-sensitive.");
			}

			ArrayList<MetricNamePos> metricNames = new ArrayList<MetricNamePos>();
			
			// get metric names listed in the SNP Input Table, and makes sure no column name
			//  is duplicated
			Hashtable<String, Integer> columnNamesToIndexHash = new Hashtable<String, Integer>();
			Hashtable<String, Integer> metricNamesToIndexHash = new Hashtable<String, Integer>();
			for (int i = 0; i < header.length; i++) {
				String columnName = header[i];
				if (!columnNamesToIndexHash.containsKey(columnName)){
					columnNamesToIndexHash.put(columnName, i);
				}else{
					throw new PriorityPrunerException("The column name '" + columnName + "' " +
						" is duplicated in the SNP Input Table.");
				}
				if (i >= 8){
					metricNamesToIndexHash.put(columnName, i);
				}
			}
			
			// loops through all Metric-objects (the metrics specified in
			// the command line) and saves their name and position in the SNP
			// input file, to metricNames.
			// If gone through whole header in SNP input file without finding
			// current metric, an exception will get thrown
			for (int j = 0; j < this.options.getMetrics().size(); j++) {
				String metric =  this.options.getMetrics().get(j).getName();
				if (columnNamesToIndexHash.containsKey(metric)){
					int index = columnNamesToIndexHash.get(metric);
					metricNames.add(new MetricNamePos(metric, index));
				}else{
					throw new PriorityPrunerException(
						"The metric column '" +  this.options.getMetrics().get(j).getName() +
					"' does not exist in the SNP Input Table. Note that column names are case-sensitive.");
				}
			}

			// reads and parses every line in the SNP input file, and checks and
			// stores this information
			while ((line = reader.readLine()) != null) {
				String[] splitString = line.split(delim);
				lineNum += 1;
				
				// check that there are the expected number of columns
				if (splitString.length != header.length) {
					throw new PriorityPrunerException(
						"On line " + lineNum + " of " + filePath + ", expected " + header.length +
						" columns, but found " + splitString.length + ".");
				} 

				// parse snpname
				String snpName = new String(splitString[0]);
				if (snpName.length() == 0){
					throw new PriorityPrunerException("Invalid name on line " + lineNum +
							" in SNP Input Table.");
				}
				
				// parse chr
				String chr = new String(splitString[1]);
				chr = chr.toUpperCase();
				if (chr.length() == 0){
					throw new PriorityPrunerException("Invalid chr on line " + lineNum +
							" in SNP Input Table.");
				}
				// checks that no unsupported chromosomes are specified in the SNP
				// input file
				if (chr.equals("M")
						|| chr.equals("MT")
						|| chr.equals("CHRM")
						|| chr.equals("Y")
						|| chr.equals("CHRY")
						|| chr.equals("24")) {
					throw new PriorityPrunerException(
							"Chromosome " + chr + " not supported. Found on line " + lineNum + ".");
				}
				
				// in case chromosome X is encountered it gets recoded to "23"
//				if (chr.equals("X")
//						|| chr.equals("CHRX")
//						|| chr.equals("23")) {
//					chr = "23";
//				}
				
				// stores parsed chromosomes
				if (!chromsomeHash.contains(chr)) {
					chromsomeHash.add(chr);
				}
				
				// skip the line if this isn't the chromosome we're filtering on
				if (this.options.getChr() != null &&
						!this.options.getChr().equals(chr)){
					continue;
				}
				
				// parse pos
				int pos;
				try {
					pos = Integer.parseInt(new String(splitString[2]));
				} catch (NumberFormatException e) {
					throw new PriorityPrunerException(
						"Invalid pos at line "
									+ lineNum + " in SNP Input Table. Integer value expected.");
				}

				// parse allele1
				String allele1 = new String(splitString[3].toUpperCase());
				if (allele1.length() == 0){
					throw new PriorityPrunerException("Invalid allele1 on line " + lineNum +
							" in SNP Input Table.");
				}
				
				// parse allele2
				String allele2 = new String(splitString[4].toUpperCase());
				if (allele2.length() == 0){
					throw new PriorityPrunerException("Invalid allele2 on line " + lineNum +
							" in SNP Input Table.");
				}
				
				// parse p-value
				double pValue;
				try {
					pValue = Double.parseDouble(new String(splitString[5]));
				} catch (NumberFormatException e) {
					throw new PriorityPrunerException(
							"Invalid p at line "
									+ lineNum + " in SNP Input Table. Decimal value expected.");
				}

//				// parse numAssays
//				int numAssays;
//				try {
//					numAssays = Integer.parseInt(new String(splitString[6]));
//					if (numAssays != 1 && numAssays != 2) {
//						throw new NumberFormatException();
//					}
//				} catch (NumberFormatException e) {
//					throw new PriorityPrunerException(
//							"Invalid num_assays at line "
//									+ lineNum + " in SNP Input Table. '1' or '2' expected.");
//				}

				// parse force include flag
				boolean forceInclude;
				if (splitString[6].equals("0")) {
					forceInclude = false;
				} else if (splitString[6].equals("1")) {
					forceInclude = true;
				} else {
					throw new PriorityPrunerException(
							"Invalid force_include at line "
								+ lineNum + " in SNP Input Table. '0' or '1' expected.");
				}

				// parse the design score
				double designScore;
				try {
					designScore = Double.parseDouble(new String(splitString[7]));
				} catch (NumberFormatException e) {
					throw new PriorityPrunerException(
							"Invalid design_score at line "
								+ lineNum + " in SNP Input Table. Decimal value expected.");
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
							"Invalid metric at line "
								+ lineNum + " in SNP Input Table. Decimal value expected.");
				}

				// creates new SnpInfo-object with the parsed info from one
				// line in SNP input file
				SnpInfo snp = new SnpInfo(snpName, chr, pos, allele1,
						allele2, pValue, forceInclude, designScore, metrics);

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
						throw new PriorityPrunerException(
								"Duplicated SNP \""
										+ snpName
										+ "\" at line "
										+ lineNum
										+ " in SNP Input Table. " 
										+ "The combination of snpname, chr, pos, allele1/allele2 must be unique.");
					}
				}
				// if not a duplicate, add the SNP to following data
				// structures:
				snpInfoList.add(snp);
				snpNameToSnpInfoListHash.put(snpName, snpInfoList);
				snps.add(snp);

			}
			if (this.options.getChr() != null
					&& !chromsomeHash.contains(this.options.getChr().toUpperCase())) {
				throw new PriorityPrunerException(
						"The SNP Input Table does not contain any SNPs at chromosome "
								+ this.options.getChr()
								+ " No SNPs to analyze.");
			}
		} catch (FileNotFoundException e) {
			throw new PriorityPrunerException("Could not open file for reading:\n\n"
					+ e.getMessage());
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file for reading:\n\n"
					+ e.getMessage());
		} catch (PriorityPrunerException e){
			throw e;
		}finally{
			if (reader != null){
				try{
					reader.close();
				}catch(IOException e){	
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 * Parses the force include file, in case that option is chosen in command
	 * line.
	 * 
	 * @throws PriorityPrunerException
	 *             if problems are encountered during parsing
	 */
//	private void parseForceInclude() throws PriorityPrunerException {
//		BufferedReader reader;
//		try {
//			reader = new BufferedReader(new FileReader(
//					options.getForceIncludeFilePath()));
//			int index = 1;
//			while (reader.ready()) {
//				String[] splitString = reader.readLine().split(delim);
//				// expecting one SNP per line, comments are allowed and won't be
//				// parsed
//				if (splitString.length < 1) {
//					reader.close();
//					throw new PriorityPrunerException(
//							"Missing SNP name at line: "
//									+ (index)
//									+ " in file \""
//									+ filePath
//									+ "\". \nPlease update file and rerun program.");
//				}
//				String snpName = new String(splitString[0]);
//
//				forceIncludeHashSet.add(snpName);
//				index++;
//			}
//			reader.close();
//		} catch (IOException e) {
//			throw new PriorityPrunerException("Could not open file: "
//					+ e.getMessage());
//		}
//	}

//	/**
//	 * Sets the "force include"-flag in each SnpInfo-object according to the
//	 * force include file.
//	 */
//	private void setForceInclude() {
//		for (SnpInfo snp : snps) {
//			if (forceIncludeHashSet.contains(snp.getSnpName())) {
//				snp.setForceInclude(true);
//			}
//		}
//	}

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


	
	
	/**
	 * Sorts the parsed SNPs according to the parameters defined by the user in
	 * command line.
	 */
	private void sortSnps() {
		// sort by force include then p-value, if a specific chromosome is
		// selected
		if (options.getChr() != null){
			if (options.getSortByForceIncludeAndPValue()) {
				Collections.sort(snps, new ForceIncludePValueSorter());
				// sort by p-value then name, if a specific chromosome is selected
			} else if (!options.getSortByForceIncludeAndPValue()) {
				Collections.sort(snps);
			}
		}else{
			// sort by force include, then by chromosome and p-value, if all
			// chromosomes are selected
			if (options.getSortByForceIncludeAndPValue()) {
				Collections.sort(snps, new FIChrPValueSorter());
				// sort by chromosome and then p-value, if all chromosomes are
				// selected
			} else if (!options.getSortByForceIncludeAndPValue()) {
				Collections.sort(snps, new ChrPValueSorter());
			}
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
	 * This sorter sorts SNPs by chromosome and base pair position. (It is being
	 * used by SnpListFile when all chromosomes are included in the run.)
	 */
	private class ChrPValueSorter implements Comparator<SnpInfo> {

		@Override
		public int compare(SnpInfo x, SnpInfo y) {
			if (x.getChr().equals(y.getChr())) {
				if ((Double.compare(x.getPValue(), y.getPValue())) == 0) {
					// the name check was just added for comparison reasons
					return (x.getSnpName().compareTo(y.getSnpName()));
				} else {
					return Double.compare(x.getPValue(), y.getPValue());
				}
			} else {
				// first trying to parse the chromosomes as integers, otherwise they'll
				// get compared as Strings
				int xInt;
				int yInt;
				try {
					xInt = Integer.parseInt(x.getChr());
					yInt = Integer.parseInt(y.getChr());
					return Integer.valueOf(xInt).compareTo(Integer.valueOf(yInt));
				} catch (NumberFormatException e) {
					return x.getChr().compareTo(y.getChr());
				}
			}
		}
	}
	
	/**
	 * This sorter sorts SNPs in regards to if they are force included, then by
	 * chromosome and position. (It is being used by SnpListFile when all
	 * chromosomes are included in the run.)
	 */
	private class FIChrPValueSorter implements Comparator<SnpInfo> {

		@Override
		public int compare(SnpInfo x, SnpInfo y) {
			if (x.getForceInclude() == y.getForceInclude()) {
				if (x.getChr().equals(y.getChr())) {
					if ((Double.compare(x.getPValue(), y.getPValue())) == 0) {
						// the name check was just added for comparison reasons
						return (x.getSnpName().compareTo(y.getSnpName()));
					} else {
						return Double.compare(x.getPValue(), y.getPValue());
					}
				} else {
					// tries to parse the chromosomes as integers, otherwise they'll
					// get compared as Strings
					int xInt;
					int yInt;
					try {
						xInt = Integer.parseInt(x.getChr());
						yInt = Integer.parseInt(y.getChr());
						return Integer.valueOf(xInt).compareTo(
								Integer.valueOf(yInt));
					} catch (NumberFormatException e) {
						return x.getChr().compareTo(y.getChr());
					}
				}
			} else {
				return Boolean.valueOf(y.getForceInclude()).compareTo(
						Boolean.valueOf(x.getForceInclude()));
			}
		}
	}
	

	/**
	 * Sorts SNPs in regards to if they are force included and then by p-value. (It
	 * is used by SnpListFile when evaluating single chromosomes.)
	 */
	private class ForceIncludePValueSorter implements Comparator<SnpInfo> {

		@Override
		public int compare(SnpInfo x, SnpInfo y) {
			if (x.getForceInclude() == y.getForceInclude()) {
				if ((Double.compare(x.getPValue(), y.getPValue())) == 0) {
					// the name check was just added for comparison reasons
					return (x.getSnpName().compareTo(y.getSnpName())); 
				} else {
					return Double.compare(x.getPValue(), y.getPValue());
				}
			} else {
				return Boolean.valueOf(y.getForceInclude()).compareTo(
						Boolean.valueOf(x.getForceInclude()));
			}
		}
	}
	
	/**
	 * Sorts SnpInfo-objects in regards to chromosome and position. (Used by
	 * SnpListFile.)
	 */
	private class PosSorter implements Comparator<SnpInfo> {

		@Override
		public int compare(SnpInfo x, SnpInfo y) {
			if (x.getChr().equals(y.getChr())) {
				return Integer.valueOf(x.getPos()).compareTo(
						Integer.valueOf(y.getPos()));
			} else {
				// tries to parse the chromsomes as integers, otherwise they'll
				// get compared as Strings
				int xInt;
				int yInt;
				try {
					xInt = Integer.parseInt(x.getChr());
					yInt = Integer.parseInt(y.getChr());
					return Integer.valueOf(xInt).compareTo(Integer.valueOf(yInt));
				} catch (NumberFormatException e) {
					return x.getChr().compareTo(y.getChr());
				}
			}
		}
	}
}
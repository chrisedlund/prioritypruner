/**
 * 
 */
package edu.usc.scrc.PriorityPruner;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

/**
 * @author cedlund
 * This class parses a PLINK-formatted sample list file. This is useful for the
 * user in case they want to keep (filter in) or exclude (filter out) a set of 
 * samples from their analysis. See 
 * http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#keep. 
 * 
 * The file should be two columns (family ID, individual ID) separated by any 
 * number of white space delimiters.
 * In order to be consistent with the behavior of PLINK, error checking is minimal.
 * That means the samples in this file are not required to be present in the main 
 * genotype data file. Duplicate entries are allowed. An empty file is allowed.
 * More than two columns are allowed (only the first two are parsed, the rest 
 * are ignored). 
 */
public class PlinkSampleListFile {
	
	// path of the file to parse
	private String filePath = null;
	
	// hash set of the family ID and individual ID concatenated with a space delimiter
	private HashSet<String> sampleHashSet = new HashSet<String>();
	
	/***
	 * Constructor 
	 * @param filePath 
	 *    Path of the sample list file
	 * @throws PriorityPrunerException 
	 */
	public PlinkSampleListFile (String filePath) throws PriorityPrunerException{
		this.filePath = filePath;
		parseSampleListFile();
	}
	
	/**
	 * Sets the boolean value to true for the option in effect (remove, keep or
	 * keep_random). If the option is remove or keep, the method also parses the
	 * associated file.
	 * 
	 * @throws PriorityPrunerException
	 *             if the associated file isn't found
	 */
	private void parseSampleListFile() throws PriorityPrunerException {

		String filePath = "";
		BufferedReader reader = null;

		// parses the file with the individuals to remove/keep
		try {
			reader = new BufferedReader(new FileReader(this.filePath));
			int index = 1;
			while (reader.ready()) {
				
				// split the line on white space (consecutive delimiters okay)
				String[] splitString = reader.readLine().split("\\s+");

				if (splitString.length < 2) {
					throw new PriorityPrunerException(
							"Invalid number of columns specified in "
									+ this.filePath
									+ " at line "
									+ (index)
									+ ". Expected at least two columns, but found: "
									+ splitString.length + ".");
				}
				String famID = new String(splitString[0]);
				String indID = new String(splitString[1]);

				if (famID.equals("") || indID.equals("")) {
					throw new PriorityPrunerException(
							"Invalid formatting in file \""
									+ filePath
									+ "\" at line: "
									+ (index)
									+ "\nPlease check that values are separated by whitespace characters only.");
				}
				// stores family ID + individual ID together as a unique key, to
				// enable easy and accurate retrieval of these individuals
				this.sampleHashSet.add(famID + " " + indID);
				index++;
			}

		} catch (IOException e) {
			throw new PriorityPrunerException("Could not open file: "
					+ e.getMessage());
		} finally{
			
			if (reader != null){
				try {
					reader.close();
				} catch (IOException e) {
					// couldn't close the file, but just do nothing
				}
			}
		}
	}
	
	/***
	 * @return True if the sample file contains a Family ID, Individual ID pair. 
	 * Otherwise returns false.
	 */
	public boolean contains(String familyId, String individualId){
		return this.sampleHashSet.contains(familyId + " " + individualId);
	}
	
	
}

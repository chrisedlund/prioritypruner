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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * This class handles writing the main PriorityPruner results file, which 
 * contains each SNP name, chromosome, position, alleles, whether the SNP
 * is tagged by selected (kept) SNPs, whether the SNP is selected (kept)
 * the best tag SNP for the SNP (highest r^2), and the r^2 between the SNP
 * and its best tag. 
 * 
 * The file is tab-delimited and in plain text format. One header line and 
 * one line per SNP in the SNP list file inputted by the user
 */
public class ResultsFile {

	private CommandLineOptions options = null;
	
	public ResultsFile(SnpListFile snpListFile, CommandLineOptions options) throws PriorityPrunerException{
		this.options = options;
		createResultsFile(snpListFile);
	}
	
	/**
	 * Creates output file for pruning results. The file will be placed in file
	 * path specified by the user. It will overwrite already existing file with
	 * same name.
	 * 
	 * @throws PriorityPrunerException
	 *             if new file couldn't be created
	 */
	private void createResultsFile(SnpListFile snpListFile) throws PriorityPrunerException {

		BufferedWriter writer = null;
		
		//CommandLineOptions options = CommandLineOptions.getInstance();
		
		try {
			// creates output file
			//DecimalFormat df = new DecimalFormat("0.00##");
			
			File outputFile = new File(this.options.getOutputPrefix() + ".results");
			writer = new BufferedWriter(new FileWriter(outputFile));
			
			LogWriter.getLogger().info("Writing pruning results to [ " + this.options.getOutputPrefix() + ".results" + " ]");
			// writes to log and output files
			
			writer.write("name" + "\t" + "chr" + "\t" + "pos" + "\t" + "a1" + "\t" + "a2" + "\t" 
						+ "tagged" + "\t" + "selected" + "\t" + "best_tag" + "\t" + "r^2" + "\n" );
			for (SnpInfo snp : snpListFile.getSnps()) {
						
				String bestTag = "NA";
				String r2 = "NA";
				if (snp.getTaggedByList().size() > 0){
					Double bestR2 = new Double(-1);
					for (SnpR2Pair tag: snp.getTaggedByList()){
						if (tag.getrSquared() > bestR2){
							bestR2 = tag.getrSquared();
							bestTag = tag.getSnp().getSnpName();
						}
					}
					r2 = bestR2.toString();
					//primaryTag = snp.getTaggedByList().
				}
				
				if (snp.getSnpGenotypes().isValid() || snp.getForceInclude()){
					writer.write(snp.getSnpName() + "\t"
							+ snp.getChr() + "\t" 
							+ snp.getPos() + "\t"
							+ snp.getAllele1() + "\t" 
							+ snp.getAllele2() + "\t"
							+ (snp.getTagged() ? "1" : "0") + "\t" 
							+ (snp.getPicked() ? "1" : "0") + "\t"
							+ bestTag + "\t" 
							+ r2 + "\n");
				}
			}
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not create file: "
					+ e.getMessage()
					+ "\nPlease check that correct file path is provided.");
		}finally{
			if (writer != null){
				try{
					writer.close();
				}catch(IOException e){e.printStackTrace();}
			}
		}
	}
	
}

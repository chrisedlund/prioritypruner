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
import java.util.StringJoiner;

/***
 * This class writes the Linkage Disequilibrium (LD) file to disk, which
 * contains pairwise LD calculations done in the Pruner class. For every 
 * SNP selected (kept) by the Pruner, a pairwide LD calculation was performed
 * to determine the list of other SNPs tagged by the selected SNP. 
 * 
 * Thus, this file contains a header line, including the following columns:
 * index SNP name (the selected/kept SNP), index SNP chromosome, index SNP
 * position, index SNP allele 1, index SNP allele 2, the partner (tagged) SNP name,
 * partner SNP chromosome, partner SNP position, partner SNP allele 1, partner SNP
 * allele 2, pairwise r^2, and pairwise D'. 
 * 
 * Because this file can be quite large, we don't want to store all the LD
 * information in memory before writing the file. For that reason, this class
 * exposes a method to write LD calculations, which can be called by the Pruner 
 * class. 
 * 
 * @author Chris Edlund
 *
 */
public class LinkageDisequilibriumFile {

	private BufferedWriter ldWriter = null;
	
	public LinkageDisequilibriumFile() throws PriorityPrunerException{
		createLdFile();
	}
	
	/***
	 * Closes this file.
	 */
	public void close(){
		if (ldWriter != null){
			try{
				ldWriter.close();
			}catch(IOException e){e.printStackTrace();}
		}
	}
	
	
	/**
	 * Creates LD file to store information calculated in SnpWorkUnit. The file
	 * will be placed in file path specified by the user. It will overwrite
	 * already existing file with same name. The file is only created and the
	 * header line written. The file is kept open.
	 * 
	 * @throws PriorityPrunerException
	 *             if new file couldn't be created
	 */
	private void createLdFile() throws PriorityPrunerException {


		CommandLineOptions options = CommandLineOptions.getInstance();
		
		try {
			// creates output file
			LogWriter.getLogger().info("Writing LD metrics to [ " + options.getOutputPrefix() + ".ld" + " ]");
			File outputFile = new File(options.getOutputPrefix() + ".ld");
			this.ldWriter = new BufferedWriter(new FileWriter(outputFile));
			
			ldWriter.write("index_snp_name" + "\t" + "index_snp_chr" + "\t" + "index_snp_pos" + "\t"
			+ "index_snp_a1" + "\t" + "index_snp_a2" + "\t" + "partner_snp_name" + "\t" + 
					"partner_snp_chr" + "\t" + "partner_snp_pos" + "\t" + "partner_snp_a1" + "\t"
			+ "partner_snp_a2" + "\t" + "r^2" + "\t" + "D'" + "\n");
			
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not create file: "
					+ e.getMessage());
			// closing of the file will be handled by the caller of this class
		} 
		
	}
	
	/***
	 * Writes one row of LD data to the file
	 * @throws IOException 
	 */
	public void writeLdRow(String indexSnpName, String indexChr, int indexPos, String indexAllele1, String indexAllele2,
						String partnerSnpName, String partnerChr, int partnerPos, String partnerAllele1, String partnerAllele2,
						double partnerR2, double partnerDPrime) throws IOException{
		StringJoiner joiner = new StringJoiner("\t");
		joiner.add(indexSnpName);
		joiner.add(indexChr);
		joiner.add(Integer.toString(indexPos));
		joiner.add(indexAllele1);
		joiner.add(indexAllele2);
		joiner.add(partnerSnpName);
		joiner.add(partnerChr);
		joiner.add(Integer.toString(partnerPos));
		joiner.add(partnerAllele1);
		joiner.add(partnerAllele2);
		joiner.add(Double.toString(partnerR2));
		joiner.add(Double.toString(partnerDPrime));
		
		this.ldWriter.write(joiner.toString() + "\n");
	}
	
}

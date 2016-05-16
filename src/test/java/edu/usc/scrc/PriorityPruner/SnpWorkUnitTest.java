package edu.usc.scrc.PriorityPruner;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/***
 * This class tests LD calculation of the SnpWorkUnit class. The SnpWorkUnit class takes an array 
 * of compressed genotypes (SnpGenotypes), an index snp, and calculates LD between the index SNP 
 * and all the other snps in the array. Tt also takes a keptFounders array, which is an array of 
 * Individual objects in the same order as the genotypes in each SnpGenotypes object. This only
 * matters for calculating LD in chrX, since sex of each individual matters.
 * 
 * Thus the purpose of this test is to make sure SnpWorkUnit is 
 *  1) Calculating pairwise r2 correctly
 *  2) Calculating D' correctly
 *  3) Returning the expected number of results (# of polymorphic markers)
 * 
 * @author cedlund
 *
 */
// the SnpWorkUnit class takes an array of compressed genotypes (SnpGenotypes), an index snp,
// and calculates LD between the index SNP and all the other snps in the array

//  it also takes a keptFounders array, which is an array of Individual objects
//   in the same order as the genotypes in each SnpGenotypes object. this only matters for
//  calculating LD in chrX, since sex of each individual matters

// we can do all unit tests by parsing a transposed PLINK data set using the TPlink class
//  which will contain SNPs for all test cases:
//   1) 1 monomorphic snp vs 1 polymorphic snp
//   2) 1 monomorphic snp vs 1 monomorphic snp
//   3) chr X SNPs
//   4) polymorphic diploid snps
public class SnpWorkUnitTest {


	@Before
	public void setUp() throws Exception {
		
	
	}

	@After
	public void tearDown() throws Exception {
	}

	/***
	 * This tests the LD calculations for a 50-kb region on chromosom 12 using genotypes
	 * from 1000 Genomes phase 3, YRI population. It loads a transposed PLINK file, SNP list file,
	 * and expected output generated using Haploview. For every pair of SNPs in the region (MAF>0),
	 * this verifies that r^2 and D' are within 0.005 of the true value from Haploview. 
	 *  
	 * @throws PriorityPrunerException
	 * @throws IOException
	 */
	@Test
	public void testDiploidLdCalculations() throws PriorityPrunerException, IOException {
		
		// set up default options
		CommandLineOptions options = new CommandLineOptions();
		
		// get test resource file paths
		ClassLoader classLoader = getClass().getClassLoader();
		String filePathTPedDiploid = classLoader.getResource("pp_1kgp3_yri_chr12_test.tped").getPath();
		String filePathTFamDiploid = classLoader.getResource("pp_1kgp3_yri_chr12_test.tfam").getPath();
		String snpListFilePathDiploid = classLoader.getResource("pp_1kgp3_yri_chr12_test.snp_input.txt").getPath();
		String haploviewLdFilePathDiploid = classLoader.getResource("pp_1kgp3_yri_chr12.haploview.ld").getPath();
		
		// create a new SnpListFile (only polymorphic SNPs)
		SnpListFile snpListFileDiploid = new SnpListFile(snpListFilePathDiploid, 0, options);
		
		// create a new Genotypes file by parsing transposed PLINK data
		Genotypes genotypesDiploid = new TPlink(filePathTPedDiploid, filePathTFamDiploid, snpListFileDiploid, null, options);
		
		// we need to first calculate MAF and missingness so that each SNP gets it's alleles assigned
		for (SnpGenotypes g : genotypesDiploid.getSnpGenotypes()) {
			g.calculateMafHweMissingPercentCompressed(genotypesDiploid.getKeptFounders());
		}
		
		// parse expected LD results from Haploview
		HashMap<String,Result> haploviewResultsDiploid  = parseHaplotypeLdFile(haploviewLdFilePathDiploid);
		
		// run validation tests
		validateR2AndDPrime(snpListFileDiploid, genotypesDiploid, haploviewResultsDiploid);
	}
	
	/***
	 * This tests the LD calculations for a 50-kb region on chromosome X using genotypes
	 * from 1000 Genomes phase 3, YRI population. It loads a transposed PLINK file, SNP list file,
	 * and expected output generated using Haploview. For every pair of SNPs in the region (MAF>0),
	 * this verifies that r^2 and D' are within 0.005 of the true value from Haploview. 
	 *  
	 * @throws PriorityPrunerException
	 * @throws IOException
	 */
	@Test
	public void testHaploidLdCalculations() throws PriorityPrunerException, IOException{
		
		// set up default options
		CommandLineOptions options = new CommandLineOptions();
		
		// get test resource file paths
		ClassLoader classLoader = getClass().getClassLoader();
		String filePathTPedHaploid = classLoader.getResource("pp_1kgp3_yri_chrX_test.tped").getPath();
		String filePathTFamHaploid = classLoader.getResource("pp_1kgp3_yri_chrX_test.tfam").getPath();
		String snpListFilePathHaploid = classLoader.getResource("pp_1kgp3_yri_chrX_test.snp_input.txt").getPath();
		String haploviewLdFilePathHaploid = classLoader.getResource("pp_1kgp3_yri_chrX.haploview.ld").getPath();		
		
		// create a new SnpListFile (only polymorphic SNPs)
		SnpListFile snpListFileHaploid = new SnpListFile(snpListFilePathHaploid, 0, options);
		
		// create a new Genotypes file by parsing transposed PLINK data
		Genotypes genotypesHaploid = new TPlink(filePathTPedHaploid, filePathTFamHaploid, snpListFileHaploid, null, options);
		
		// we need to first calculate MAF and missingness so that each SNP gets it's alleles assigned
		for (SnpGenotypes g : genotypesHaploid.getSnpGenotypes()) {
			g.calculateMafHweMissingPercentCompressed(genotypesHaploid.getKeptFounders());
		}
		
		// parse expected LD results from Haploview
		HashMap<String,Result> haploviewResultsHaploid = parseHaplotypeLdFile(haploviewLdFilePathHaploid);
		
		// run validation tests
		validateR2AndDPrime(snpListFileHaploid, genotypesHaploid, haploviewResultsHaploid);
	}
	
	/***
	 * This is a helper function method to loop through every pair of SNPs in a SnpListFile with 
	 * corresponding Genotypes and Haploview results to compare r2 and D'
	 * @param snpListFile
	 * @param genotypes
	 * @param haploviewResults
	 * @throws PriorityPrunerException
	 */
	private void validateR2AndDPrime(SnpListFile snpListFile, Genotypes genotypes, HashMap<String,Result> haploviewResults) throws PriorityPrunerException{
		
		// set up genotypes array of all SNPs
		ArrayList<SnpGenotypes> genotypesList = new ArrayList<SnpGenotypes>();
		for (SnpInfo snpInfo: snpListFile.getSnps()){
			if (snpInfo.getSnpGenotypes().isValid()){
				genotypesList.add(snpInfo.getSnpGenotypes());
			}
		}
		
		
		// for each snp, calculate LD with all other SNPs
		int numValidated = 0;
		for (int i = 0; i<  snpListFile.getSnps().size(); i++){
			SnpInfo snpInfo = snpListFile.getSnps().get(i);
			if (snpInfo.getSnpGenotypes().isValid()){
				SnpWorkUnit snpWorkUnit = new SnpWorkUnit(snpInfo.getSnpName(), genotypesList, i, 
						genotypes.getIndividuals());
				snpWorkUnit.performWork();
				
				// verify the LD calculations and number of results returned
				// we know that all SNPs in the array are within 500kb (the default distance) and all SNPs are polymorphic, so the expected
				//  number of results is equal to the number of genotypes
				numValidated += verifyLdCalculations(snpWorkUnit, haploviewResults, genotypesList.size());
				
			}
		}
		System.out.println(numValidated + " comparisons validated.");
				
	}
	
	/***
	 * This is a helper function to compare the LD calculations from one SnpWorkUnit with LD calculations
	 * from Haploview. 
	 * @param snpWorkUnit
	 * @param haploviewResults
	 * @return
	 */
	private int verifyLdCalculations(SnpWorkUnit snpWorkUnit, HashMap<String,Result> haploviewResults, int expectedNumResults){
		for (Result snpWorkUnitResult: snpWorkUnit.getResults()){
			Result haploviewResult = haploviewResults.get(
					snpWorkUnitResult.getSnpName() + " "  + snpWorkUnitResult.getPartnerSnpName());
			if (haploviewResult == null && snpWorkUnitResult.getSnpName().equals(snpWorkUnitResult.getPartnerSnpName())){
				haploviewResult = new Result(snpWorkUnitResult.getSnpName(), 0, "", 0, snpWorkUnitResult.getPartnerSnpName(), 0, "", 0, 1, 1, null, null);
			}
			//System.out.println(snpWorkUnitResult.getSnpName() + " " + snpWorkUnitResult.getPartnerSnpName() + " " + String.format("%.3f", snpWorkUnitResult.getRSquared()) + " " + haploviewResult.getRSquared());
			
			// assert r-squared difference between PriorityPruner and Haploview is less than 0.005
			assert(Math.abs(snpWorkUnitResult.getRSquared() - haploviewResult.getRSquared()) < 0.005 );
			
			// assert d-prime difference between PriorityPruner and Haploview is less than 0.005
			assert(Math.abs(snpWorkUnitResult.getDPrime() - haploviewResult.getDPrime()) < 0.005 );
		}
		
		//verify we got the correct number of results back
		assertEquals(snpWorkUnit.getResults().size(), expectedNumResults);
		
		return snpWorkUnit.getResults().size();
	}
	
	/***
	 * This method parses and stores Haploview LD calculations into a HashMap. Since Haploview only
	 * outputs each pair only once, we store SNP pair key twice for convenience. Once as "snp1 snp2", 
	 * then again as "snp2 snp1".
	 * @param filePath
	 * @return
	 * @throws IOException
	 */
	private HashMap<String,Result> parseHaplotypeLdFile(String filePath) throws IOException{
		String delim = "[\\s|\\t]";
		BufferedReader reader = new BufferedReader(new FileReader(filePath));
		HashMap<String,Result> snpNamePairToResultHash = new HashMap<String,Result>();
		// skip header line
		reader.readLine();
		
		// read each line and store results in hashmap for both "snp1 snp2" and "snp2 snp1"
		while (reader.ready()) {
			String[] splitString = reader.readLine().split(delim);
			String snpName1 = splitString[0];
			String snpName2 = splitString[1];
			double dprime = Double.parseDouble(splitString[2]);
			double r2 = Double.parseDouble(splitString[4]);
			snpNamePairToResultHash.put(snpName1 + " " + snpName2, 
					new Result(snpName1, 0, "", 0, snpName2, 0, "", 0, r2, dprime, null, null));
			if (snpName1 != snpName2){
				snpNamePairToResultHash.put(snpName2 + " " + snpName1, 
						new Result(snpName2, 0, "", 0, snpName1, 0, "", 0, r2, dprime, null, null));				
			}
		}
		reader.close();
		return snpNamePairToResultHash;
	}
}

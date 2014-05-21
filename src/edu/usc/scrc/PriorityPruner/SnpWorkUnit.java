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

import java.util.ArrayList;

/**
 * This class is used to find all the SNPs (in the given search range) that are
 * in LD. The code for calculating LD is adapted from Haploview (Barrett et
 * al.).
 * 
 * Barrett JC, Fry B, Maller J, Daly MJ. Haploview: analysis and visualization of LD and haplotype maps. 
 * Bioinformatics. 2005 Jan 15 [PubMed ID: 15297300]
 */
public class SnpWorkUnit {

	// genotypes (AB variety)
	private final int AA = 0;
	private final int AB = 1;
	private final int BA = 2;
	private final int BB = 3;

	// index in currentGenotypes for this SNP
	private int referenceSNPIndex;

	// LD results for this SNP
	private ArrayList<Result> results = new ArrayList<Result>();

	private String snpName;
	private SnpInfo snpInfo;
	private ArrayList<SnpGenotypes> currentGenotypes = new ArrayList<SnpGenotypes>();
	//private ArrayList<Integer> founderIndices = new ArrayList<Integer>();
	//private ArrayList<String> subjectSexes = new ArrayList<String>();
	private ArrayList<Individual> keptFounders;
	private double minMaf;
	private double minimumHardyWeinbergPvalue;
	private double minimumGenotypePercentage;
	private boolean indexSnpPassed = true;

	// variables for r^2 calculation
	private double[] known = new double[5];
	private double[] numHaps = new double[4];
	private double[] probHaps = new double[4];
	private double const_prob = -1;
	private int unknownDH = -1;
	private final double LN10 = Math.log(10);
	private final double TOLERANCE = 0.00000001;

	/**
	 * Constructor for SnpWorkUnit. Initializes instance variables and initiates
	 * conversion of genotypes from compressed format to the original format
	 * used for LD calculation.
	 * 
	 * @param snpName
	 *            name of index SNP
	 * @param genotypesList
	 *            list of SnpGenotypes-object within the current pruning window
	 * @param referenceSNPIndex
	 *            the position of the index SNP in genotypesList
	 * @param keptFounders
	 *            founder subjects in the genotypes list. Right now
	 *            we only support founders, so at the moment this corresponds to
	 *            an index over all subjects chosen in this run
	 * @param minMaf
	 *            minimum allele frequency, as defined in command line
	 * @param minimumHardyWeinbergPvalue
	 *            minimum Hardy-Weinberg p-value, as defined in command line
	 * @param minimumGenotypePercentage
	 *            minimum genotype percentage, as defined in command line
	 * @throws PriorityPrunerException
	 *             if problem are encountered during conversion (decompression)
	 *             of genotypes
	 */
	public SnpWorkUnit(String snpName, ArrayList<SnpGenotypes> genotypesList,
			int referenceSNPIndex, ArrayList<Individual> keptFounders, double minMaf,
			double minimumHardyWeinbergPvalue, double minimumGenotypePercentage)
			throws PriorityPrunerException {

		this.snpName = snpName;
		this.currentGenotypes = genotypesList;
		this.referenceSNPIndex = referenceSNPIndex + 1;
		this.keptFounders = keptFounders;
		this.minMaf = minMaf;
		this.minimumHardyWeinbergPvalue = minimumHardyWeinbergPvalue;
		this.minimumGenotypePercentage = minimumGenotypePercentage;
		this.snpInfo = currentGenotypes.get(referenceSNPIndex).getSnpInfo();
		//convertGenotypes(currentGenotypes);
	}

	/**
	 * Initiates all the calculations necessary for this SNP work unit.
	 * 
	 * @throws PriorityPrunerException
	 *             if SNP is missing, didn't pass all filters, or problems are
	 *             encountered during the getMafHweMissingPercent and
	 *             linkageToChrom methods
	 */
	public void performWork() throws PriorityPrunerException {

		if (snpInfo == null) {
			throw new PriorityPrunerException(
					"Couldn't find index SNP \""
							+ snpName
							+ "\". Make sure that data in tped and SNP input file matches.");
		} else {
			// if calculations of maf, hwe, and missing percentage have not
			// already been done for current SNP - calculate this!
			for (SnpGenotypes genotypes : currentGenotypes) {
				if (!genotypes.getCalculationDone()) {
					//getMafHweMissingPercent(genotype, founderIndices,subjectSexes);
					getMafHweMissingPercentCompressed(genotypes);
				}
			}
			
			// add "4" to heterozygous genotypes
			//linkageToChrom();
		}

		// initiates calculations, will return false if the index SNP didn't
		// pass the user defined filters
		if (!calculateResults()) {
			indexSnpPassed = false;
		}

		// deletes the LD-formats, since they're memory consuming
		for (SnpGenotypes genotype : currentGenotypes) {
			genotype.eraseLdFormatGenotypes();
		}
	}

	/**
	 * Converts compressed genotypes stored as bytes to the original LD-format,
	 * used by the SnpWorkUnit at the moment. This might be updated in future
	 * versions, so that the SnpWorkUnit could work with the compressed format
	 * and save on memory.
	 * 
	 * @param currentGenotypes
	 *            the SNPs in the current pruning window
	 * @throws PriorityPrunerException
	 *             if an invalid genotype is encountered during conversion
	 */
/*	private void convertGenotypes(ArrayList<SnpGenotypes> currentGenotypes)
			throws PriorityPrunerException {

		// goes through all SNPs in pruning window
		for (SnpGenotypes currentSnpGenotype : currentGenotypes) {

			byte[] byteArray = currentSnpGenotype.getGenotypes();
			int length = byteArray.length;
			byte[][] genotypesInt = new byte[length * 4][2];
			int index = 0;

			for (byte b : byteArray) {
				int a = (int) (b < 0 ? 256 + b : b);
				// each byte stores 4 genotypes
				for (int i = 0; i < 4; i++) {
					if ((192 & a) == 0) {
						genotypesInt[index][0] = 0;
						genotypesInt[index][1] = 0;
						a = a << 2;
						index++;

					} else if ((192 & a) == 64) {
						genotypesInt[index][0] = 1;
						genotypesInt[index][1] = 1;
						a = a << 2;
						index++;

					} else if ((192 & a) == 128) {
						genotypesInt[index][0] = 2;
						genotypesInt[index][1] = 2;
						a = a << 2;
						index++;

					} else if ((192 & a) == 192) {
						genotypesInt[index][0] = 1;
						genotypesInt[index][1] = 2;
						a = a << 2;
						index++;
					} else {
						throw new PriorityPrunerException(
								"Invalid genotype \""
										+ currentSnpGenotype.getAllele1()
										+ currentSnpGenotype.getAllele2()
										+ "\" found in SNP \""
										+ currentSnpGenotype.getSnpName()
										+ "\". \nPlease update information about this SNP in the input files and rerun program.");
					}
				}
			}
			// set the LD-format in SnpGenotypes-object
			currentSnpGenotype.setLdFormatGenotypes(genotypesInt);
		}
	}
*/
	/**
	 * Method ported from Haploview for LD calculation. Adds "4" to each allele
	 * in a heterozygous genotype stored in the LD-format.
	 */
/*	private void linkageToChrom() throws PriorityPrunerException {

		// for each individual, and each SNP in pruning window,
		// check if alleles are either homozygous or if one of them are
		// missing, else add "4".
		for (int i = 0; i < keptFounders.size(); i++) {
			for (int j = 0; j < currentGenotypes.size(); j++) {
				byte[][] geno = currentGenotypes.get(j).getLdFormatGenotypes();
				byte thisMarkerA = geno[i][0];
				byte thisMarkerB = geno[i][1];
				if (thisMarkerA == thisMarkerB || thisMarkerA == 0
						|| thisMarkerB == 0) {
				} else {
					// if heterozygous, add "4"
					geno[i][0] = (byte) (thisMarkerA + 4);
					geno[i][1] = (byte) (thisMarkerB + 4);
				}
			}
		}
	}*/

	/**
	 * Checks if a SNP passes all of the user-defined filters.
	 * 
	 * @param genotypes
	 *            the SNP to be evaluated, stored as a SnpGenotypes-object
	 * @return flag indicating if SNP is valid - returns true if it is, false
	 *         otherwise.
	 */
	private boolean isSnpValid(SnpGenotypes genotypes) {

		if ((genotypes.getMaf() >= minMaf)
				&& (((1 - genotypes.getMissingPercent())) >= minimumGenotypePercentage)
				&& (genotypes.getHwePvalue() >= this.minimumHardyWeinbergPvalue)) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Method ported from Haploview to calculate LD. Loops through each SNP and
	 * checks that it passes the user's filters for maf, hwe, and missing
	 * genotype percentage. If the SNP passes, LD calculation is initiated and a
	 * Result-object gets created when calculation is successful. The result is
	 * added to a result list.
	 * 
	 * @return flag indicating if index SNP is valid
	 */
	private boolean calculateResults() {

		// gets genotypes of index SNP
		SnpGenotypes referenceGenotypes = currentGenotypes
				.get(referenceSNPIndex - 1);

		if (isSnpValid(referenceGenotypes)) {
			for (SnpGenotypes genotypes : currentGenotypes) {
				if (isSnpValid(genotypes)) {
					// if both index SNP and the second SNP are valid (i.e.,
					// they passed the user defined filter for maf, hwe and
					// missing genotype percentage), LD calculation for this
					// pair is initiated
					LdResult ldResult = calculateLdResultCompressed(referenceGenotypes,
							genotypes);
					if (ldResult != null && ldResult.getRSquared() > 1) {
						ldResult.setRSquared(1);
					}
					// if the result from the calculation is not null and the
					// r^2-value is a valid number, the information gets stored
					// in a Result-object.
					if (ldResult != null
							&& !Double.isNaN(ldResult.getRSquared())) {
						Result result = new Result(snpName,
								referenceGenotypes.getMaf(), snpInfo.getChr(),
								snpInfo.getPos(), genotypes.getSnpName(),
								genotypes.getMaf(), genotypes.getSnpInfo()
										.getChr(), genotypes.getSnpInfo()
										.getPos(), ldResult.getRSquared(),
								ldResult.getDPrime(),
								referenceGenotypes.getSnpInfo(),
								genotypes.getSnpInfo());
						// valid results get added to a result list
						addResult(result);
					}
				}
			}
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Method ported from Haploview. Calculates maf, hwe and missing genotype
	 * percentage for the current SNP.
	 * 
	 * @param genotypes
	 *            the SNP to be evaluated, stored as a SnpGenotypes-object
	 * @param founderIndices
	 *            indices of founder subjects
	 * @param subjectSexes
	 *            gender of the subjects
	 * @throws PriorityPrunerException
	 *             if an invalid genotype is encountered in a SNP
	 */

	private void getMafHweMissingPercentCompressed(SnpGenotypes genotypes)
			throws PriorityPrunerException {
		
		//byte allele1 = 0;
		//byte allele2 = 0;
		int numAllele1 = 0;
		int numAllele2 = 0;
		int numMissing = 0;
		int numChromosomes = 0;
		int[] founderHomCount = new int[5];
		int founderHetCount = 0;

		// initialize founderHomCount array to contain all zeros
		for (int i = 0; i < 5; i++) {
			founderHomCount[i] = 0;
		}

		for (int f = 0; f < keptFounders.size(); f++) {
			// at the moment, "i" will equal "f", since we only support
			// founders. This might be updated in future versions.
			Individual founder = keptFounders.get(f);
			
			// counts genotypes for hwe (only for diploids)
			boolean haploid  = genotypes.getSnpInfo().isChrX() && founder.getSex() == Individual.Sex.MALE;
			
			byte genotype = genotypes.getByteGenotype(f);
			
			if (!haploid) {
				if (genotype != 0) {
					numChromosomes+= 2;
					if (genotype == 3) {
						founderHetCount++;
						numAllele1++;
						numAllele2++;
					} else if (genotype == 1){
						founderHomCount[1]++;
						numAllele1+=2;
					} else if (genotype == 2){
						founderHomCount[2]++;
						numAllele2+=2;
					}
				}else{
					numMissing+= 2;
				}
			}else{
				
				if (genotype != 0) {
					numChromosomes+= 1;
					if (genotype == 1){
						numAllele1+=1;
					} else if (genotype == 2){
						numAllele2+=1;
					}
				}else{
					numMissing+= 1;
				}
			}
		}

		// sets the values calculated, as well as a flag indicating that
		// calculation is completed, in the SnpGenotypes-object

		if (numAllele1 > numAllele2) {
			genotypes.setMaf((double) numAllele2
					/ (double) (numAllele1 + numAllele2));
			genotypes.setMinorAllele((byte)2);
			genotypes.setMajorAllele((byte)1);
		} else {
			genotypes.setMaf((double) numAllele1
					/ (double) (numAllele1 + numAllele2));
			genotypes.setMinorAllele((byte)1);
			genotypes.setMajorAllele((byte)2);
		}
		genotypes.setMissingPercent((double) numMissing / numChromosomes);
		genotypes.setHwePvalue(this.getHweValue(founderHomCount,
				founderHetCount));
		genotypes.setCalculationDone();
		//System.out.println(genotypes.getSnpName() + "\t" + genotypes.getMissingPercent() + "\t" + genotypes.getHwePvalue() + "\t" + genotypes.getMaf());
	}

/*
	*//**
	 * Method ported from Haploview. Calculates maf, hwe and missing genotype
	 * percentage for the current SNP.
	 * 
	 * @param genotypes
	 *            the SNP to be evaluated, stored as a SnpGenotypes-object
	 * @param founderIndices
	 *            indices of founder subjects
	 * @param subjectSexes
	 *            gender of the subjects
	 * @throws PriorityPrunerException
	 *             if an invalid genotype is encountered in a SNP
	 *//*

	private void getMafHweMissingPercent(SnpGenotypes genotypes,
			ArrayList<Integer> founderIndices, ArrayList<String> subjectSexes)
			throws PriorityPrunerException {

		byte allele1 = 0;
		byte allele2 = 0;
		int numAllele1 = 0;
		int numAllele2 = 0;
		int numMissing = 0;
		int numChromosomes = 0;
		int[] founderHomCount = new int[5];
		int founderHetCount = 0;

		// initialize founderHomCount array to contain all zeros
		for (int i = 0; i < 5; i++) {
			founderHomCount[i] = 0;
		}


		for (int f = 0; f < founderIndices.size(); f++) {
			// at the moment, "i" will equal "f", since we only support
			// founders. This might be updated in future versions.
			int i = founderIndices.get(f);

			boolean haploid = genotypes.getSnpInfo().isChrX() && subjectSexes.get(i).equals("1");
			
			// counts genotypes for hwe (only for diploids)
			if (!haploid) {

				byte a1 = genotypes.getLdFormatGenotypes()[i][0];
				byte a2 = genotypes.getLdFormatGenotypes()[i][1];

				if (a1 >= 5) {
					a1 -= 4;
				}
				if (a2 >= 5) {
					a2 -= 4;
				}
				if (a1 != 0 && a2 != 0) {
					if (a1 != a2) {
						founderHetCount++;
					} else {
						founderHomCount[a1]++;
					}
				}
			}

			for (int j = 0; j < 2; j++) {
				byte tempAllele = genotypes.getLdFormatGenotypes()[i][j];
				if (tempAllele >= 5) {
					tempAllele = (byte) (tempAllele - 4);
				}
				try {
					numChromosomes++;
					if (tempAllele == 0) {
						numMissing++;
					} else if (allele1 == 0) {
						allele1 = tempAllele;
						numAllele1++;
					} else if (allele1 == tempAllele) {
						numAllele1++;
					} else if (allele2 == 0) {
						allele2 = tempAllele;
						numAllele2++;
					} else if (allele2 == tempAllele) {
						numAllele2++;
					} else {
						throw new PriorityPrunerException(
								"Invalid genotype \""
										+ allele1
										+ allele2
										+ "\" in SNP \""
										+ genotypes.getSnpName()
										+ "\". Please update input files and rerun program.");
					}
				} catch (Exception e) {
					throw new PriorityPrunerException(
							"Invalid genotype \""
									+ allele1
									+ allele2
									+ "\" in SNP \""
									+ genotypes.getSnpName()
									+ "\". Please update input files and rerun program.");
				}

				// we only count the first allele for male haploids
				if (haploid) {
					break;
				}
			}
		}

		// sets the values calculated, as well as a flag indicating that
		// calculation is completed, in the SnpGenotypes-object

		if (numAllele1 > numAllele2) {
			genotypes.setMaf((double) numAllele2
					/ (double) (numAllele1 + numAllele2));
			genotypes.setMinorAllele(allele2);
			genotypes.setMajorAllele(allele1);
		} else {
			genotypes.setMaf((double) numAllele1
					/ (double) (numAllele1 + numAllele2));
			genotypes.setMinorAllele(allele1);
			genotypes.setMajorAllele(allele2);
		}
		genotypes.setMissingPercent((double) numMissing / numChromosomes);
		genotypes.setHwePvalue(this.getHweValue(founderHomCount,
				founderHetCount));
		genotypes.setCalculationDone();
		//System.out.println(genotypes.getSnpName() + "\t" + genotypes.getMissingPercent() + "\t" + genotypes.getHwePvalue() + "\t" + genotypes.getMaf());
	}*/

	/**
	 * Method ported from Haploview. Calculates r^2 and D' between two SNPs.
	 * 
	 * @param genotypes1
	 *            genotypes for SNP 1
	 * @param genotypes2
	 *            genotypes for SNP 2
	 * @param founderIndices
	 *            indices of founder subjects
	 * @return LdResult-object
	 */
	private LdResult calculateLdResultCompressed(SnpGenotypes genotypes1,
			SnpGenotypes genotypes2) {
		int doublehet = 0;
		int[][] twoMarkerHaplos = new int[3][3];
		int count;
		double loglike, oldloglike, rsq, num, tmp, denom, denom1, denom2, dprime;

		// if comparing a SNP with itself just return r^2=1, D'=1
		if (genotypes1 == genotypes2) {
			return new LdResult(1, 1);
		}
		// initialize twoMarkerHaplos matrix to contain all zeros
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				twoMarkerHaplos[i][j] = 0;
			}
		}

		int[] marker1num = new int[5];
		int[] marker2num = new int[5];

		marker1num[0] = 0;
		marker1num[genotypes1.getMajorAllele()] = 1;
		marker1num[genotypes1.getMinorAllele()] = 2;
		marker2num[0] = 0;
		marker2num[genotypes2.getMajorAllele()] = 1;
		marker2num[genotypes2.getMinorAllele()] = 2;

		byte a1 = 0, a2 = 0, b1 = 0, b2 = 0;

		// iterate through all chromosomes in dataset
		for (int f = 0; f < keptFounders.size(); f++) {
			// at the moment, "i" will equal "f", since we only support founders
			//int i = founderIndices.get(f);
			Individual founder = keptFounders.get(f);

			byte genotype1 = genotypes1.getByteGenotype(f);
			byte genotype2 = genotypes2.getByteGenotype(f);
			if (genotype1 == 0){
				a1 = 0;
				b1 = 0;
			}else if (genotype1 == 1){
				a1 = 1;
				b1 = 1;
			}else if (genotype1 == 2){
				a1 = 2;
				b1 = 2;
			}else if (genotype1 == 3){
				a1 = 1;
				b1 = 2;
			}
			if (genotype2 == 0){
				a2 = 0;
				b2 = 0;
			}else if (genotype2 == 1){
				a2 = 1;
				b2 = 1;
			}else if (genotype2 == 2){
				a2 = 2;
				b2 = 2;
			}else if (genotype2 == 3){
				a2 = 1;
				b2 = 2;
			}
			
			// assign alleles for each of a pair of chromosomes at a marker
			// to four variables
			boolean haploid = snpInfo.isChrX() && founder.getSex() == Individual.Sex.MALE;
			if (haploid) {

				// haploid (x chrom/male)
				if (a1 != 0 && a2 != 0) {
					twoMarkerHaplos[marker1num[a1]][marker2num[a2]]++;
				}
			} else {
				// diploid
				if (genotype1 == 0 || genotype2 == 0) {
					// skip missing data
				} else if (genotype1 == 3 && genotype2 == 3) {
					doublehet++;
					// find doublehets and resolved haplotypes
				} else if (genotype1 == 3) {
					twoMarkerHaplos[1][marker2num[a2]]++;
					twoMarkerHaplos[2][marker2num[a2]]++;
				} else if (genotype2 == 3) {
					twoMarkerHaplos[marker1num[a1]][1]++;
					twoMarkerHaplos[marker1num[a1]][2]++;
				} else {
					twoMarkerHaplos[marker1num[a1]][marker2num[a2]]++;
					twoMarkerHaplos[marker1num[b1]][marker2num[b2]]++;
				}
			}
		}

		// another monomorphic marker check
		int r1 = twoMarkerHaplos[1][1] + twoMarkerHaplos[1][2];
		int r2 = twoMarkerHaplos[2][1] + twoMarkerHaplos[2][2];
		int c1 = twoMarkerHaplos[1][1] + twoMarkerHaplos[2][1];
		int c2 = twoMarkerHaplos[1][2] + twoMarkerHaplos[2][2];

		if ((r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0) && doublehet == 0) {
			return null;
		}

		double pA1, pB1, pA2, pB2;
		int total_chroms;

		known[AA] = twoMarkerHaplos[1][1];
		known[AB] = twoMarkerHaplos[1][2];
		known[BA] = twoMarkerHaplos[2][1];
		known[BB] = twoMarkerHaplos[2][2];
		unknownDH = doublehet;

		total_chroms = (int) (known[AA] + known[AB] + known[BA] + known[BB] + 2 * unknownDH);

		pA1 = (known[AA] + known[AB] + unknownDH) / (double) total_chroms;
		pB1 = 1 - pA1;
		pA2 = (known[AA] + known[BA] + unknownDH) / (double) total_chroms;
		pB2 = 1 - pA2;
		const_prob = 0.1;

		probHaps[AA] = const_prob;
		probHaps[AB] = const_prob;
		probHaps[BA] = const_prob;
		probHaps[BB] = const_prob;

		count_haps(0);
		estimate_p();

		// now we have an initial reasonable guess at p we can
		// start the EM - let the fun begin
		const_prob = 0.0;
		count = 1;
		loglike = -999999999;

		while (count < 1000) {
			oldloglike = loglike;
			count_haps(count);
			loglike = (known[AA] * Math.log(probHaps[AA]) + known[AB]
					* Math.log(probHaps[AB]) + known[BA]
					* Math.log(probHaps[BA]) + known[BB]
					* Math.log(probHaps[BB]))
					/ LN10
					+ ((double) unknownDH * Math.log(probHaps[AA]
							* probHaps[BB] + probHaps[AB] * probHaps[BA]))
					/ LN10;
			if (Math.abs(loglike - oldloglike) < TOLERANCE) {
				break;
			}
			estimate_p();
			count++;
		}

		num = probHaps[AA] * probHaps[BB] - probHaps[AB] * probHaps[BA];

		if (num < 0) {

			// flip matrix so we get the positive D'
			// flip AA with AB and BA with BB
			tmp = probHaps[AA];
			probHaps[AA] = probHaps[AB];
			probHaps[AB] = tmp;

			tmp = probHaps[BB];
			probHaps[BB] = probHaps[BA];
			probHaps[BA] = tmp;

			// flip frequency of second allele
			// done in this slightly asinine way because of a compiler
			// bugz0r in the dec-alpha version of java
			// which causes it to try to parallelize the swapping operations
			// and mis-schedules them
			pA2 = pA2 + pB2;
			pB2 = pA2 - pB2;
			pA2 = pA2 - pB2;

			// flip counts in the same fashion as p's
			tmp = numHaps[AA];
			numHaps[AA] = numHaps[AB];
			numHaps[AB] = tmp;

			tmp = numHaps[BB];
			numHaps[BB] = numHaps[BA];
			numHaps[BA] = tmp;

			// num has now undergone a sign change
			num = probHaps[AA] * probHaps[BB] - probHaps[AB] * probHaps[BA];

			// flip known array for likelihood computation
			tmp = known[AA];
			known[AA] = known[AB];
			known[AB] = tmp;

			tmp = known[BB];
			known[BB] = known[BA];
			known[BA] = tmp;
		}

		denom1 = (probHaps[AA] + probHaps[BA]) * (probHaps[BA] + probHaps[BB]);
		denom2 = (probHaps[AA] + probHaps[AB]) * (probHaps[AB] + probHaps[BB]);

		if (denom1 < denom2) {
			denom = denom1;
		} else {
			denom = denom2;
		}
		dprime = num / denom;

		// add computation of r^2 = (D^2)/p(1-p)q(1-q)
		rsq = num * num / (pA1 * pB1 * pA2 * pB2);

		return new LdResult(rsq, dprime);
	}
	
	
	/**
	 * Method ported from Haploview. Calculates r^2 and D' between two SNPs.
	 * 
	 * @param genotypes1
	 *            genotypes for SNP 1
	 * @param genotypes2
	 *            genotypes for SNP 2
	 * @param founderIndices
	 *            indices of founder subjects
	 * @return LdResult-object
	 */
/*	private LdResult calculateLdResult(SnpGenotypes genotypes1,
			SnpGenotypes genotypes2) {
		int doublehet = 0;
		int[][] twoMarkerHaplos = new int[3][3];
		int count;
		double loglike, oldloglike, rsq, num, tmp, denom, denom1, denom2, dprime;

		// if comparing a SNP with itself just return r^2=1, D'=1
		if (genotypes1 == genotypes2) {
			return new LdResult(1, 1);
		}
		// initialize twoMarkerHaplos matrix to contain all zeros
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				twoMarkerHaplos[i][j] = 0;
			}
		}

		int[] marker1num = new int[5];
		int[] marker2num = new int[5];

		marker1num[0] = 0;
		marker1num[genotypes1.getMajorAllele()] = 1;
		marker1num[genotypes1.getMinorAllele()] = 2;
		marker2num[0] = 0;
		marker2num[genotypes2.getMajorAllele()] = 1;
		marker2num[genotypes2.getMinorAllele()] = 2;

		byte a1, a2, b1, b2;

		// iterate through all chromosomes in dataset
		for (int f = 0; f < keptFounders.size(); f++) {
			// at the moment, "i" will equal "f", since we only support founders
			Individual founder = keptFounders.get(f);
			//int i = founderIndices.get(f);

			// assign alleles for each of a pair of chromosomes at a marker
			// to four variables
			boolean haploid = snpInfo.isChrX() && founder.getGender().equals("1");
			if (haploid) {

				// haploid (x chrom/male)
				a1 = genotypes1.getLdFormatGenotypes()[f][0];
				a2 = genotypes2.getLdFormatGenotypes()[f][0];

				if (a1 != 0 && a2 != 0) {
					twoMarkerHaplos[marker1num[a1]][marker2num[a2]]++;
				}
			} else {
				// diploid
				a1 = genotypes1.getLdFormatGenotypes()[f][0];
				a2 = genotypes2.getLdFormatGenotypes()[f][0];
				b1 = genotypes1.getLdFormatGenotypes()[f][1];
				b2 = genotypes2.getLdFormatGenotypes()[f][1];

				if (a1 == 0 || a2 == 0 || b1 == 0 || b2 == 0) {
					// skip missing data
				} else if (((a1 >= 5 || b1 >= 5) && (a2 >= 5 || b2 >= 5))
						|| (a1 >= 5 && a2 != b2) || (a2 >= 5 && a1 != b1)) {
					doublehet++;
					// find doublehets and resolved haplotypes
				} else if (a1 >= 5 || b1 >= 5) {
					twoMarkerHaplos[1][marker2num[a2]]++;
					twoMarkerHaplos[2][marker2num[a2]]++;
				} else if (a2 >= 5 || b2 >= 5) {
					twoMarkerHaplos[marker1num[a1]][1]++;
					twoMarkerHaplos[marker1num[a1]][2]++;
				} else {
					twoMarkerHaplos[marker1num[a1]][marker2num[a2]]++;
					twoMarkerHaplos[marker1num[b1]][marker2num[b2]]++;
				}
			}
		}

		// another monomorphic marker check
		int r1 = twoMarkerHaplos[1][1] + twoMarkerHaplos[1][2];
		int r2 = twoMarkerHaplos[2][1] + twoMarkerHaplos[2][2];
		int c1 = twoMarkerHaplos[1][1] + twoMarkerHaplos[2][1];
		int c2 = twoMarkerHaplos[1][2] + twoMarkerHaplos[2][2];

		if ((r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0) && doublehet == 0) {
			return null;
		}

		double pA1, pB1, pA2, pB2;
		int total_chroms;

		known[AA] = twoMarkerHaplos[1][1];
		known[AB] = twoMarkerHaplos[1][2];
		known[BA] = twoMarkerHaplos[2][1];
		known[BB] = twoMarkerHaplos[2][2];
		unknownDH = doublehet;

		total_chroms = (int) (known[AA] + known[AB] + known[BA] + known[BB] + 2 * unknownDH);

		pA1 = (known[AA] + known[AB] + unknownDH) / (double) total_chroms;
		pB1 = 1 - pA1;
		pA2 = (known[AA] + known[BA] + unknownDH) / (double) total_chroms;
		pB2 = 1 - pA2;
		const_prob = 0.1;

		probHaps[AA] = const_prob;
		probHaps[AB] = const_prob;
		probHaps[BA] = const_prob;
		probHaps[BB] = const_prob;

		count_haps(0);
		estimate_p();

		// now we have an initial reasonable guess at p we can
		// start the EM - let the fun begin
		const_prob = 0.0;
		count = 1;
		loglike = -999999999;

		while (count < 1000) {
			oldloglike = loglike;
			count_haps(count);
			loglike = (known[AA] * Math.log(probHaps[AA]) + known[AB]
					* Math.log(probHaps[AB]) + known[BA]
					* Math.log(probHaps[BA]) + known[BB]
					* Math.log(probHaps[BB]))
					/ LN10
					+ ((double) unknownDH * Math.log(probHaps[AA]
							* probHaps[BB] + probHaps[AB] * probHaps[BA]))
					/ LN10;
			if (Math.abs(loglike - oldloglike) < TOLERANCE) {
				break;
			}
			estimate_p();
			count++;
		}

		num = probHaps[AA] * probHaps[BB] - probHaps[AB] * probHaps[BA];

		if (num < 0) {

			// flip matrix so we get the positive D'
			// flip AA with AB and BA with BB
			tmp = probHaps[AA];
			probHaps[AA] = probHaps[AB];
			probHaps[AB] = tmp;

			tmp = probHaps[BB];
			probHaps[BB] = probHaps[BA];
			probHaps[BA] = tmp;

			// flip frequency of second allele
			// done in this slightly asinine way because of a compiler
			// bugz0r in the dec-alpha version of java
			// which causes it to try to parallelize the swapping operations
			// and mis-schedules them
			pA2 = pA2 + pB2;
			pB2 = pA2 - pB2;
			pA2 = pA2 - pB2;

			// flip counts in the same fashion as p's
			tmp = numHaps[AA];
			numHaps[AA] = numHaps[AB];
			numHaps[AB] = tmp;

			tmp = numHaps[BB];
			numHaps[BB] = numHaps[BA];
			numHaps[BA] = tmp;

			// num has now undergone a sign change
			num = probHaps[AA] * probHaps[BB] - probHaps[AB] * probHaps[BA];

			// flip known array for likelihood computation
			tmp = known[AA];
			known[AA] = known[AB];
			known[AB] = tmp;

			tmp = known[BB];
			known[BB] = known[BA];
			known[BA] = tmp;
		}

		denom1 = (probHaps[AA] + probHaps[BA]) * (probHaps[BA] + probHaps[BB]);
		denom2 = (probHaps[AA] + probHaps[AB]) * (probHaps[AB] + probHaps[BB]);

		if (denom1 < denom2) {
			denom = denom1;
		} else {
			denom = denom2;
		}
		dprime = num / denom;

		// add computation of r^2 = (D^2)/p(1-p)q(1-q)
		rsq = num * num / (pA1 * pB1 * pA2 * pB2);

		return new LdResult(rsq, dprime);
	}*/

	/**
	 * Method ported from Haploview.
	 * 
	 * @param em_round
	 */
	private void count_haps(int em_round) {
		/*
		 * only the double heterozygote [AB][AB] results in ambiguous
		 * reconstruction, so we'll count the obligates then tack on the
		 * [AB][AB] for clarity
		 */
		numHaps[AA] = known[AA];
		numHaps[AB] = known[AB];
		numHaps[BA] = known[BA];
		numHaps[BB] = known[BB];

		if (em_round > 0) {
			numHaps[AA] += unknownDH
					* (probHaps[AA] * probHaps[BB])
					/ ((probHaps[AA] * probHaps[BB]) + (probHaps[AB] * probHaps[BA]));
			numHaps[BB] += unknownDH
					* (probHaps[AA] * probHaps[BB])
					/ ((probHaps[AA] * probHaps[BB]) + (probHaps[AB] * probHaps[BA]));
			numHaps[AB] += unknownDH
					* (probHaps[AB] * probHaps[BA])
					/ ((probHaps[AA] * probHaps[BB]) + (probHaps[AB] * probHaps[BA]));
			numHaps[BA] += unknownDH
					* (probHaps[AB] * probHaps[BA])
					/ ((probHaps[AA] * probHaps[BB]) + (probHaps[AB] * probHaps[BA]));
		}
	}

	/**
	 * Method ported from Haploview.
	 */
	private void estimate_p() {
		double total = numHaps[AA] + numHaps[AB] + numHaps[BA] + numHaps[BB]
				+ (4.0 * const_prob);

		probHaps[AA] = (numHaps[AA] + const_prob) / total;
		if (probHaps[AA] < 1e-10) {
			probHaps[AA] = 1e-10;
		}
		probHaps[AB] = (numHaps[AB] + const_prob) / total;
		if (probHaps[AB] < 1e-10) {
			probHaps[AB] = 1e-10;
		}
		probHaps[BA] = (numHaps[BA] + const_prob) / total;
		if (probHaps[BA] < 1e-10) {
			probHaps[BA] = 1e-10;
		}
		probHaps[BB] = (numHaps[BB] + const_prob) / total;
		if (probHaps[BB] < 1e-10) {
			probHaps[BB] = 1e-10;
		}
	}

	/**
	 * Method ported from Haploview. Used for calculating HWE.
	 */
	private double getHweValue(int[] parentHom, int parentHet)
			throws PriorityPrunerException {
		// ie: 11 13 31 33 -> homA =1 homB = 1 parentHet=2
		int homA = 0;
		int homB = 0;
		double pvalue = 0;
		for (int i = 0; i < parentHom.length; i++) {
			if (parentHom[i] != 0) {
				if (homA > 0) {
					homB = parentHom[i];
				} else {
					homA = parentHom[i];
				}
			}
		}

		// calculate p value from homA, parentHet and homB
		if (homA + parentHet + homB <= 0) {
			pvalue = 0;
		} else {
			pvalue = hwCalculate(homA, parentHet, homB);
		}
		return pvalue;
	}

	/**
	 * Ported from Haploview. Calculates exact two-sided Hardy-Weinberg p-value.
	 * Parameters are number of genotypes, number of rare alleles observed and
	 * number of heterozygotes observed.
	 * 
	 * (c) 2003 Jan Wigginton, Goncalo Abecasis
	 */
	private double hwCalculate(int obsAA, int obsAB, int obsBB)
			throws PriorityPrunerException {

		int diplotypes = obsAA + obsAB + obsBB;
		int rare = (obsAA * 2) + obsAB;
		int hets = obsAB;

		// make sure "rare" allele is really the rare allele
		if (rare > diplotypes) {
			rare = 2 * diplotypes - rare;
		}

		// make sure numbers aren't screwy
		if (hets > rare) {
			throw new PriorityPrunerException("HW test: " + hets
					+ "heterozygotes but only " + rare + "rare alleles.");
		}
		double[] tailProbs = new double[rare + 1];
		for (int z = 0; z < tailProbs.length; z++) {
			tailProbs[z] = 0;
		}

		// start at midpoint
		// all the casting is to make sure we don't overflow ints if there are
		// 10's of 1000's of inds
		int mid = (int) ((double) rare * (double) (2 * diplotypes - rare) / (double) (2 * diplotypes));

		// check to ensure that midpoint and rare alleles have same parity
		if (((rare & 1) ^ (mid & 1)) != 0) {
			mid++;
		}

		int het = mid;
		int hom_r = (rare - mid) / 2;
		int hom_c = diplotypes - het - hom_r;

		// calculate probability for each possible observed heterozygote
		// count up to a scaling constant, to avoid underflow and overflow
		tailProbs[mid] = 1.0;
		double sum = tailProbs[mid];

		for (het = mid; het > 1; het -= 2) {
			tailProbs[het - 2] = (tailProbs[het] * het * (het - 1.0))
					/ (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
			sum += tailProbs[het - 2];
			// 2 fewer hets for next iteration -> add one rare and one common
			// homozygote
			hom_r++;
			hom_c++;
		}

		het = mid;
		hom_r = (rare - mid) / 2;
		hom_c = diplotypes - het - hom_r;

		for (het = mid; het <= rare - 2; het += 2) {
			tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c)
					/ ((het + 2.0) * (het + 1.0));
			sum += tailProbs[het + 2];
			// 2 more hets for next iteration -> subtract one rare and one
			// common homozygote
			hom_r--;
			hom_c--;
		}

		for (int z = 0; z < tailProbs.length; z++) {
			tailProbs[z] /= sum;
		}

		double top = tailProbs[hets];
		for (int i = hets + 1; i <= rare; i++) {
			top += tailProbs[i];
		}

		double otherSide = tailProbs[hets];
		for (int i = hets - 1; i >= 0; i--) {
			otherSide += tailProbs[i];
		}

		if (top > 0.5 && otherSide > 0.5) {
			return 1.0;
		} else {
			if (top < otherSide) {
				return top * 2;
			} else {
				return otherSide * 2;
			}
		}
	}

	// public getters and setters for private fields of this class

	public void addResult(Result result) {
		results.add(result);
	}

	public SnpInfo getSnpInfo() {
		return snpInfo;
	}

	public void setSnpInfo(SnpInfo snpInfo) {
		this.snpInfo = snpInfo;
	}

	public ArrayList<Result> getResults() {
		return results;
	}

	public void setResults(ArrayList<Result> results) {
		this.results = results;
	}

	public boolean getIndexSnpPassed() {
		return indexSnpPassed;
	}
}
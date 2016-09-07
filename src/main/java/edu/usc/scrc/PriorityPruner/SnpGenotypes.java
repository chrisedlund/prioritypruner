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
 * This class stores all available information for a certain SNP. The genotypes
 * of this SNP is initially stored in a compressed format that later gets
 * converted by the SnpWorkUnit, to the format used for LD calculation. In order
 * to save memory, the decompressed representations get discarded as soon as
 * they've been used.
 */
public class SnpGenotypes {

	private String snpName;
	private SnpInfo snpInfo;
	private String allele1;
	private String allele2;
	private byte[] genotypes;
	private byte[][] ldFormatGenotypes;
	private double maf;
	private double a1Freq;


	private double missingPercent;
	private Byte minorAllele;
	private Byte majorAllele;
	//private int numMendelianErrors;
	private double hwePvalue;
	private boolean calculationDone = false;
	
	private boolean valid = true;
	


	/**
	 * Constructor for SnpGenotypes. Converts the original list of genotypes to
	 * the compressed version mentioned above.
	 * 
	 * @param snpName
	 *            name of this SNP
	 * @param snpInfo
	 *            the corresponding SnpInfo-object
	 * @param allele1
	 *            name of the first allele
	 * @param allele2
	 *            name of the second allele
	 * @param genotypes
	 *            original list of genotypes represnted as Strings
	 * @throws PriorityPrunerException
	 *             if an invalid genotype is encountered during the compression
	 */
	public SnpGenotypes(String snpName, SnpInfo snpInfo, String allele1,
			String allele2, String[] genotypes) throws PriorityPrunerException {
		this.snpName = snpName;
		this.snpInfo = snpInfo;
		this.allele1 = allele1;
		this.allele2 = allele2;
		this.genotypes = compressGenotypes(genotypes);
	}

	/**
	 * Converts the original list of genotypes, as provided by the user in the
	 * tped file, to a compressed format. In this format one byte is used
	 * to represent 4 genotypes (8 alleles) in following way: 00 - both alleles
	 * are missing, 01 - alleles are homozygous for allele 1, 10 - alleles are
	 * homozygous for allele 2, 11 - alleles are heterozygous.
	 * 
	 * @param genotypes
	 *            original list of genotypes stored as an array or Strings
	 * @return compressed version of the genotypes, stored as an array of bytes
	 * @throws PriorityPrunerException
	 *             if an invalid genotype is encountered during the conversion
	 */
	private byte[] compressGenotypes(String[] genotypes)
			throws PriorityPrunerException {

		byte[] genotypesByte = new byte[(genotypes.length / 8) + 1];
		int byteIndex = 0;
		byte byteValue = (byte) 0;

		for (int i = 0; i < genotypes.length; i += 8) {
			for (int j = 0; j < 8; j += 2) {
				if ((i + j) == genotypes.length) { 
					// shifts genotypes of the last byte over
					if (j == 2) {
						byteValue = (byte) (byteValue << 6);
					} else if (j == 4) {
						byteValue = (byte) (byteValue << 4);
					} else if (j == 6) {
						byteValue = (byte) (byteValue << 2);
					}
					break;
				}
				String alleleA = genotypes[i + j];
				String alleleB = genotypes[i + j + 1];
				byteValue = (byte) (byteValue << 2);

				// either allele is missing, gets set to 00 as default by bit
				// shifting
				if (alleleA.equals("0") || alleleB.equals("0")){
					
					// if only one of the alleles is missing
					if (!alleleA.equals(alleleB)){
						throw new PriorityPrunerException("Invalid genotype: "
								+ alleleA + alleleB + " found for locus "
								+ this.getSnpName());
					}
				    // alleles are homozygous for allele 1
				}else if (alleleA.equals(allele1) && alleleB.equals(allele1)) {
					byteValue += (byte) (1);

					// alleles are homozygous for allele 2
				} else if (alleleA.equals(allele2) && alleleB.equals(allele2)) {
					byteValue += (byte) (2);

					// alleles are heterozygous
				} else if ((alleleA.equals(allele1) && alleleB.equals(allele2))
						|| (alleleA.equals(allele2) && alleleB.equals(allele1))) {
					byteValue += (byte) (3);
				} else {
					throw new PriorityPrunerException("Invalid genotype: "
							+ alleleA + alleleB + " found for locus "
							+ this.getSnpName());
				}
			}
			genotypesByte[byteIndex] = byteValue;
			byteValue = (byte) 0;
			byteIndex++;
		}
		return genotypesByte;
	}

	/**
	 * Sets the genotypes stored in LD-format to null. This is done after LD
	 * calculation is performed by the SnpWorkUnit.
	 */
	public void eraseLdFormatGenotypes() {
		this.ldFormatGenotypes = null;
	}

	
	// public getters and setters for private fields of this class

	public String getSnpName() {
		return snpName;
	}

	public void setSnpName(String snpName) {
		this.snpName = snpName;
	}

	public SnpInfo getSnpInfo() {
		return snpInfo;
	}

	public String getAllele1() {
		return allele1;
	}

	public String getAllele2() {
		return allele2;
	}

	public byte[] getGenotypes() {
		return genotypes;
	}

	public void setGenotypes(byte[] genotypes) {
		this.genotypes = genotypes;
	}

	public byte[][] getLdFormatGenotypes() {
		return ldFormatGenotypes;
	}

	public void setLdFormatGenotypes(byte[][] ldFormatGenotypes) {
		this.ldFormatGenotypes = ldFormatGenotypes;
	}

	public double getMaf() {
		return maf;
	}

	public void setMaf(double maf) {
		this.maf = maf;
	}

	public double getMissingPercent() {
		return missingPercent;
	}

	public void setMissingPercent(double missingPercent) {
		this.missingPercent = missingPercent;
	}

	public Byte getMinorAllele() {
		return minorAllele;
	}

	public void setMinorAllele(Byte minorAllele) {
		this.minorAllele = minorAllele;
	}

	public Byte getMajorAllele() {
		return majorAllele;
	}

	public void setMajorAllele(Byte majorAllele) {
		this.majorAllele = majorAllele;
	}

//	public int getNumMendelianErrors() {
//		return numMendelianErrors;
//	}
//
//	public void setNumMendelianErrors(int numMendelianErrors) {
//		this.numMendelianErrors = numMendelianErrors;
//	}

//	public double getHwePvalue() {
//		return hwePvalue;
//	}
//
//	public void setHwePvalue(double hwePvalue) {
//		this.hwePvalue = hwePvalue;
//	}

	public boolean getCalculationDone() {
		return calculationDone;
	}

	public void setCalculationDone() {
		this.calculationDone = true;
	}
	
	public byte getByteGenotype(int sampleIndex){
		int compressedIndex = sampleIndex / 4;
		int offset = (3 - (sampleIndex % 4)) * 2;
		return (byte) ( (genotypes[compressedIndex] >> offset) & 3);
	}
	
	public byte getIntegerA1(int sampleIndex){
		int compressedIndex = sampleIndex / 4;
		int offset = (3 - (sampleIndex % 4)) * 2;
		byte genotype = (byte) ( (genotypes[compressedIndex] >> offset) & 3);
		
		if (genotype == 0){
			return 0;
		}else if (genotype == 1){
			return 1;
		}else if (genotype == 2){
			return 2;
		}else if (genotype == 3){
			return 1;
		}else{
			return 0;
		}
	}
	
	public byte getIntegerA2(int sampleIndex){
		int compressedIndex = sampleIndex / 4;
		int offset = (3 - (sampleIndex % 4)) * 2;
		byte genotype = (byte) ( (genotypes[compressedIndex] >> offset) & 3);
		
		if (genotype == 0){
			return 0;
		}else if (genotype == 1){
			return 1;
		}else if (genotype == 2){
			return 2;
		}else if (genotype == 3){
			return 2;
		}else{
			return 0;
		}
	}
	
	
	/**
	 * Method ported from Haploview. Calculates maf, hwe and missing genotype
	 * percentage for the current SNP.
	 * 
	 * @throws PriorityPrunerException
	 *             if an invalid genotype is encountered in a SNP
	 */

	public void calculateMafHweMissingPercentCompressed(ArrayList<Individual> keptFounders)
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
			boolean haploid  = this.getSnpInfo().isChrX() && founder.getSex() == Individual.Sex.MALE;
			
			byte genotype = this.getByteGenotype(f);
			
			if (!haploid) {
				numChromosomes+= 2;
				if (genotype != 0) {
					
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
				numChromosomes+= 1;
				if (genotype != 0) {
					
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
			this.setMaf((double) numAllele2
					/ (double) (numAllele1 + numAllele2));
			this.setMinorAllele((byte)2);
			this.setMajorAllele((byte)1);
		} else {
			this.setMaf((double) numAllele1
					/ (double) (numAllele1 + numAllele2));
			this.setMinorAllele((byte)1);
			this.setMajorAllele((byte)2);
		}
		this.a1Freq = (double) numAllele1 / (double) (numAllele1 + numAllele2);
		this.setMissingPercent((double) numMissing / numChromosomes);
//		this.setHwePvalue(this.getHweValue(founderHomCount,
//				founderHetCount));
		this.setCalculationDone();
		//System.out.println(genotypes.getSnpName() + "\t" + genotypes.getMissingPercent() + "\t" + genotypes.getHwePvalue() + "\t" + genotypes.getMaf());
	}

	/**
	 * Method ported from Haploview. Used for calculating HWE.
	 */
//	private double getHweValue(int[] parentHom, int parentHet)
//			throws PriorityPrunerException {
//		// ie: 11 13 31 33 -> homA =1 homB = 1 parentHet=2
//		int homA = 0;
//		int homB = 0;
//		double pvalue = 0;
//		for (int i = 0; i < parentHom.length; i++) {
//			if (parentHom[i] != 0) {
//				if (homA > 0) {
//					homB = parentHom[i];
//				} else {
//					homA = parentHom[i];
//				}
//			}
//		}
//
//		// calculate p value from homA, parentHet and homB
//		if (homA + parentHet + homB <= 0) {
//			pvalue = 0;
//		} else {
//			pvalue = hwCalculate(homA, parentHet, homB);
//		}
//		return pvalue;
//	}

	/**
	 * Ported from Haploview. Calculates exact two-sided Hardy-Weinberg p-value.
	 * Parameters are number of genotypes, number of rare alleles observed and
	 * number of heterozygotes observed.
	 * 
	 * (c) 2003 Jan Wigginton, Goncalo Abecasis
	 */
//	private double hwCalculate(int obsAA, int obsAB, int obsBB)
//			throws PriorityPrunerException {
//
//		int diplotypes = obsAA + obsAB + obsBB;
//		int rare = (obsAA * 2) + obsAB;
//		int hets = obsAB;
//
//		// make sure "rare" allele is really the rare allele
//		if (rare > diplotypes) {
//			rare = 2 * diplotypes - rare;
//		}
//
//		// make sure numbers aren't screwy
//		if (hets > rare) {
//			throw new PriorityPrunerException("HW test: " + hets
//					+ "heterozygotes but only " + rare + "rare alleles.");
//		}
//		double[] tailProbs = new double[rare + 1];
//		for (int z = 0; z < tailProbs.length; z++) {
//			tailProbs[z] = 0;
//		}
//
//		// start at midpoint
//		// all the casting is to make sure we don't overflow ints if there are
//		// 10's of 1000's of inds
//		int mid = (int) ((double) rare * (double) (2 * diplotypes - rare) / (double) (2 * diplotypes));
//
//		// check to ensure that midpoint and rare alleles have same parity
//		if (((rare & 1) ^ (mid & 1)) != 0) {
//			mid++;
//		}
//
//		int het = mid;
//		int hom_r = (rare - mid) / 2;
//		int hom_c = diplotypes - het - hom_r;
//
//		// calculate probability for each possible observed heterozygote
//		// count up to a scaling constant, to avoid underflow and overflow
//		tailProbs[mid] = 1.0;
//		double sum = tailProbs[mid];
//
//		for (het = mid; het > 1; het -= 2) {
//			tailProbs[het - 2] = (tailProbs[het] * het * (het - 1.0))
//					/ (4.0 * (hom_r + 1.0) * (hom_c + 1.0));
//			sum += tailProbs[het - 2];
//			// 2 fewer hets for next iteration -> add one rare and one common
//			// homozygote
//			hom_r++;
//			hom_c++;
//		}
//
//		het = mid;
//		hom_r = (rare - mid) / 2;
//		hom_c = diplotypes - het - hom_r;
//
//		for (het = mid; het <= rare - 2; het += 2) {
//			tailProbs[het + 2] = (tailProbs[het] * 4.0 * hom_r * hom_c)
//					/ ((het + 2.0) * (het + 1.0));
//			sum += tailProbs[het + 2];
//			// 2 more hets for next iteration -> subtract one rare and one
//			// common homozygote
//			hom_r--;
//			hom_c--;
//		}
//
//		for (int z = 0; z < tailProbs.length; z++) {
//			tailProbs[z] /= sum;
//		}
//
//		double top = tailProbs[hets];
//		for (int i = hets + 1; i <= rare; i++) {
//			top += tailProbs[i];
//		}
//
//		double otherSide = tailProbs[hets];
//		for (int i = hets - 1; i >= 0; i--) {
//			otherSide += tailProbs[i];
//		}
//
//		if (top > 0.5 && otherSide > 0.5) {
//			return 1.0;
//		} else {
//			if (top < otherSide) {
//				return top * 2;
//			} else {
//				return otherSide * 2;
//			}
//		}
//	}
	
	/**
	 * Checks if a SNP passes all of the user-defined filters.
	 * 
	 */
//	public void checkSnpValid(double minMaf, double minimumGenotypePercentage, double minimumHardyWeinbergPvalue) {
//
//		if ((this.getMaf() >= minMaf)
//				&& (((1 - this.getMissingPercent())) >= minimumGenotypePercentage)
//				&& (this.getHwePvalue() >= minimumHardyWeinbergPvalue)) {
//			this.valid = true;
//		} else {
//			this.valid = false;
//		}
//	}
	

	public boolean isValid() {
		return valid;
	}
	
	public void setValid(boolean valid) {
		this.valid = valid;
	}
	
	public double getA1Freq() {
		return a1Freq;
	}
}
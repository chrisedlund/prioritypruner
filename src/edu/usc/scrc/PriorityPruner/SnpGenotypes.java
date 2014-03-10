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
	private double missingPercent;
	private Byte minorAllele;
	private Byte majorAllele;
	private int numMendelianErrors;
	private double hwePvalue;
	private boolean calculationDone = false;

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

				// alleles are homozygous for allele 1
				if (alleleA.equals(allele1) && alleleB.equals(allele1)) {
					byteValue += (byte) (1);

					// alleles are homozygous for allele 2
				} else if (alleleA.equals(allele2) && alleleB.equals(allele2)) {
					byteValue += (byte) (2);

					// alleles are heterozygous
				} else if ((alleleA.equals(allele1) && alleleB.equals(allele2))
						|| (alleleA.equals(allele2) && alleleB.equals(allele1))) {
					byteValue += (byte) (3);

					// alleles are missing, gets set to 00 as default by bit
					// shifting
				} else if (alleleA.equals("0") && alleleB.equals("0")) {
				} else {
					throw new PriorityPrunerException("Invalid genotype: "
							+ alleleA + alleleB + " found in SNP \""
							+ this.getSnpName()
							+ "\".\nPlease redefine this genotype.");
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

	public int getNumMendelianErrors() {
		return numMendelianErrors;
	}

	public void setNumMendelianErrors(int numMendelianErrors) {
		this.numMendelianErrors = numMendelianErrors;
	}

	public double getHwePvalue() {
		return hwePvalue;
	}

	public void setHwePvalue(double hwePvalue) {
		this.hwePvalue = hwePvalue;
	}

	public boolean getCalculationDone() {
		return calculationDone;
	}

	public void setCalculationDone() {
		this.calculationDone = true;
	}
}
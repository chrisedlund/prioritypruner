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
 * This class represents a result record between a pair of SNPs. The data from
 * this object is displayed in the LD output file. The class implements
 * Comparable interface so that objects of this type can be sorted by the
 * position of the second SNP (ParterPos).
 */
public class Result implements Comparable<Object> {

	// variables for SNP 1 (index SNP)
	private String snpName;
	private double snpMaf;
	private String snpChr;
	private int snpPos;
	private SnpInfo indexSnp;

	// variables for SNP 2
	private String partnerSnpName;
	private double partnerMaf;
	private String partnerChr;
	private int partnerPos;
	private SnpInfo partnerSnp;

	// r^2 calculated between between SNP 1 and SNP 2
	private double rSquared;
	// D' calculated between between SNP 1 and SNP 2
	private double dPrime;

	/**
	 * Constructor for Result.
	 * 
	 * @param snpName
	 *            SNP 1 name
	 * @param snpMaf
	 *            SNP 1 MAF
	 * @param snpChr
	 *            SNP 1 Chr
	 * @param snpPos
	 *            SNP 1 Pos
	 * @param partnerSnpName
	 *            SNP 2 Name
	 * @param partnerMaf
	 *            SNP 2 MAF
	 * @param partnerChr
	 *            SNP 2 Chr
	 * @param partnerPos
	 *            SNP 2 Pos
	 * @param rSquared
	 *            r^2 calculated between SNP 1 and SNP 2
	 * @param dPrime
	 *            D' calculated between SNP 1 and SNP 2
	 * @param indexSnp
	 *            SnpInfo for SNP 1
	 * @param partnerSnp
	 *            SnpInfo for SNP 2
	 */
	public Result(String snpName, double snpMaf, String snpChr, int snpPos,
			String partnerSnpName, double partnerMaf, String partnerChr,
			int partnerPos, double rSquared, double dPrime, SnpInfo indexSnp,
			SnpInfo partnerSnp) {
		this.snpName = snpName;
		this.snpMaf = snpMaf;
		this.snpChr = snpChr;
		this.snpPos = snpPos;
		this.partnerSnpName = partnerSnpName;
		this.partnerMaf = partnerMaf;
		this.partnerChr = partnerChr;
		this.partnerPos = partnerPos;
		this.rSquared = rSquared;
		this.dPrime = dPrime;
		this.indexSnp = indexSnp;
		this.partnerSnp = partnerSnp;
	}

	/**
	 * Allows sorting of Result-objects, by ascending SNP 2 Pos.
	 */
	@Override
	public int compareTo(Object obj) {
		return Integer.valueOf(partnerPos).compareTo(
				Integer.valueOf(((Result) obj).getPartnerPos()));
	}

	// public getters and setters for private fields of this class

	public String getSnpName() {
		return snpName;
	}

	public double getSnpMaf() {
		return snpMaf;
	}

	public String getSnpChr() {
		return snpChr;
	}

	public int getSnpPos() {
		return snpPos;
	}

	public String getPartnerSnpName() {
		return partnerSnpName;
	}

	public double getPartnerMaf() {
		return partnerMaf;
	}

	public String getPartnerChr() {
		return partnerChr;
	}

	public int getPartnerPos() {
		return partnerPos;
	}

	public double getRSquared() {
		return rSquared;
	}

	public void setRSquared(double value) {
		this.rSquared = value;
	}

	public double getDPrime() {
		return dPrime;
	}

	public SnpInfo getIndexSnp() {
		return indexSnp;
	}

	public SnpInfo getPartnerSnp() {
		return partnerSnp;
	}
}
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
 * This class stores information about a certain SNP. When information from the
 * tped and SNP input file is matched, an object of this class will be stored in
 * a SnpGenotypes-object, which stores all information about a single SNP.
 */
public class SnpInfo implements Comparable<SnpInfo> {

	private String snpName;
	private String chr;
	private int pos;
	private String allele1;
	private String allele2;
	private double pValue;
	// design score for this SNP provided in the SNP input file
	private double designScore=Double.MAX_VALUE;
	private boolean forceInclude = false;
	//private int numBeadTypes;
	private boolean tagged = false;
	private boolean picked = false;
	private int pickOrder = -1;
	private boolean chrX = false;
	private double score;
	private int sortedByPosIndex;
	private boolean inTped = false;
	private SnpGenotypes snpGenotypes;
	private ArrayList<SnpR2Pair> taggedByList = new ArrayList<SnpR2Pair>();
	// metric weights defined for this SNP in SNP input file
	private double[] metrics;
	//private CommandLineOptions options = CommandLineOptions.getInstance();

	/**
	 * Constructor for SnpInfo. Makes quick checks to determine if this SNP has
	 * passed the design score minimum and if it's pseudo autosomal or not.
	 * 
	 * @param snpName
	 *            name of this SNP
	 * @param chr
	 *            chromosome of this SNP
	 * @param pos
	 *            base pair position
	 * @param allele1
	 *            first allele
	 * @param allele2
	 *            second allele
	 * @param pValue
	 *            associated p-value
	 * @param numAssays
	 *            associated number of bead types as defined in SNP input file
	 * @param forceInclude
	 *            flag showing whether or not to force include this SNP
	 * @param designScore
	 *            design score
	 * @param metrics
	 *            metric weights defined for this SNP in SNP input file
	 * @throws PriorityPrunerException
	 *             if unable to retrieve data from CommandLineOptions
	 */
	public SnpInfo(String snpName, String chr, int pos, String allele1,
			String allele2, double pValue,
			boolean forceInclude, double designScore, double[] metrics)
			throws PriorityPrunerException {
		this.snpName = snpName;
		this.chr = chr;
		this.pos = pos;
		this.metrics = metrics;
		this.allele1 = allele1;
		this.allele2 = allele2;
		this.pValue = pValue;
		//this.numBeadTypes = numAssays;
		this.forceInclude = forceInclude;
		this.designScore = designScore;

		if (chr.toUpperCase().equals("X")
				|| chr.toUpperCase().equals("CHRX")
				|| chr.toUpperCase().equals("23")) {
			this.chrX = true;
		}
		
	}

	/**
	 * Allows sorting of SnpInfo-objects, by ascending p-values. Implements
	 * comparison by SNP name in case p-value is the same.
	 */
	@Override
	public int compareTo(SnpInfo obj) {
		if ((Double.compare(this.pValue, obj.pValue)) == 0) {
			return (this.getSnpName().compareTo(obj.getSnpName()));
		} else {
			return (Double.compare(this.pValue, obj.pValue));
		}
	}

	/**
	 * Allows the storing of a SNP that has tagged this SNP, and the associated
	 * r^2 value. In cases where the SNP is tagged by a surrogate and the r^2
	 * value is unknown, "-1" is entered (in Pruner).
	 * 
	 * @param snpInfo
	 *            the SNP tagging this SNP
	 * @param rSquared
	 *            associated r^2 value
	 */
	public void addTaggedBy(SnpInfo snpInfo, double rSquared) {
		this.taggedByList.add(new SnpR2Pair(snpInfo, rSquared));
	}

	// public getters and setters for private fields of this class

	public String getSnpName() {
		return snpName;
	}

	public String getChr() {
		return chr;
	}

	public int getPos() {
		return pos;
	}

	public String getAllele1() {
		return allele1;
	}

	public void setAllele1(String value) {
		this.allele1 = value;
	}

	public String getAllele2() {
		return allele2;
	}

	public void setAllele2(String value) {
		this.allele2 = value;
	}

	public double getPValue() {
		return pValue;
	}

	public double getDesignScore() {
		return designScore;
	}
	
	public void setDesignScore(double value) {
		this.designScore=value;
	}
	
	public boolean getForceInclude() {
		return forceInclude;
	}

	public void setForceInclude(boolean value) {
		this.forceInclude = value;
	}

//	public int getNumBeadTypes() {
//		return numBeadTypes;
//	}
//	
//	public void setNumBeadTypes(int value) {
//		this.numBeadTypes=value;
//	}

	public boolean getTagged() {
		return tagged;
	}

	public void setTagged(boolean tagged) {
		this.tagged = tagged;
	}

	public boolean getPicked() {
		return picked;
	}

	public void setPicked(boolean picked) {
		this.picked = picked;
	}

	public int getPickOrder() {
		return pickOrder;
	}

	public void setPickOrder(int pickOrder) {
		this.pickOrder = pickOrder;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}

	public int getSortedByPosIndex() {
		return sortedByPosIndex;
	}

	public void setSortedByPosIndex(int sortedByPosIndex) {
		this.sortedByPosIndex = sortedByPosIndex;
	}

	public boolean getPassDesignScore() throws PriorityPrunerException {
		if (this.designScore >= CommandLineOptions.getInstance().getMinDesignScore()) {
			return true;
		}
		return false;
	}

	public SnpGenotypes getSnpGenotypes() {
		return snpGenotypes;
	}

	public void setSnpGenotypes(SnpGenotypes snpGenotypes) {
		this.snpGenotypes = snpGenotypes;
	}

	public ArrayList<SnpR2Pair> getTaggedByList() {
		return taggedByList;
	}

	public double[] getMetrics() {
		return metrics;
	}

	public boolean getInTped() {
		return inTped;
	}

	public void setInTped(boolean inTped) {
		this.inTped = inTped;
	}
	
	public boolean isChrX() {
		return chrX;
	}
	


}
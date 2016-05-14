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
 * This class is used to easily store a SNP together with an associated r^2
 * value. It's used by Pruner to store information about a SNP that is tagging
 * another and at what r^2 value this is done. The information is stored in a
 * "TaggedBy-list" in the SnpInfo-object corresponding to the SNP that is being
 * tagged.
 */
public class SnpR2Pair {

	private SnpInfo snp;
	private double rSquared;

	/**
	 * Constructor for SnpR2Pair.
	 * 
	 * @param snp
	 *            tagging SNP
	 * @param rSquared
	 *            associated r^2 value
	 */
	public SnpR2Pair(SnpInfo snp, double rSquared) {
		this.snp = snp;
		this.rSquared = rSquared;
	}

	// public getters and setters for private fields of this class

	public SnpInfo getSnp() {
		return snp;
	}

	public void setSnp(SnpInfo snp) {
		this.snp = snp;
	}

	public double getrSquared() {
		return rSquared;
	}

	public void setrSquared(double rSquared) {
		this.rSquared = rSquared;
	}
}
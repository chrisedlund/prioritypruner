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
 * This class stores the LD information for each result calculated between two
 * SNPs.
 */
public class LdResult {

	private double rSquared;
	private double dPrime;

	/**
	 * Constructor for LdResult.
	 * 
	 * @param rSquared
	 *            r^2-value for this LdResult-object
	 * @param dPrime
	 *            D'-value for this LdResult-object
	 */
	public LdResult(double rSquared, double dPrime) {
		this.rSquared = rSquared;
		this.dPrime = dPrime;
	}

	// public getters and setter for private fields of this class

	public double getRSquared() {
		return rSquared;
	}

	public void setRSquared(double rSquared) {
		this.rSquared = rSquared;
	}

	public double getDPrime() {
		return dPrime;
	}

	public void setDPrime(double dPrime) {
		this.dPrime = dPrime;
	}
}
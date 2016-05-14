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
 * This class stores an r^2 threshold together with its associated p-value. The
 * class implements the Comparable interface so that objects of this type can be
 * sorted by their p-value.
 */
public class R2Threshold implements Comparable<R2Threshold> {

	private double pValue;
	private double r2Threshold;

	/**
	 * Constructor for R2Threshold.
	 * 
	 * @param pValue
	 *            p-value associated with current threshold
	 * @param r2Threshold
	 *            threshold
	 */
	public R2Threshold(double pValue, double r2Threshold) {
		this.pValue = pValue;
		this.r2Threshold = r2Threshold;
	}

	/**
	 * Allows sorting of R2Threshold-objects, by ascending p-values.
	 */
	@Override
	public int compareTo(R2Threshold o) {
		return Double.compare(pValue, o.getPValue());
	}

	// public getters and setters for private fields of this class

	public double getPValue() {
		return pValue;
	}

	public double getR2Threshold() {
		return r2Threshold;
	}
}
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
 * This class stores the desired number of surrogates together with an
 * associated p-value. The class implements the Comparable interface so that
 * objects of this type can be sorted by their p-value.
 */
public class SurrogateThreshold implements Comparable<SurrogateThreshold> {

	private double pValue;
	private int numSurrogates;

	/**
	 * Constructor for SurrogateThreshold.
	 * 
	 * @param pValue
	 *            p-value for the associated number of surrogates
	 * @param numSurrogates
	 *            number of surrogates
	 */
	public SurrogateThreshold(double pValue, int numSurrogates) {
		this.pValue = pValue;
		this.numSurrogates = numSurrogates;
	}

	/**
	 * Allows sorting of SurrogateThreshold-objects, by ascending p-values.
	 */
	@Override
	public int compareTo(SurrogateThreshold o) {
		return Double.compare(pValue, o.getPValue());
	}

	// public getters and setters for private fields of this class

	public double getPValue() {
		return pValue;
	}

	public int getNumSurrogates() {
		return numSurrogates;
	}
}
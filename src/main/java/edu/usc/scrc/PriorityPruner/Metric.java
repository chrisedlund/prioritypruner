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
 * This class stores the name and weight of a metric.
 */
public class Metric {

	// name of metric, has to have a matching counterpart in the SNP input file
	private String name;
	private double weight;

	/**
	 * Constructor for Metric.
	 * 
	 * @param name
	 *            name of metric
	 * @param weight
	 *            weight of metric
	 */
	public Metric(String name, double weight) {
		this.name = name;
		this.weight = weight;
	}

	
	// public getters and setters for private fields of this class

	public String getName() {
		return name;
	}

	public void setColumnName(String columnName) {
		this.name = columnName;
	}

	public double getWeight() {
		return weight;
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}
}
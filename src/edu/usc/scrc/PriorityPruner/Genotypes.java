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
 * This is the superclass of the different input formats that eventually will be
 * implemented. Right now TPlink is the only class extending it. This class
 * stores all parsed information about SNPs and individuals.
 */
public class Genotypes {
	// SNPs stored as SnpGenotypes-object
	protected ArrayList<SnpGenotypes> snpGenotypes = new ArrayList<SnpGenotypes>();
	// gender of chosen individuals
	protected ArrayList<String> subjectSexes = new ArrayList<String>();
	// all individuals, each object contains a flag indicating if it's chosen or
	// not
	protected ArrayList<Individual> individuals = new ArrayList<Individual>();
	// since we only support founders at the moment, this is basically just an
	// index over chosen individuals
	protected ArrayList<Integer> founderIndices = new ArrayList<Integer>();

	// public getters and setter for private fields of this class

	public ArrayList<SnpGenotypes> getSnpGenotypes() {
		return snpGenotypes;
	}

	public ArrayList<String> getSubjectSexes() {
		return subjectSexes;
	}

	public ArrayList<Individual> getIndividuals() {
		return individuals;
	}

	public ArrayList<Integer> getFounderIndices() {
		return founderIndices;
	}
}
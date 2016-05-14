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
 * This class stores the family information for a certain individual. It also
 * has a "keep-flag" indicating if the individual should be included in the
 * current run.
 */
public class Individual {
	
	public enum Sex {
	    MALE, FEMALE, UNKNOWN 
	}
	
	private String familyID;
	private String individualID;
	private String momID;
	private String dadID;
	//private String gender;
	private Individual mom;
	private Individual dad;
	private boolean keep;
	private Sex sex;

	/**
	 * Constructor for Individual.
	 * 
	 * @param familyID
	 *            family identification
	 * @param individualID
	 *            individual identification
	 * @param dadID
	 *            identification string for father
	 * @param momID
	 *            identification string for mother
	 * @param gender
	 *            gender of individual
	 * @throws PriorityPrunerException 
	 */
	public Individual(String familyID, String individualID, String dadID,
			String momID, String sexString) throws PriorityPrunerException {
		this.familyID = familyID;
		this.individualID = individualID;
		this.dadID = dadID;
		this.momID = momID;
		if (sexString.equals("1")){
			this.sex = Sex.MALE;
		}else if (sexString.equals("2")){
			this.sex = Sex.FEMALE;
		}else{
			throw new PriorityPrunerException("Invalid sex " + sexString + " for [ " + familyID + " " + individualID + " ]");
		}
		//this.gender = sexString;
	}

	// public getters and setter for private fields of this class


	public String getFamilyID() {
		return familyID;
	}

	public void setFamilyID(String familyID) {
		this.familyID = familyID;
	}

	public String getIndividualID() {
		return individualID;
	}

	public void setIndividualID(String individualID) {
		this.individualID = individualID;
	}

	public String getDadID() {
		return dadID;
	}

	public void setDadID(String dadID) {
		this.dadID = dadID;
	}

	public String getMomID() {
		return momID;
	}

	public void setMomID(String momID) {
		this.momID = momID;
	}

//	public String getGender() {
//		return gender;
//	}
//
//	public void setGender(String gender) {
//		this.gender = gender;
//	}

	public Individual getMom() {
		return mom;
	}

	public void setMom(Individual mom) {
		this.mom = mom;
	}

	public Individual getDad() {
		return dad;
	}

	public void setDad(Individual dad) {
		this.dad = dad;
	}

	public boolean getKeep() {
		return keep;
	}

	public void setKeep(boolean keep) {
		this.keep = keep;
	}

	public Sex getSex() {
		return sex;
	}
}
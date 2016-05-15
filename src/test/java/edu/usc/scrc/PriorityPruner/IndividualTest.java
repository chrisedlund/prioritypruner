package edu.usc.scrc.PriorityPruner;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class IndividualTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	/***
	 * Test that passing an string value of "1" for sex sets the Sex of Individual to Individual.Sex.Male 
	 * @throws PriorityPrunerException
	 */
	@Test
	public void testValidMaleSex() throws PriorityPrunerException{
		Individual ind = new Individual("", "", "", "", "1");
		assertEquals(ind.getSex(), Individual.Sex.MALE);
	}

	/***
	 * Test that passing an string value of "2" for sex sets the Sex of Individual to Individual.Sex.Female 
	 */
	@Test
	public void testValidFemaleSex() throws PriorityPrunerException {
		Individual ind = new Individual("", "", "", "", "2");
		assertEquals(ind.getSex(), Individual.Sex.FEMALE);
	}
	
	/***
	 * Test that passing an string value not equal to "1" or "2" for sex throws a PriorityPrunerException
	 */
	@Test(expected = PriorityPrunerException.class)
	public void testInvalidSex() throws Exception {
		Individual ind = new Individual("", "", "", "", "0");	
	}
	
	/***
	 * Test all getters for this class, that they equal what we pass them
	 * @throws PriorityPrunerException
	 */
	@Test
	public void testGetters() throws PriorityPrunerException{
		Individual ind = new Individual("famid","indid","dadid","momid","1");
		assertEquals(ind.getMomID(),"momid");
		assertEquals(ind.getDadID(),"dadid");
		assertEquals(ind.getFamilyID(),"famid");
		assertEquals(ind.getIndividualID(),"indid");
		assertEquals(ind.getSex(),Individual.Sex.MALE);
				
		// test keep is initialize to true
		assertEquals(ind.getKeep(),true);
	}

}

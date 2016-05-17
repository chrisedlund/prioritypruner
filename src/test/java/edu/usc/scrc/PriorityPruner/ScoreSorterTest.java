package edu.usc.scrc.PriorityPruner;

import static org.junit.Assert.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

/***
 * This test verifies that ScoreSorter is sorting Result objects 
 * by these members:
 *  1) partner SNP forceInclude (true before false)
 *  2) partner SNP score (higher before lower)
 *  3) partner SNP name (reverse lexicographic order) 
 * @author cedlund
 *
 */
public class ScoreSorterTest {

	ArrayList<Result> expectedResults = null;
	
	@Before
	public void setUp() throws Exception {
		
		// set up a small number of partner snps (we don't need index snp since it's not use
		// unneeded parameters are set to 0 or ""
		// the order of the SnpInfo are in the same order as the Result objects and are listed
		//  in the order they should end up in (p1-p8) in their parent result object
		// order of SnpInfo parameters are:
		//   snpName, chr, pos, allele1, allele2, pValue, forceInclude, designScore, metrics
		SnpInfo p1 = new SnpInfo("b", "", 0, "", "", 0, true, 0, null); p1.setScore(0.8);
		SnpInfo p2 = new SnpInfo("a", "", 0, "", "", 0, true, 0, null); p2.setScore(0.8);
	
		SnpInfo p3 = new SnpInfo("f", "", 0, "", "", 0, true, 0, null); p3.setScore(0.7);
		SnpInfo p4 = new SnpInfo("e", "", 0, "", "", 0, true, 0, null); p4.setScore(0);
		
		SnpInfo p5 = new SnpInfo("i", "", 0, "", "", 0, false, 0, null); p5.setScore(1.0);
		SnpInfo p6 = new SnpInfo("h", "", 0, "", "", 0, false, 0, null); p6.setScore(0.5);
		SnpInfo p7 = new SnpInfo("g", "", 0, "", "", 0, false, 0, null); p7.setScore(0.5);
		SnpInfo p8 = new SnpInfo("j", "", 0, "", "", 0, false, 0, null); p8.setScore(0);
		
		// set up result objects to hold the SnpInfo objects
		// order of params are: snpName, snpMaf, snpChr, snpPos, partnerSnpName, partnerMaf, 
		//                        partnerChr, partnerPos, rSquared, double dPrime, indexSnp, partnerSnp
		// unneeded parameters are set to 0 or ""
		
		Result r1 = new Result("", 0, "", 0, p1.getSnpName(), 0, "", 0, 0, 0, null, p1);
		Result r2 = new Result("", 0, "", 0, p2.getSnpName(), 0, "", 0, 0, 0, null, p2);
		Result r3 = new Result("", 0, "", 0, p3.getSnpName(), 0, "", 0, 0, 0, null, p3);
		Result r4 = new Result("", 0, "", 0, p4.getSnpName(), 0, "", 0, 0, 0, null, p4);
		Result r5 = new Result("", 0, "", 0, p5.getSnpName(), 0, "", 0, 0, 0, null, p5);
		Result r6 = new Result("", 0, "", 0, p6.getSnpName(), 0, "", 0, 0, 0, null, p6);
		Result r7 = new Result("", 0, "", 0, p7.getSnpName(), 0, "", 0, 0, 0, null, p7);
		Result r8 = new Result("", 0, "", 0, p8.getSnpName(), 0, "", 0, 0, 0, null, p8);
		
		expectedResults = new ArrayList<Result>(Arrays.asList(r1,r2,r3,r4,r5,r6,r7,r8));
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testSorter() {
		// copy expected results into new array and randomly shuffle, then sort, then compare to expected
		//  do this 100 times - the result should be the same each time
		ScoreSorter sorter = new ScoreSorter();
		for (int i =0; i<100; i++){
			// set up temp array
			ArrayList<Result> tempList = new ArrayList<Result>(expectedResults.size());
			for (Result result: expectedResults){
				tempList.add(result);
			}
			
			// shuffle
			Collections.shuffle(tempList);
			
			// sort with ScoreSorter
			Collections.sort(tempList, sorter);
			
			// verify order is correct
			for (int j=0; j<tempList.size(); j++){
				assertSame(tempList.get(j), expectedResults.get(j));
			}
			
		}
		
		
	}

}

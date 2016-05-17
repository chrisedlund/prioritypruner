package edu.usc.scrc.PriorityPruner;

import static org.junit.Assert.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.StringWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;

import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Level;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

import edu.usc.scrc.PriorityPruner.ExitDeniedSecurityManager.ExitSecurityException;
public class PriorityPrunerTest {

	SecurityManager originalSecurityManager;
	ArrayList<Level> loggerLevels = new ArrayList<Level>();
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	/***
	 * This is called before each test. Set up the ExitDeniedSecurityManager to intercept system.exit() calls
	 * @throws Exception
	 */
	@Before
	public void setUp() throws Exception {
		this.originalSecurityManager = System.getSecurityManager();
		System.setSecurityManager(new ExitDeniedSecurityManager());
		
		// save logger levels for restoring later
//		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
//		loggers.add(LogManager.getRootLogger());
//		for ( Logger logger : loggers ) {
//			loggerLevels.add(logger.getLevel());
//		}
		
	}

	@After
	public void tearDown() throws Exception {
		System.setSecurityManager(this.originalSecurityManager);
		
		// restore loggers
		// save logger levels for restoring later
//		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
//		loggers.add(LogManager.getRootLogger());
//		int i=0;
//		for ( Logger logger : loggers ) {
//			logger.setLevel(loggerLevels.get(i));
//			i++;
//		}
		
		
	}

	/***
	 * The TemporaryFolder Rule allows creation of files and folders that should be deleted when the test method finishes (whether it passes or fails).
	 */
	@Rule
	public TemporaryFolder tempFolder= new TemporaryFolder();
	  
	
	/***
	 * v0.1.2 Regression Test 1 
	 * This test verifies the output for a given set of complex command line options and input data, that the output
	 * data including main results and ld calculations are identical to v0.1.2
	 * @throws IOException 
	 * @throws URISyntaxException 
	 */
	@Test
	public void regressionTest_v0_1_2_test1() throws IOException, URISyntaxException{
		
		// get paths to test resource files
		ClassLoader classLoader = getClass().getClassLoader();
		String tpedPath = classLoader.getResource("regression_v0.1.2/pp_1kgp3_yri_chr22-X_regression_test.tped").getPath();
		String tfamPath = classLoader.getResource("regression_v0.1.2/pp_1kgp3_yri_chr22-X_regression_test.tfam").getPath();
		String snpTablePath = classLoader.getResource("regression_v0.1.2/pp_1kgp3_yri_chr22-X_regression_test.snp_table.txt").getPath(); 
		String keepPath = classLoader.getResource("regression_v0.1.2/test1_keep.txt").getPath();
		
		
		System.out.println("Creating temporary files in: " + tempFolder.getRoot());
		
		// get paths to temporary output files
		File resultsFile = tempFolder.newFile("test1.results");
		File ldFile = tempFolder.newFile("test1.ld");
		File logFile = tempFolder.newFile("test1.log");
		
		
		File outputPrefix = new File(tempFolder.getRoot(), "test1");
		
		// build command line options
		StringBuilder command = new StringBuilder("");
		command.append("--tped " + tpedPath);
		command.append(" --tfam " + tfamPath);
		command.append(" --snp_table " + snpTablePath);
		command.append(" --r2t 1 0.3 --r2t 0.01 0.5"); // r2 thresholds
		command.append(" --min_design_score 0.2");
		command.append(" --max_distance 30000");
		command.append(" --min_snp_callrate 0.9");
		command.append(" --out " + outputPrefix.getPath()); // output prefix
		command.append(" --ld"); // create ld output file
		command.append(" --st 1 0 --st 0.01 2");
		command.append(" --metric p 10");
		command.append(" --metric metric1 5");
		command.append(" --metric metric2 3");
		command.append(" --keep " + keepPath); // subset of individuals to keep
		
		// run PriorityPruner
		try{
			PriorityPruner.main(command.toString().split(" "));
			fail("Expected system.exit to be called from PriorityPruner");
		}catch(ExitSecurityException e) {
			// the system.exit call from PriorityPruner was intercepted; now make sure exit code was 0
			int status = e.getStatus();
            assertEquals(0, status);
		}

		// turn off all log4j loggers so we don't get errors
		
//		List<Logger> loggers = Collections.<Logger>list(LogManager.getCurrentLoggers());
//		loggers.add(LogManager.getRootLogger());
//		for ( Logger logger : loggers ) {
//		    logger.setLevel(Level.OFF);
//		}
//		
		
		// compare to expected results
		runPriorityPrunerAndVerify(command.toString(), "regression_v0.1.2/test1_expected", resultsFile, ldFile, logFile);
		
	}
	

	/***
	 * v0.1.2 Regression Test 2
	 * This test verifies the output for a given set of complex command line options and input data, that the output
	 * data including main results and ld calculations are identical to v0.1.2
	 * @throws IOException 
	 * @throws URISyntaxException 
	 */
	@Test
	public void regressionTest_v0_1_2_test2() throws IOException, URISyntaxException{

		// get paths to test resource files
		ClassLoader classLoader = getClass().getClassLoader();
		String tpedPath = classLoader.getResource("regression_v0.1.2/pp_1kgp3_yri_chr22-X_regression_test.tped").getPath();
		String tfamPath = classLoader.getResource("regression_v0.1.2/pp_1kgp3_yri_chr22-X_regression_test.tfam").getPath();
		String snpTablePath = classLoader.getResource("regression_v0.1.2/pp_1kgp3_yri_chr22-X_regression_test.snp_table.txt").getPath(); 
		String removePath = classLoader.getResource("regression_v0.1.2/test2_remove.txt").getPath();
		
		
		System.out.println("Creating temporary files in: " + tempFolder.getRoot());
		
		// get paths to temporary output files, that should be deleted automatically on test end
		File resultsFile = tempFolder.newFile("test2.results");
		File ldFile = tempFolder.newFile("test2.ld");
		File logFile = tempFolder.newFile("test2.log");
		
		File outputPrefix = new File(tempFolder.getRoot(), "test2");
		
		// build command line options
		StringBuilder command = new StringBuilder("");
		command.append("--tped " + tpedPath);
		command.append(" --tfam " + tfamPath);
		command.append(" --snp_table " + snpTablePath);
		command.append(" --r2 0.3"); // r2 thresholds
		command.append(" --out " + outputPrefix); // output prefix
		command.append(" --ld"); // create ld output file
		command.append(" --st 0.01 3");
		command.append(" --metric p 10");
		command.append(" --max_distance 30000");
		command.append(" --remove " + removePath); // subset of individuals to remove
		command.append(" --do_not_pick_force_included_first");
		command.append(" --no_surrogates_for_force_included_snps");
		
		// run PriorityPruner
		runPriorityPrunerAndVerify(command.toString(), "regression_v0.1.2/test2_expected", resultsFile, ldFile, logFile);
		
	}
	
	
	private void runPriorityPrunerAndVerify(String command, String expectedOutputPrefix, File resultsFile, File ldFile, File logFile) throws UnsupportedEncodingException, URISyntaxException, IOException{
		// run PriorityPruner
		try{
			PriorityPruner.main(command.split(" "));
			fail("Expected system.exit to be called from PriorityPruner");
		}catch(ExitSecurityException e) {
			// the system.exit call from PriorityPruner was intercepted; now make sure exit code was 0
			int status = e.getStatus();
            assertEquals(0, status);
		}
		
		// compare to expected results
		try{
			compareToExpectedResults(expectedOutputPrefix, resultsFile, ldFile, logFile);
		}catch (ExitSecurityException e) {}
	}
	
	
	/***
	 * Compares output from PriorityPruner with expected output in test resource folder
	 * @param expectedOutputPrefix
	 * @param resultsFile
	 * @param ldFile
	 * @param logFile
	 * @throws URISyntaxException
	 * @throws UnsupportedEncodingException
	 * @throws IOException
	 */
	private void compareToExpectedResults(String expectedOutputPrefix, File resultsFile, File ldFile, File logFile) throws URISyntaxException, UnsupportedEncodingException, IOException{
		
		ClassLoader classLoader = getClass().getClassLoader();
		
		// get path to expected results file
		java.net.URL expectedResultsUrl = classLoader.getResource(expectedOutputPrefix + ".results");
		java.nio.file.Path expectedResultsPath = java.nio.file.Paths.get(expectedResultsUrl.toURI());
		
		// get path to test results file
		java.nio.file.Path testResultsPath = java.nio.file.Paths.get(resultsFile.getPath());
				
		// check the .results file equals the expected
		String expectedResults = new String(java.nio.file.Files.readAllBytes(expectedResultsPath), "UTF8");
		String testResults = new String(java.nio.file.Files.readAllBytes(testResultsPath), "UTF8");
		
		//System.out.println("Expected results string is of length: " + expectedResults.length());
		//System.out.println("Test results string is of length: " + testResults.length());
		
		assertEquals(expectedResults, testResults);
		
		// get path to expected ld file
		java.net.URL expectedLdUrl = classLoader.getResource(expectedOutputPrefix + ".ld.gz");
		java.nio.file.Path expectedLdPath = java.nio.file.Paths.get(expectedResultsUrl.toURI());
		
		// get path to test ld file
		java.nio.file.Path testLdPath = java.nio.file.Paths.get(ldFile.getPath());
		
		// we need to gunzip the expected ld file since it's too large to store uncompressed in the code repository
		FileInputStream in = new FileInputStream(expectedLdUrl.getPath());
		GZIPInputStream  ungzippedStream = new GZIPInputStream(in);
		Reader reader = new InputStreamReader(ungzippedStream, "UTF-8");
		Writer writer = new StringWriter();
		char[] buffer = new char[10240];
		for (int length = 0; (length = reader.read(buffer)) > 0;) {
			writer.write(buffer, 0, length);
		}
		String expectedLd = writer.toString();
		writer.close();
		reader.close();
		ungzippedStream.close();
		in.close();
		
		// check the .ld file equals the expected
		//String expectedLd = new String(java.nio.file.Files.readAllBytes(expectedLdPath), "UTF8");
		String testLd = new String(java.nio.file.Files.readAllBytes(testLdPath), "UTF8");
		//System.out.println("Expected string is of length: " + expectedLd.length());
		//System.out.println("Test string is of length: " + testLd.length());
		assertEquals(expectedLd, testLd);
		
	}
	
	

	
}

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

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;


/**
 * This class is responsible for pruning and prioritizing the SNPs defined in
 * the input files. It initiates the parsing of these files, goes through
 * current SNPs, checks that valid properties are defined, initiates LD
 * calculation, and based on these LD results - picks and tags SNPs. It also
 * creates output files for the LD calculations and pruning results.
 */
public class Pruner {

	private SnpListFile snpListFile;
	private LinkageDisequilibriumFile ldFile;
	private Genotypes genotypes;
	private int pickOrder = 1;
	private CommandLineOptions options = CommandLineOptions.getInstance();
	//private BufferedWriter ldWriter = null;

	/**
	 * Constructor for Pruner. Initiates parsing of genotype and family data by
	 * creating SnpListFile- and TPlink-objects. It creates an output file for
	 * the LD results calculated by SnpWorkUnit and initiates the pruning.
	 * 
	 * @param snpListFile
	 *            The SnpListFile containing the SNPs to prune
	 * @param ldFile
	 *            The LinkageDisequilibriumFile to write to (if not null)
	 *            
	 * @throws PriorityPrunerException
	 *             if problems are encountered during parsing or pruning
	 */
	public Pruner(Genotypes genotypes, SnpListFile snpListFile, 
			LinkageDisequilibriumFile ldFile) throws PriorityPrunerException {
		
		this.snpListFile = snpListFile;
		this.ldFile = ldFile;
		
		// parse the list of samples defined by the --keep or --remove options if specified
//		PlinkSampleListFile keepRemoveSamples = null;
//		if (options.getKeep() != null){
//			keepRemoveSamples = new PlinkSampleListFile(options.getKeep());
//		}else if (options.getRemove() != null){
//			keepRemoveSamples = new PlinkSampleListFile(options.getRemove());
//		}
		
		this.genotypes = genotypes;
//		this.genotypes = new TPlink(options.getTped(), options.getTfam(),
//				this.snpListFile, keepRemoveSamples);
		
		try{
			//checkSnpsAreInGenotypeFile();
			
			//LogWriter.getLogger().info("Calculating SNP statistics");
			// calculate maf, call rate, hwe for all snps
			// already been done for current SNP - calculate this!
			int numFailMaf = 0;
			//int numFailHwe = 0;
			int numFailCallRate = 0;
			int numFail = 0;
			for (SnpGenotypes g : this.genotypes.getSnpGenotypes()) {
				//getMafHweMissingPercent(genotype, founderIndices,subjectSexes);
				g.calculateMafHweMissingPercentCompressed(genotypes.getKeptFounders());
				boolean fail = false;
				if (g.getMaf() < options.getMinMaf()){
					numFailMaf++;
					fail = true;
				}
//				if (g.getHwePvalue() < options.getMinHwe()){
//					numFailHwe++;
//					fail = true;
//				}
				if ( (1 - g.getMissingPercent()) <  options.getMinSnpCallRate()){
					numFailCallRate++;
					fail = true;
				}
				if (fail){
					g.setValid(false);
					numFail++;
				}
				
				//g.checkSnpValid(options.getMinMaf(),options.getMinSnpCallRate(), options.getMinHwe());
			}
			if (numFailMaf >0){
				LogWriter.getLogger().info(numFailMaf + " SNPs failed frequency test ( MAF < " + options.getMinMaf() + " )");
			}
//			if (numFailHwe >0){
//				LogWriter.getLogger().info(numFailMaf + " SNPs failed HWE test ( p < " + options.getMinHwe() + " )");
//			}
			if (numFailCallRate >0){
				LogWriter.getLogger().info(numFailCallRate + " SNPs failed callrate test ( callrate < " + options.getMinSnpCallRate() + " )");
			}
			if (numFail > 0){
				LogWriter.getLogger().info("After filtering, there are " + (this.genotypes.getSnpGenotypes().size() - numFail) + " SNPs");
			}
			
			
//			if (options.isOutputLDTable()){
//				createLdFile();
//			}
			startPruning();
		}catch (PriorityPrunerException e){
			throw e;
		}
//		finally{
//			if (ldWriter != null){
//				try{
//					ldWriter.close();
//				}catch(IOException e){e.printStackTrace();}
//			}
//		}
	}

//	/**
//	 * Creates LD file to store information calculated in SnpWorkUnit. The file
//	 * will be placed in file path specified by the user. It will overwrite
//	 * already existing file with same name.
//	 * 
//	 * @throws PriorityPrunerException
//	 *             if new file couldn't be created
//	 */
//	private void createLdFile() throws PriorityPrunerException {
//
//
//		
//		try {
//			// creates output file
//			LogWriter.getLogger().info("Writing LD metrics to [ " + options.getOutputPrefix() + ".ld" + " ]");
//			File outputFile = new File(options.getOutputPrefix() + ".ld");
//			ldWriter = new BufferedWriter(new FileWriter(outputFile));
//			
//			ldWriter.write("index_snp_name" + "\t" + "index_snp_chr" + "\t" + "index_snp_pos" + "\t"
//			+ "index_snp_a1" + "\t" + "index_snp_a2" + "\t" + "partner_snp_name" + "\t" + 
//					"partner_snp_chr" + "\t" + "partner_snp_pos" + "\t" + "partner_snp_a1" + "\t"
//			+ "partner_snp_a2" + "\t" + "r^2" + "\t" + "D'" + "\n");
//			
//		} catch (IOException e) {
//			throw new PriorityPrunerException("Could not create file: "
//					+ e.getMessage());
//		}
//		
//	}
	
	
//	/**
//	 * Checks that all SNPs defined in the SNP Input Table have corresponding data in the
//	 * genotype dataset.
//	 * 
//	 * @throws PriorityPrunerException
//	 *             if SNP not found in the genotype dataset
//	 */
//	private void checkSnpsAreInGenotypeFile() throws PriorityPrunerException{
//		// checks that corresponding SNP is provided in tped file
//		for (SnpInfo snp : snpListFile.getSnps()) {
//			if (!snp.getInTped()) {
//				throw new PriorityPrunerException(snp.getSnpName()
//						+ " at chromsome " + snp.getChr()
//						+ ":" + snp.getPos()
//						+ " with alleles: " + snp.getAllele1() + " "
//						+ snp.getAllele2() + ", not found in genotype dataset.");
//			}
//		}
//	}

	/**
	 * Initiates pruning of SNPs defined in the SNP input file. Loops through
	 * provided SNPs and checks that chromosome is valid and that information
	 * about the SNP was also in tped file. Prunes the current SNP if it's not
	 * already picked and either force included or not tagged.
	 * 
	 * @throws PriorityPrunerException
	 *             if invalid information is encountered
	 * @throws IOException 
	 */
	private void startPruning() throws PriorityPrunerException {

		//int prunedSnpIndex = 0;

		// loops through all SNPs in the SNP Input Table in order of ascending p-value
		for (SnpInfo snp : snpListFile.getSnps()) {

			// prunes current SNP if it's not already picked and either force
			// included or not tagged
			if (!snp.getPicked() && (!snp.getTagged() || snp.getForceInclude())) {
				LogWriter.getLogger().debug(
					"\nUsing index SNP " + snp.getSnpName() + " with p-value: " + snp.getPValue());
				prune(snp);
				//prunedSnpIndex++;

				// prints out how many SNPs that currently are processed, if
				// verbose option is entered
//				if (prunedSnpIndex % 1000 == 0) {
//					LogWriter.getLogger().debug(
//							"\nProcessed " + prunedSnpIndex + " SNPs\n");
//				}
				
			} else {
				LogWriter.getLogger().debug(
						"\nSkipping " + snp.getSnpName()
								+ " - already tagged or picked");
			}
		}

		LogWriter.getLogger().debug(
				"--------------------------------------------------");

		// creates the pruning results file
		//createResultsFile();
	}

//	/**
//	 * Creates output file for pruning results. The file will be placed in file
//	 * path specified by the user. It will overwrite already existing file with
//	 * same name.
//	 * 
//	 * @throws PriorityPrunerException
//	 *             if new file couldn't be created
//	 */
//	private void createResultsFile() throws PriorityPrunerException {
//
//		BufferedWriter writer = null;
//		
//		try {
//			// creates output file
//			//DecimalFormat df = new DecimalFormat("0.00##");
//			
//			File outputFile = new File(options.getOutputPrefix() + ".results");
//			writer = new BufferedWriter(new FileWriter(outputFile));
//			
//			LogWriter.getLogger().info("Writing pruning results to [ " + options.getOutputPrefix() + ".results" + " ]");
//			// writes to log and output files
//			
//			writer.write("name" + "\t" + "chr" + "\t" + "pos" + "\t" + "a1" + "\t" + "a2" + "\t" 
//						+ "tagged" + "\t" + "selected" + "\t" + "best_tag" + "\t" + "r^2" + "\n" );
//			for (SnpInfo snp : snpListFile.getSnps()) {
//						
//				String bestTag = "NA";
//				String r2 = "NA";
//				if (snp.getTaggedByList().size() > 0){
//					Double bestR2 = new Double(-1);
//					for (SnpR2Pair tag: snp.getTaggedByList()){
//						if (tag.getrSquared() > bestR2){
//							bestR2 = tag.getrSquared();
//							bestTag = tag.getSnp().getSnpName();
//						}
//					}
//					r2 = bestR2.toString();
//					//primaryTag = snp.getTaggedByList().
//				}
//				
//				if (snp.getSnpGenotypes().isValid() || snp.getForceInclude()){
//					writer.write(snp.getSnpName() + "\t"
//							+ snp.getChr() + "\t" 
//							+ snp.getPos() + "\t"
//							+ snp.getAllele1() + "\t" 
//							+ snp.getAllele2() + "\t"
//							+ (snp.getTagged() ? "1" : "0") + "\t" 
//							+ (snp.getPicked() ? "1" : "0") + "\t"
//							+ bestTag + "\t" 
//							+ r2 + "\n");
//				}
//			}
//		} catch (IOException e) {
//			throw new PriorityPrunerException("Could not create file: "
//					+ e.getMessage()
//					+ "\nPlease check that correct file path is provided.");
//		}finally{
//			if (writer != null){
//				try{
//					writer.close();
//				}catch(IOException e){e.printStackTrace();}
//			}
//		}
//	}

	/**
	 * Calculates scaled scores for potential surrogates, based on their minimum
	 * and maximum values of chosen metrics.
	 * 
	 * @param potentialSurrogates
	 *            potential surrogates for current SNP, represented as
	 *            Result-objects
	 */
	private void recalculateScores(ArrayList<Result> potentialSurrogates) {

		// for each metric, find min and max
		double[] metricMins = new double[options.getMetrics().size()];
		double[] metricMaxes = new double[options.getMetrics().size()];

		for (int i = 0; i < options.getMetrics().size(); i++) {
			double min = Double.MAX_VALUE;
			double max = Double.MIN_VALUE;

			for (Result result : potentialSurrogates) {
				if (result.getPartnerSnp().getMetrics()[i] > max) {
					max = result.getPartnerSnp().getMetrics()[i];
				}
				if (result.getPartnerSnp().getMetrics()[i] < min) {
					min = result.getPartnerSnp().getMetrics()[i];
				}
			}
			metricMins[i] = min;
			metricMaxes[i] = max;
		}
		// now calculate the scaled score for each SNP
		for (Result result : potentialSurrogates) {
			double score = 0;
			for (int i = 0; i < options.getMetrics().size(); i++) {
				if (metricMaxes[i] == metricMins[i]) {
					score += options.getMetrics().get(i).getWeight();
				} else {
					score += options.getMetrics().get(i).getWeight()
							* ((result.getPartnerSnp().getMetrics()[i] - metricMins[i]) / (metricMaxes[i] - metricMins[i]));
				}
			}
			result.getPartnerSnp().setScore(score);
		}
	}

	/**
	 * Method responsible for the pruning procedure. Based on current index SNP
	 * and window size, a list of SNPs to calculate LD with is created. Together
	 * with index SNP and thresholds for maf, hwe, and genotype percentage, LD
	 * calculation is initiated by this method, and performed by the
	 * SnpWorkUnit. According to r^2 value, design score and current thresholds,
	 * both already picked and potential surrogates are identified and stored.
	 * 
	 * If the index SNP doesn't pass the design score minimum, the default is
	 * that it's skipped and a new index SNP gets chosen. If the surrogate
	 * option "-s" is entered in the command line, an additional surrogate gets
	 * added instead. If index SNP pass design score minimum, it gets picked and
	 * tagged.
	 * 
	 * The scores of the potential surrogates are calculated and sorted, and
	 * depending on the number of surrogates still needed, surrogates get picked
	 * and tagged. Finally, SNPs within the window are tagged if their r^2 value
	 * is equal to or above the current r^2 threshold.
	 * 
	 * @param indexSnp
	 *            current index SNP
	 * @throws PriorityPrunerException
	 *             if problems are encountered during LD calculation
	 * @throws IOException 
	 */
	private void prune(SnpInfo indexSnp) throws PriorityPrunerException {

		// first check if the index SNP passes the design score or is force-included
		if (indexSnp.getDesignScore() < options.getMinDesignScore() && 
				!indexSnp.getForceInclude()){
			LogWriter.getLogger().debug("Skipping index SNP " + indexSnp.getSnpName() + 
					"- design score is less than threshold.");
			return;
		}
		
		// check if the index SNP passes maf, hwe and call rate thresholds
		if (!indexSnp.getSnpGenotypes().isValid() && !indexSnp.getForceInclude()){
			LogWriter.getLogger().debug("Skipping index SNP " + indexSnp.getSnpName() + 
					"- does not pass MAF, HWE or call rate threshold.");
			return;
		}
		
		int startPos = (int) (indexSnp.getPos() - options.getMaxDistance());
		int endPos = (int) (indexSnp.getPos() + options.getMaxDistance());

		if (startPos < 1) {
			startPos = 1;
		}

		SnpInfo startSnpInfo = indexSnp;
		SnpInfo endSnpInfo = indexSnp;

		// find the starting SNP
		while (startSnpInfo.getSortedByPosIndex() > 0) {
			SnpInfo checkSnpInfo = snpListFile.getSnpsSortedByChrPos().get(
					startSnpInfo.getSortedByPosIndex() - 1);
			if (!checkSnpInfo.getChr().equals(indexSnp.getChr())
					|| checkSnpInfo.getPos() < startPos) {
				break;
			} else {
				startSnpInfo = checkSnpInfo;
			}
		}

		// find the ending SNP
		while (endSnpInfo.getSortedByPosIndex() < snpListFile
				.getSnpsSortedByChrPos().size() - 1) {
			SnpInfo checkSnpInfo = snpListFile.getSnpsSortedByChrPos().get(
					endSnpInfo.getSortedByPosIndex() + 1);
			if (!checkSnpInfo.getChr().equals(indexSnp.getChr())
					|| checkSnpInfo.getPos() > endPos) {
				break;
			} else {
				endSnpInfo = checkSnpInfo;
			}
		}

		int referenceSNPIndex = -1;
		ArrayList<SnpGenotypes> genotypesList = new ArrayList<SnpGenotypes>();

		// get the list of SNPs to calculate LD with
		for (int i = startSnpInfo.getSortedByPosIndex(); i < endSnpInfo
				.getSortedByPosIndex() + 1; i++) {
			SnpInfo snpInfo = snpListFile.getSnpsSortedByChrPos().get(i);
			// if genotypes couldn't be found -- this shouldn't ever happen
			if (snpInfo.getSnpGenotypes() == null) {
				throw new PriorityPrunerException("Could not find genotypes for " + 
					snpInfo.getSnpName());
			}
			
			// if the snp is valid or if this is the index snp, add it to the genotype list
			if (snpInfo.getSnpGenotypes().isValid() || indexSnp == snpInfo){
				genotypesList.add(snpInfo.getSnpGenotypes());
			}
			
			if (indexSnp == snpInfo) {
				// the position of the index SNP in the list
				referenceSNPIndex = genotypesList.size() - 1;
			}
			

		}
		// if index SNP wasn't found - this shouldn't ever happen
		if (referenceSNPIndex < 0) {
			throw new PriorityPrunerException("Could not find genotypes for index SNP: " + 
					indexSnp.getSnpName());
		}
		
		
		// calling SnpWorkUnit to do LD calculations
		SnpWorkUnit snpWorkUnit = new SnpWorkUnit(indexSnp.getSnpName(),
				genotypesList, referenceSNPIndex,
				genotypes.getKeptFounders());

		snpWorkUnit.performWork();


		// determine which r^2 threshold to use		
		// if no threshold is defined for associated p-value, an exception gets
		// thrown
		double r2Threshold = -1;
		for (R2Threshold threshold : options.getSortedR2Thresholds()) {
			if (indexSnp.getPValue() <= threshold.getPValue()) {
				r2Threshold = threshold.getR2Threshold();
				LogWriter.getLogger().debug(
						"Defined r^2-threshold: " + r2Threshold);
				break;
			}
		}
		if (r2Threshold < 0) {
			throw new PriorityPrunerException(
					"No r-squared threshold defined for p-value of " + indexSnp.getSnpName() + ": "
							+ indexSnp.getPValue());
		}



		// pick index SNP
		indexSnp.setPicked(true);
		indexSnp.setTagged(true);
		indexSnp.setPickOrder(pickOrder);
		pickOrder++;
		LogWriter.getLogger().debug("Selecting index SNP " + indexSnp.getSnpName());
		
		// pick surrogates if necessary
		pickSurrogates(indexSnp, r2Threshold, snpWorkUnit.getResults());
		
		// tags SNPs within the pruning window if their r^2-value are equal to
		// or above the current r^2-threshold
		int numTagged = 0;
		for (Result result : snpWorkUnit.getResults()) {

			if (this.ldFile != null){
				try{
					// prints to LD output file
					this.ldFile.writeLdRow(indexSnp.getSnpName(),indexSnp.getChr(),indexSnp.getPos(),
							indexSnp.getAllele1(),indexSnp.getAllele2(),result.getPartnerSnpName(), 
							result.getPartnerChr(),result.getPartnerPos(),result.getPartnerSnp().getAllele1(),
							result.getPartnerSnp().getAllele2(),result.getRSquared(),result.getDPrime());
				}catch(IOException e){
					throw new PriorityPrunerException("Could not write to LD table: " + e.getMessage());
				}
			}
			if (result.getRSquared() >= r2Threshold) {
				result.getPartnerSnp().setTagged(true);
				numTagged++;
				result.getPartnerSnp().addTaggedBy(indexSnp,
						result.getRSquared());
			}
		}

		LogWriter.getLogger().debug("Marking "+ numTagged + " SNP(s) as tagged.");
	}
	
	private void pickSurrogates(SnpInfo indexSnp, double r2Threshold, ArrayList<Result> results){
		
		DecimalFormat decimal = new DecimalFormat("##.00");
		
		// determine number of surrogates needed
		// if user specified that surrogates shouldn't be added for force
		// included SNPs - set numSurrogates at 0
		int numSurrogates = 0;
		if (indexSnp.getForceInclude()
				&& !options.getAddSurrogatesForForceIncludedSnps()) {
			numSurrogates = 0;
		}else{
			for (SurrogateThreshold threshold : options.getSortedSurrogateThresholds()) {
				if (indexSnp.getPValue() < threshold.getPValue()) {
					numSurrogates = threshold.getNumSurrogates();
					break;
				}
			}
		}
		
		LogWriter.getLogger().debug("Surrogates needed: " + numSurrogates);
		
		if (numSurrogates == 0){
			return;
		}
		
		// count how many surrogates the index SNP already has (picked by previous
		// index SNPs), and how many potential surrogates that are available.
		// also write LD results to file
		ArrayList<SnpInfo> surrogatesPicked = new ArrayList<SnpInfo>();
		ArrayList<Result> potentialSurrogateResults = new ArrayList<Result>();
		for (int i = 0; i < results.size(); i++) {
			Result result = results.get(i);

			// adds to list of already picked surrogates
			if (result.getPartnerSnp().getPicked()
					&& result.getRSquared() >= r2Threshold
					&& !result.getPartnerSnp().equals(indexSnp)) {
				surrogatesPicked.add(result.getPartnerSnp());
			}
			
			// adds to list of potential surrogates
			if (!result.getPartnerSnp().getPicked()
					&& result.getRSquared() >= r2Threshold
					&& (result.getPartnerSnp().getDesignScore() >= options.getMinDesignScore()
						|| result.getPartnerSnp().getForceInclude())
					&& !result.getPartnerSnp().equals(indexSnp)) {
				potentialSurrogateResults.add(result);
			}
		}
		LogWriter.getLogger().debug("Surrogates available: " + potentialSurrogateResults.size());
		
		
		// calculates scaled scores for potential surrogates and sorts them
		if (options.getMetrics().size() > 0){
			recalculateScores(potentialSurrogateResults);
			Collections.sort(potentialSurrogateResults, new ScoreSorter());
		}
		// sort by r-squared
		else{
			
			Collections.sort(potentialSurrogateResults, new RSquaredSorter());
		}
		
		// now pick additional surrogates
		for (Result result : potentialSurrogateResults) {
			if (surrogatesPicked.size() < numSurrogates) {
				result.getPartnerSnp().setPicked(true);
				result.getPartnerSnp().setTagged(true);
				result.getPartnerSnp().setPickOrder(pickOrder);
				pickOrder++;
				surrogatesPicked.add(result.getPartnerSnp());

				LogWriter.getLogger().debug(
							"Selecting surrogate "
									+ result.getPartnerSnp().getSnpName()
							+ " with r^2: " + decimal.format(result.getRSquared()));
			}
		}

	}
	
	/***
	 * Class to sort Result objects by force include status, then 
	 * pairwise r-squared, then SNP name (to prevent randomness) 
	 *
	 */
	private class RSquaredSorter implements Comparator<Result> {
		@Override
		public int compare(Result x, Result y) {
			if (x.getPartnerSnp().getForceInclude() == y.getPartnerSnp()
					.getForceInclude()) {
						if (x.getRSquared() == y.getRSquared()) {
							return (y.getPartnerSnp().getSnpName().compareTo(x
									.getPartnerSnp().getSnpName()));
						} else {
							return (Double.compare(y.getRSquared(), x.getRSquared()));
						}
			} else {
				return (Boolean.valueOf(y.getPartnerSnp().getForceInclude())
						.compareTo(Boolean.valueOf(x.getPartnerSnp()
								.getForceInclude())));
			}
		}
	}
}
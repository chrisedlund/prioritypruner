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

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;

/**
 * This class is responsible for pruning and prioritizing the SNPs defined in
 * the input files. It initiates the parsing of these files, goes through
 * current SNPs, checks that valid properties are defined, initiates LD
 * calculation, and based on these LD results - picks and tags SNPs. It also
 * creates output files for the LD calculations and pruning results.
 */
public class Pruner {

	private SnpListFile snpListFile;
	private Genotypes genotypes;
	private int pickOrder = 1;
	private CommandLineOptions options = CommandLineOptions.getInstance();
	private PrintStream ldPrintStream;
	private PrintStream printStreamChrOutput;

	/**
	 * Constructor for Pruner. Initiates parsing of genotype and family data by
	 * creating SnpListFile- and TPlink-objects. It creates an output file for
	 * the LD results calculated by SnpWorkUnit and initiates the pruning.
	 * 
	 * @throws PriorityPrunerException
	 *             if problems are encountered during parsing or pruning
	 */
	public Pruner() throws PriorityPrunerException {
		// initiates parsing of the SNP input file
		this.snpListFile = new SnpListFile(options.getSnpFilePath(),
				options.getNumMetrics());
		// initiates parsing of tped and tfam files
		this.genotypes = new TPlink(options.getTped(), options.getTfam(),
				this.snpListFile);
		createLdFile();
		init();
	}

	/**
	 * Creates LD file to store information calculated in SnpWorkUnit. The file
	 * will be placed in file path specified by the user. It will overwrite
	 * already existing file with same name.
	 * 
	 * @throws PriorityPrunerException
	 *             if new file couldn't be created
	 */
	private void createLdFile() throws PriorityPrunerException {

		File ldFile;
		FileOutputStream ldOutputStream;
		try {
			ldFile = new File(options.getOutputFile() + ".ld");
			ldOutputStream = new FileOutputStream(ldFile);

			if (!ldFile.exists()) {
				ldFile.createNewFile();
			}
			ldPrintStream = new PrintStream(ldOutputStream);
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not create file: "
					+ e.getMessage()
					+ "\nPlease check that correct file path is provided.");
		}
		this.ldPrintStream.println("Index SNP	Partner SNP	Distance	r^2	D'");
	}

	/**
	 * Initiates pruning of SNPs defined in the SNP input file. Loops through
	 * provided SNPs and checks that chromosome is valid and that information
	 * about the SNP was also in tped file. Prunes the current SNP if it's not
	 * already picked and either force included or not tagged.
	 * 
	 * @throws PriorityPrunerException
	 *             if invalid information is encountered
	 */
	private void init() throws PriorityPrunerException {

		int prunedSnpIndex = 0;

		// goes through all SNPs in the SNP input file
		for (SnpInfo snp : snpListFile.getSnps()) {
			// checks that no unsupported chromosomes are specified in the SNP
			// input file
			if (snp.getChr().toUpperCase().equals("M")
					|| snp.getChr().toUpperCase().equals("MT")
					|| snp.getChr().toUpperCase().equals("CHRM")
					|| snp.getChr().toUpperCase().equals("CHR_M")
					|| snp.getChr().toUpperCase().equals("Y")
					|| snp.getChr().toUpperCase().equals("CHRY")
					|| snp.getChr().toUpperCase().equals("CHR_Y")
					|| snp.getChr().toUpperCase().equals("24")) {
				throw new PriorityPrunerException(
						"Chromosome \""
								+ snp.getChr()
								+ "\" is not supported. \nPlease update input files and rerun program");
			}
			// checks that corresponding SNP is provided in tped file
			if (!snp.getInTped()) {
				throw new PriorityPrunerException("SNP \"" + snp.getSnpName()
						+ "\" at chromsome " + snp.getChr()
						+ ", base pair position " + snp.getPos()
						+ " with alleles: " + snp.getAllele1() + " "
						+ snp.getAllele2() + ", couldn't be found in file \""
						+ options.getTped()
						+ "\".\nPlease update input files and rerun program.");
			}
			// prunes current SNP if it's not already picked and either force
			// included or not tagged
			if (!snp.getPicked() && (!snp.getTagged() || snp.getForceInclude())) {
				LogWriter.getLogger().debug(
						"\nUsing index SNP \"" + snp.getSnpName()
								+ "\" with p-value: " + snp.getPValue());
				prune(snp);
				prunedSnpIndex++;

				// prints out how many SNPs that currently are processed, if
				// verbose option is entered
				if (prunedSnpIndex % 1000 == 0) {
					LogWriter.getLogger().debug(
							"\nProcessed " + prunedSnpIndex + " SNPs\n");
				}
			} else {
				LogWriter.getLogger().debug(
						"\nSkipping SNP \"" + snp.getSnpName()
								+ "\" - already tagged/picked");
			}
		}
		ldPrintStream.close();
		LogWriter.getLogger().debug(
				"--------------------------------------------------");

		// creates the pruning results file
		createReslutFile();
	}

	/**
	 * Creates output file for pruning results. The file will be placed in file
	 * path specified by the user. It will overwrite already existing file with
	 * same name.
	 * 
	 * @throws PriorityPrunerException
	 *             if new file couldn't be created
	 */
	private void createReslutFile() throws PriorityPrunerException {

		File outputFile;
		FileOutputStream outputStream;
		try {
			// creates output file
			outputFile = new File(options.getOutputFile() + ".results");
			outputStream = new FileOutputStream(outputFile);
			if (!outputFile.exists()) {
				outputFile.createNewFile();
			}
			printStreamChrOutput = new PrintStream(outputStream);
			// writes to log and output files
			printStreamChrOutput
					.println("SNP	Chr	Pos		A1	A2	P-value		BeadTypes Tagged Picked	PickOrder");
			LogWriter.getLogger().debug(
					"\nSNP	Chr	Pos	A1	A2	Score	Tagged	Picked");
			for (SnpInfo snp : snpListFile.getSnps()) {
				LogWriter.getLogger().debug(
						snp.getSnpName() + "\t" + snp.getChr() + "\t"
								+ snp.getPos() + "\t" + snp.getAllele1() + "\t"
								+ snp.getAllele2() + "\t" + snp.getScore()
								+ "\t" + snp.getTagged() + "\t"
								+ snp.getPicked());
				printStreamChrOutput.println(snp.getSnpName() + "\t"
						+ snp.getChr() + "\t" + snp.getPos() + "\t"
						+ snp.getAllele1() + "\t" + snp.getAllele2() + "\t"
						+ snp.getPValue() + "\t\t" + snp.getNumBeadTypes()
						+ "\t" + snp.getTagged() + "\t" + snp.getPicked()
						+ "\t" + snp.getPickOrder());
			}
			printStreamChrOutput.close();
		} catch (IOException e) {
			throw new PriorityPrunerException("Could not create file: "
					+ e.getMessage()
					+ "\nPlease check that correct file path is provided.");
		}
	}

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
	 */
	private void prune(SnpInfo indexSnp) throws PriorityPrunerException {

		int startPos = (int) (indexSnp.getPos() - options.getHalfWindowSize());
		int endPos = (int) (indexSnp.getPos() + options.getHalfWindowSize());

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
			if (snpInfo.getSnpGenotypes() == null) {
				continue;
			}
			genotypesList.add(snpInfo.getSnpGenotypes());
			if (indexSnp == snpInfo) {
				// the position of the index SNP in the list
				referenceSNPIndex = genotypesList.size() - 1;
			}
		}
		// if index SNP wasn't found
		if (referenceSNPIndex < 0) {
			if (indexSnp.getDesignScore() >= options
					.getAbsoluteMinDesignScore()) {
				indexSnp.setTagged(true);
				indexSnp.setPicked(true);
				indexSnp.setPickOrder(pickOrder);
				pickOrder++;
				LogWriter.getLogger().debug(
						"Could not find index SNP \"" + indexSnp.getSnpName()
								+ "\" at chr \"" + indexSnp.getChr() + ":"
								+ indexSnp.getPos()
								+ "\" in input files. Picking anyways.");
			} else {
				LogWriter
						.getLogger()
						.debug("Could not find index SNP \""
								+ indexSnp.getSnpName()
								+ "\" at chr \""
								+ indexSnp.getChr()
								+ ":"
								+ indexSnp.getPos()
								+ "\" in input files. It could not be picked due to low design score.");
			}
			return;
		}
		// calling SnpWorkUnit to do LD calculations
		SnpWorkUnit snpWorkUnit = new SnpWorkUnit(indexSnp.getSnpName(),
				genotypesList, referenceSNPIndex,
				genotypes.getFounderIndices(), genotypes.getSubjectSexes(),
				options.getMinMaf(), options.getMinimumHardyWeinbergPvalue(),
				options.getMinimumGenotypePercentage());

		snpWorkUnit.performWork();

		// if index SNP didn't pass all filters for MAF, HWE and missing
		// genotype percentage - skip it, and choose new index SNP
		if (!snpWorkUnit.getIndexSnpPassed() && !indexSnp.getForceInclude()) {
			LogWriter
					.getLogger()
					.debug("Skipping SNP \""
							+ indexSnp.getSnpName()
							+ "\", since it didn't pass all filters for MAF, HWE and missing genotype percentage.");
			return;
		}

		// checks the number of surrogates needed
		int numSurrogates = 0;
		for (SurrogateThreshold threshold : options
				.getSortedSurrogateThresholds()) {
			if (indexSnp.getPValue() < threshold.getPValue()) {
				numSurrogates = threshold.getNumSurrogates();
				break;
			}
		}

		// checks which r^2 threshold to use
		double r2Threshold = -1;
		for (R2Threshold threshold : options.getSortedR2Thresholds()) {
			if (indexSnp.getPValue() <= threshold.getPValue()) {
				r2Threshold = threshold.getR2Threshold();
				LogWriter.getLogger().debug(
						"Defined r^2-threshold: " + r2Threshold);
				break;
			}
		}
		// if no threshold is defined for associated p-value, an exception gets
		// thrown
		if (r2Threshold == -1) {
			throw new PriorityPrunerException(
					"No r^2-threshold defined for p-values equal to, or above: "
							+ indexSnp.getPValue()
							+ ".\nTo use the p-value dependent r^2-threshold option, please define thresholds valid for all p-values in the range. \nFor more information about r^2 options, type \"-h\".");
		}

		// first count how many surrogates we already have (picked by previous
		// index SNPs), and how many potential surrogates that are available
		ArrayList<SnpInfo> surrogatesPicked = new ArrayList<SnpInfo>();
		ArrayList<Result> potentialSurrogateResults = new ArrayList<Result>();

		for (int i = 0; i < snpWorkUnit.getResults().size(); i++) {
			Result result = snpWorkUnit.getResults().get(i);

			// prints to LD output file
			this.ldPrintStream.println(indexSnp.getSnpName()
					+ "	"
					+ result.getPartnerSnpName()
					+ "	"
					+ Math.abs(indexSnp.getPos()
							- result.getPartnerSnp().getPos()) + "	"
					+ result.getRSquared() + "	" + result.getDPrime());

			// adds to list of already picked surrogates
			if (result.getPartnerSnp().getPicked()
					&& result.getRSquared() >= r2Threshold
					&& !result.getPartnerSnp().equals(indexSnp)) {
				surrogatesPicked.add(result.getPartnerSnp());
			}
			// adds to list of potential surrogates
			if (!result.getPartnerSnp().getPicked()
					&& result.getRSquared() >= r2Threshold
					&& result.getPartnerSnp().getDesignScore() >= options
							.getAbsoluteMinDesignScore()
					&& !result.getPartnerSnp().equals(indexSnp)) {
				potentialSurrogateResults.add(result);
			}
		}
		// decides if index SNP as well as additional surrogates will be picked

		// if we can't pick the index SNP due to low design score, add
		// an additional surrogate instead (only if the
		// "AddOneAdditionalSurrogate-option" is entered in the command line,
		// the default is not to add an extra surrogate)
		if (indexSnp.getDesignScore() < options.getAbsoluteMinDesignScore()
				&& options.getOneAdditionalSurrogate()) {
			numSurrogates++;

			// if user specified that surrogates shouldn't be added for force
			// included SNPs but still added the
			// "AddOneAdditionalSurrogate-option", the later option will be
			// overridden, and no additional surrogates will be added for force
			// included SNPs
			if (indexSnp.getForceInclude()
					&& !options.getAddSurrogatesForForceIncludedSnps()) {

				numSurrogates = 0;
				LogWriter.getLogger().debug(
						"Surrogates needed: " + numSurrogates);
				LogWriter
						.getLogger()
						.debug("SNP \""
								+ indexSnp.getSnpName()
								+ "\" didn't pass the absolute minimum design score requirement. No surrogates are added since this SNP is force included.\nChoosing new index SNP.");
				return;
			} else {
				LogWriter.getLogger().debug(
						"Surrogates needed: " + numSurrogates);
				LogWriter
						.getLogger()
						.debug("SNP \""
								+ indexSnp.getSnpName()
								+ "\" didn't pass the absolute minimum design score requirement. Choosing an extra surrogate instead.");
			}

			// if user specified that no additional surrogate should be
			// added in case the index SNP didn't pass design score minimum
			// - skip this one, and choose new index SNP instead
		} else if (indexSnp.getDesignScore() < options
				.getAbsoluteMinDesignScore()
				&& !options.getOneAdditionalSurrogate()) {
			LogWriter.getLogger().debug("Surrogates needed: " + numSurrogates);
			LogWriter
					.getLogger()
					.debug("SNP \""
							+ indexSnp.getSnpName()
							+ "\" didn't pass the absolute minimum design score requirement. Choosing new index SNP.");
			return;
		} else {
			// if user specified that surrogates shouldn't be added for force
			// included SNPs - set numSurrogates to 0
			if (indexSnp.getForceInclude()
					&& !options.getAddSurrogatesForForceIncludedSnps()) {
				numSurrogates = 0;
			}
			LogWriter.getLogger().debug("Surrogates needed: " + numSurrogates);

			// picks index SNP
			if (!indexSnp.getPicked()) {
				indexSnp.setPicked(true);
				indexSnp.setTagged(true);
				indexSnp.setPickOrder(pickOrder);
				pickOrder++;
				LogWriter.getLogger().debug(
						"Picking index SNP \"" + indexSnp.getSnpName() + "\"");
			}
		}
		// calculates scaled scores for potential surrogates and sorts them
		recalculateScores(potentialSurrogateResults);
		Collections.sort(potentialSurrogateResults, new ScoreSorter());
		// now pick additional surrogates
		for (Result result : potentialSurrogateResults) {
			if (surrogatesPicked.size() < numSurrogates) {
				result.getPartnerSnp().setPicked(true);
				result.getPartnerSnp().setTagged(true);
				result.getPartnerSnp().setPickOrder(pickOrder);
				pickOrder++;
				surrogatesPicked.add(result.getPartnerSnp());
				if (result.getPartnerSnp().getForceInclude()) {
					LogWriter.getLogger().debug(
							"Picking surrogate \""
									+ result.getPartnerSnp().getSnpName()
									+ "\" with scaled score: "
									+ result.getPartnerSnp().getScore()
									+ " (force included SNP)");
					LogWriter.getLogger().debug(
							"r^2: " + result.getRSquared() + " and D': "
									+ result.getDPrime());
				} else {
					LogWriter.getLogger().debug(
							"Picking surrogate \""
									+ result.getPartnerSnp().getSnpName()
									+ "\" with scaled score: "
									+ result.getPartnerSnp().getScore());
					LogWriter.getLogger().debug(
							"r^2: " + result.getRSquared() + " and D': "
									+ result.getDPrime());
				}
			}
		}

		// tags SNPs within the pruning window if their r^2-value are equal to
		// or above the current r^2-threshold

		// if index SNP didn't get picked, the partner SNP is marked as tagged
		// by the first surrogate picked. Since r^2-value is unknown in this
		// case, a default value of -1 is entered in the "TaggedBy-list" to
		// indicate this
		if (!indexSnp.getPicked() && surrogatesPicked.size() > 0) {
			for (Result result : snpWorkUnit.getResults()) {
				if (result.getRSquared() >= r2Threshold
						&& !result.getPartnerSnp().equals(indexSnp)) {
					result.getPartnerSnp().setTagged(true);
					result.getPartnerSnp().addTaggedBy(surrogatesPicked.get(0),
							-1);
				}
			}
			// if index SNP got picked, the SNP is marked as tagged by it
		} else if (indexSnp.getPicked()) {
			for (Result result : snpWorkUnit.getResults()) {
				if (result.getRSquared() >= r2Threshold) {
					result.getPartnerSnp().setTagged(true);
					result.getPartnerSnp().addTaggedBy(indexSnp,
							result.getRSquared());
				}
			}
		}
	}
}
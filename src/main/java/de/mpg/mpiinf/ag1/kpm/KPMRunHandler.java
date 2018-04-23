package de.mpg.mpiinf.ag1.kpm;

import de.mpg.mpiinf.ag1.kpm.parsers.PriorityFileParser;
import de.mpg.mpiinf.ag1.kpm.utils.Parser;
import de.mpg.mpiinf.ag1.kpm.utils.StatisticsUtility;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.charts.IChart;
import dk.sdu.kpm.results.IKPMResultSet;
import dk.sdu.kpm.results.IKPMRunListener;
import dk.sdu.kpm.runners.BatchRunWithPerturbationParameters;
import dk.sdu.kpm.runners.BatchRunWithPerturbationRunner;
import dk.sdu.kpm.runners.BatchRunner;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class that handles the actual KPM runs, and ensures we get the results.
 * Assumes all other parameters has already been set.
 * @author Martin
 *
 */
public class KPMRunHandler implements IKPMRunListener{

    private volatile KPMSettings kpmSettings;

    private volatile Parameters params;

	public KPMRunHandler(KPMSettings settings){
		this.kpmSettings = settings;
	}

	/**
	 * Starts the runs, ensures the settings are correctly set up. Assumes all other parameters has already been set.
	 * @param params
	 */
	public void runBatch(Parameters params){
		this.params = params;
		// Register starting time
		long start, end, time;
		kpmSettings.STARTING_TIME = System.nanoTime();
		start = System.currentTimeMillis();



		System.out.println("\n*********** CREATING GRAPH ***************\n");
		// Parse the graph and matrix files
		params.PARSER = new Parser(params.GRAPH_FILE, kpmSettings.MATRIX_FILES_MAP, Parser.TAB);

		// Create KPM graph
		params.INPUT_GRAPH = params.PARSER.createGraph(kpmSettings);

		// Print graph statistics
		if (params.PRINT_GRAPH_STATS) {
			System.out.println("\n NETWORK STATS: \n");
			StatisticsUtility.printGraphStatistics(params.INPUT_GRAPH);
		}

		// Print dataset statistics
		if (params.PRINT_DATASETS_STATS) {
			System.out.println("\n DATASETS STATS: \n");
			StatisticsUtility.printMatrixStatistics(params.PARSER, params.INPUT_GRAPH);
		}

		// We set the positive and negative lists of genes in the graph
		if (params.POSITIVE_FILE != null) {
			Set<String> positiveList = PriorityFileParser.parse(params.POSITIVE_FILE);
			params.INPUT_GRAPH.setPositiveList(positiveList);
		}

		if (params.NEGATIVE_FILE != null) {
			Set<String> negativeList = PriorityFileParser.parse(params.NEGATIVE_FILE);
			params.INPUT_GRAPH.setNegativeList(negativeList);
		}

		// Add validation to the search, if any file was set.
		try{
			if(!params.VALIDATION_FILE.isEmpty()){
				params = setValidationFile(params);
			}
		}catch(FileNotFoundException e){
			System.out.println(String.format("\nValidation file could not be found: '%s'.\n\nKPM will now terminate.", params.VALIDATION_FILE));
			return;
		}catch(IOException e){
			System.out.println(String.format("\nProblems reading validation file: '%s'.\n\nKPM will now terminate.", params.VALIDATION_FILE));
			return;
		}

		// IMPORTANT:  Graph must be refreshed before each new run of an algorithm
		params.INPUT_GRAPH.refreshGraph(kpmSettings);

		kpmSettings.MAIN_GRAPH = params.INPUT_GRAPH;//new KPMGraph(params.INPUT_GRAPH);
        kpmSettings.INCLUDE_CHARTS = kpmSettings.IS_BATCH_RUN;

		// Check which strategy and algorithm to use

		switch (kpmSettings.ALGO) {
		case GREEDY:
			params.STRATEGY =  KPMStrategy.INES;
			params.ALGORITHM = KPMAlgorithm.GREEDY;
			break;

		case LCG:
			params.STRATEGY =  KPMStrategy.INES;
			params.ALGORITHM = KPMAlgorithm.ACO;
			break;

		case OPTIMAL:
			params.STRATEGY =  KPMStrategy.INES;
			params.ALGORITHM = KPMAlgorithm.OPTIMAL;
			break;

		case EXCEPTIONSUMGREEDY:
			params.STRATEGY =  KPMStrategy.GLONE;
			params.ALGORITHM = KPMAlgorithm.GREEDY;
			break;

		case EXCEPTIONSUMACO:
			params.STRATEGY = KPMStrategy.GLONE;
			params.ALGORITHM = KPMAlgorithm.ACO;
			break;

		case EXCEPTIONSUMOPTIMAL:
			params.STRATEGY = KPMStrategy.GLONE;
			params.ALGORITHM = KPMAlgorithm.OPTIMAL;
			break;
		}

		// The configurator also sets the values in the KPMParameters-class.
		KPMParameterConfigurator configurator = new KPMParameterConfigurator(params, kpmSettings);
		if(!configurator.areSettingsValid()){
			System.out.println("Settings are invalid. Stopping.");
			return;
		}



		System.out.println("\n********* RUNNING " + params.STRATEGY + " (" + params.ALGORITHM + ") " + "... *************\n");
		if(params.IS_PERTURBATION_RUN){
			runBatchWithPertubation(params);
		}else{
			runStandard(params);
		}

		end = System.currentTimeMillis();
		time = (end - start) / 1000;
		// Do stats and post-processing stuff after obtaining the results
		System.out.println("TOTAL RUNNING TIME: " + time + " seconds");

		kpmSettings.TOTAL_RUNNING_TIME = time;
	}

	private void runBatchWithPertubation(Parameters params){
		this.params = params;
		KPMStandaloneTaskMonitor monitor = new KPMStandaloneTaskMonitor();
		try {

			System.out.println("params.GRAPHS_PER_STEP = " + params.GRAPHS_PER_STEP);

			String formattedTime = new SimpleDateFormat("HH.mm.ss").format(Calendar.getInstance().getTime());
			BatchRunWithPerturbationParameters bprp = new BatchRunWithPerturbationParameters(
					params.PERTURBATION,
					params.MIN_PERCENTAGE,
					params.STEP_PERCENTAGE,
					params.MAX_PERCENTAGE,
					params.GRAPHS_PER_STEP,
					params.STRATEGY == KPMStrategy.INES,
					String.format("kpm_run_%s", formattedTime),
					monitor,
					this,
					true);
			BatchRunWithPerturbationRunner batcher = new BatchRunWithPerturbationRunner(bprp, kpmSettings);
			batcher.run();

			System.out.println("Finished runs.");
			System.out.println("Started generating results:");

		} catch (Exception e) {
			Logger.getLogger(KPMRunHandler.class.getName()).log(Level.SEVERE, null, e);
			System.out.println("An error occurred during the KPM run with perturbed graphs. Exiting.");
		}
	}

	private void runStandard(Parameters params){
		KPMStandaloneTaskMonitor monitor = new KPMStandaloneTaskMonitor();
		try {
			kpmSettings.USE_INES = params.STRATEGY == KPMStrategy.INES;
			BatchRunner batcher = new BatchRunner("standalone", monitor, this, kpmSettings);
			batcher.run();
		} catch (Exception e) {
			Logger.getLogger(KPMRunHandler.class.getName()).log(Level.SEVERE, null, e);
			System.out.println("An error occurred during the KPM run. Exiting.");
		}
	}

	private Parameters setValidationFile(Parameters params) throws IOException {

			BufferedReader br = new BufferedReader(new FileReader(params.VALIDATION_FILE));
			String line;

			List<String> validationList = new ArrayList<String>();
			while ((line = br.readLine()) != null) {
				validationList.add(line.trim());
			}
			br.close();

			kpmSettings.VALIDATION_GOLDSTANDARD_NODES = validationList;

		return params;
	}

	@Override
	public void runCancelled(String reason, String runID) {
		System.out.println("The run was cancelled. \n" + reason);
	}

	@Override
	public void runFinished(IKPMResultSet results) {
			int counter = 1;

		File resultsFolder = new File(params.RESULTS_FOLDER);
		if(resultsFolder.isDirectory()){


			File resultsChartsFolder = new File(params.RESULTS_FOLDER + "\\charts\\");

			if(!resultsChartsFolder.exists()){
				resultsChartsFolder.mkdir();
			}

			String formattedChartTime = new SimpleDateFormat("yyyy-MM-dd HH.mm.ss").format(Calendar.getInstance().getTime());

			File resultsChartsCurrentFolder = new File(params.RESULTS_FOLDER + "\\charts\\" + formattedChartTime + "\\");

			if(!resultsChartsCurrentFolder.exists()){
				resultsChartsCurrentFolder.mkdir();
			}

			for(IChart chart : results.getCharts().values()){
				String chartFileName = resultsChartsCurrentFolder.getAbsolutePath() + "\\" + String.format("chart%d.png", counter);
				System.out.println(chartFileName);
				File chartsFile = new File(chartFileName);
				try {
					chart.saveAsPng(chartsFile, 1600, 1200);
				} catch (IOException e) {
					e.printStackTrace();
				}
				counter ++;
			}

			if(counter > 0){
				System.out.println(" - Saved charts to '"+resultsChartsFolder.getAbsolutePath()+"'");
			}
		}else{
			System.out.println(" - No path found, could not save charts.");
		}
	}
}

package de.mpg.mpiinf.ag1.kpm;

import EDU.oswego.cs.dl.util.concurrent.FJTask;
import de.mpg.mpiinf.ag1.kpm.parsers.PriorityFileParser;
import de.mpg.mpiinf.ag1.kpm.utils.OutputSettings;
import de.mpg.mpiinf.ag1.kpm.utils.Parser;
import de.mpg.mpiinf.ag1.kpm.utils.StatisticsUtility;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.charts.IChart;
import dk.sdu.kpm.results.IKPMResultSet;
import dk.sdu.kpm.results.IKPMRunListener;
import dk.sdu.kpm.runners.BatchRunWithPerturbationParameters;
import dk.sdu.kpm.runners.BatchRunWithPerturbationRunner;
import dk.sdu.kpm.runners.BatchRunner;
import dk.sdu.kpm.runners.ProbabilisticRunner;

import java.io.*;
import java.nio.file.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.time.*;


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
	 * //@param paramsr (removed due to conflict)
	 */
	public void runBatch(Parameters params){
		this.params = params;
		// Register starting time
		long start, end, time;
		kpmSettings.STARTING_TIME = System.nanoTime();
		start = System.currentTimeMillis();



		System.out.println("\n*********** CREATING GRAPH ***************\n");
		// Parse the graph and matrix files
		params.PARSER = new Parser(params.GRAPH_FILE, kpmSettings.MATRIX_FILES_MAP, Parser.TAB,
                kpmSettings.PVALUE_MAP, params.IS_BINARY_MATRIX, kpmSettings.COMPARATOR, kpmSettings.PVALUE_FILES_MAP, kpmSettings.USE_DOUBLE_VALUES);

		//Test
		params.PARSER2 = new Parser(params.GRAPH_FILE, kpmSettings.MATRIX_FILES_MAP, Parser.TAB,
				kpmSettings.PVALUE_MAP, params.IS_BINARY_MATRIX, kpmSettings.COMPARATOR, kpmSettings.PVALUE_FILES_MAP, kpmSettings.USE_DOUBLE_VALUES);

		// Create KPM graph
		params.INPUT_GRAPH = params.PARSER.createGraph(kpmSettings);

		//Test
		params.INPUT_GRAPH2 = params.PARSER2.createGraph(kpmSettings);

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

            case FDR:
                params.STRATEGY = KPMStrategy.FDR;
                params.ALGORITHM = KPMAlgorithm.FDR;
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
		}
		else if(params.ALGORITHM == KPMAlgorithm.FDR){
		    runFDR(params);
        }
		else{
			runStandard(params);
		}

		end = System.currentTimeMillis();
		time = (end - start) / 1000;
		// Do stats and post-processing stuff after obtaining the results
		System.out.println("TOTAL RUNNING TIME: " + time + " seconds");

		kpmSettings.TOTAL_RUNNING_TIME = time;
	}

	private void runFDR(Parameters params){
        KPMStandaloneTaskMonitor monitor = new KPMStandaloneTaskMonitor();
	    ProbabilisticRunner runner = new ProbabilisticRunner("standalone", monitor, this, kpmSettings, params.INPUT_GRAPH2, params.RESULTS_FOLDER);
	    runner.run();
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
		LocalDateTime now = LocalDateTime.now();
			int counter = 1;

        File resultsFolder;
		try{
		     Files.createDirectories(Paths.get(params.RESULTS_FOLDER));
        } catch (IOException e){
		    e.printStackTrace();
        }

		if(Files.exists(Paths.get(params.RESULTS_FOLDER))){

		    Path resultsChartsFolder = Paths.get(params.RESULTS_FOLDER + File.separator+ "charts");
				try {
					Files.createDirectories(resultsChartsFolder);
				}
				catch (IOException ioe){
                    ioe.printStackTrace();
				}

			//String formattedChartTime = new SimpleDateFormat("yyyy-MM-dd HH.mm.ss").format(Calendar.getInstance().getTime());

			File resultsChartsCurrentFolder = new File(params.RESULTS_FOLDER + File.separator+"charts"+File.separator + now + File.separator);

			if(!resultsChartsCurrentFolder.exists()){
				resultsChartsCurrentFolder.mkdir();
			}

			for(IChart chart : results.getCharts().values()){
				String chartFileName = resultsChartsCurrentFolder.getAbsolutePath() + File.separator + String.format("chart%d.png", counter);
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
				System.out.println(" - Saved charts to '"+resultsChartsFolder.toString()+"'");
			}
		}else{
			System.out.println(" - No path found, could not save charts.");
		}

		// saving result tables to files

        Path resultTableFolder = Paths.get(params.RESULTS_FOLDER+File.separator+"tables"+File.separator+results.getKpmSettings().ALGO
				+params.STRATEGY+kpmSettings.getKpmRunID()+now);
		try{
		    Files.createDirectories(resultTableFolder);
        }catch (IOException ioe){
		    ioe.printStackTrace();
        }


        if(OutputSettings.GENERATE_SUMMARY_FILE){ //(Done)
            String resultSummary = resultTableFolder.toString()+File.separator+ OutputSettings.SUMMARY_FILE;
            StatisticsUtility.writeSummaryFile(resultSummary, results, kpmSettings);
        }
        if(OutputSettings.GENERATE_DATASETS_STATS_FILE){ //(Done)
            String resultDatasetStats = resultTableFolder.toString()+File.separator+OutputSettings.DATASETS_STATS_FILE;
            StatisticsUtility.writeDatasetsStats(resultDatasetStats, params.PARSER, kpmSettings.MAIN_GRAPH );
        }

        if(OutputSettings.GENERATE_GENE_STATS_FILE){ //Done
            String resultGeneStats = resultTableFolder.toString()+File.separator+OutputSettings.GENE_STATS_FILE;
            StatisticsUtility.writeGeneStatsFile(resultGeneStats, results, kpmSettings.MAIN_GRAPH, kpmSettings);
        }

        if(OutputSettings.GENERATE_PATHWAYS_FILE){ // Done
            String resultIndividualPathwaysFile= resultTableFolder.toString();
            StatisticsUtility.writeIndividualPathwayFiles(resultIndividualPathwaysFile, results, kpmSettings.MAIN_GRAPH, params);
            String resultPathwaysFile= resultTableFolder.toString()+File.separator+OutputSettings.PATHWAYS_FILE;
            StatisticsUtility.writePathwaysFile(resultPathwaysFile, results, kpmSettings.MAIN_GRAPH, kpmSettings);
        }

        if(OutputSettings.GENERATE_PATHWAYS_STATS_FILE) { //Done
            String resultPathwayStatsFile = resultTableFolder.toString() + File.separator + OutputSettings.PATHWAYS_STATS_FILE;
            StatisticsUtility.writePathwaysStatsFile(resultPathwayStatsFile, results, kpmSettings.MAIN_GRAPH);
        }

       // String resultFile2= resultTableFolder.toString()+File.separator+"results2.txt";
        //StatisticsUtility.writeResultsToFile2(resultFile2);
        //String resultMatrixStatisics= resultTableFolder.toString()+File.separator+"resultSummary.txt";
        StatisticsUtility.printMatrixStatistics(params.PARSER, kpmSettings.MAIN_GRAPH);
        //String resultGraphicsStatistics= resultTableFolder.toString()+File.separator+"resultSummary.txt";
        //StatisticsUtility.printGraphStatistics(resultGraphicsStatistics,);
        //StatisticsUtility.generateHistogramLValues();
	}
}

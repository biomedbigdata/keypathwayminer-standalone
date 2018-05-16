/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.main;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import de.mpg.mpiinf.ag1.kpm.KPMRunHandler;
import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.Program;
import de.mpg.mpiinf.ag1.kpm.parsers.ArgsParametersParser;
import de.mpg.mpiinf.ag1.kpm.parsers.InputFileParser;
import de.mpg.mpiinf.ag1.kpm.parsers.PropertiesFileParser;
import de.mpg.mpiinf.ag1.kpm.shortestpath.ShortestPathAlgorithms;
import dk.sdu.kpm.KPMSettings;

/**
 *
 * @author nalcaraz
 */
public class Main {

    private volatile KPMSettings kpmSettings;
    
	public static void main(String[] args) throws Exception {
		Main main = new Main();
		main.run(args);
	}

	public void run(String[] args){
        int mb = 1024*1024;
        Runtime runtime = Runtime.getRuntime();
        System.out.println("Used Memory:" + (runtime.freeMemory()) / mb);

        // Call kmpDefaultSettings
		kpmSettings = new KPMSettings();


		// parse parameters from the kpm.properties file
		PropertiesFileParser pfp = new PropertiesFileParser(kpmSettings);
		Parameters params = pfp.parse("src/main/resources/kpm.properties");

		// Setting folder prefix for testing:
		String prefix = "src/main/resources/";
		//Parameters params = new Parameters();
		/*params.DATASETS_FILE = prefix + params.DATASETS_FILE;
		params.POSITIVE_FILE = prefix + params.POSITIVE_FILE;
		params.NEGATIVE_FILE = prefix + params.NEGATIVE_FILE;
		params.RESULTS_FOLDER = prefix + params.RESULTS_FOLDER;
		params.VALIDATION_FILE = prefix + params.VALIDATION_FILE;
		kpmSettings.IS_BATCH_RUN = true;
		params.GRAPH_FILE = prefix + "sampleNetwork.sif";
		params.DATASETS_FILE_HAS_HEADER = true;
		//params.MIN_PERCENTAGE = 5;
		//kpmSettings.MIN_L = new HashMap<String, Integer>()
		kpmSettings.INCLUDE_CHARTS=true;

		//params.GRAPH_FILE = params.GRAPH_FILE; */

		try {
			start(args, prefix + "kpm.properties", params);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	

	public void start(String[] args, String propertiesFile, Parameters params) throws Exception {

		// parse parameters from the command line
	//	try {
		ArgsParametersParser argspp  = new ArgsParametersParser(kpmSettings);
		//As far as I understood, this method creates the param object but also modifies
        // kmpSettings to set up the correct values for the run.
			params = argspp.parse(args, params);
			System.out.println("Arrived here");
	//	} catch (Exception e1) {
	//		System.out.println("An error occurred. Exiting.");
	//		return;
	//	}
		// check if all input files are correct
		InputFileParser ifp = new InputFileParser(kpmSettings);
		params = ifp.parse(params);
	System.out.println(params.PROGRAM);
		if (params.PROGRAM == Program.SP) {
			System.out.println("PROGRAM SELECTED: Shortest Paths");
			ShortestPathAlgorithms.shortestPathways(params.GRAPH_FILE, params);
		} else {

			KPMRunHandler kpmHandler = new KPMRunHandler(kpmSettings);
			kpmHandler.runBatch(params);
			
			/*for(List<Integer> kl_values : KPMParameters.STATS_MAP.keySet()){
				String param = "";
				for(Integer val : kl_values){
					param +=  val + " ";
				}
				System.out.println("param = " + param);
				List<Result> results = KPMParameters.STATS_MAP.get(kl_values).getResults();

				if (results.size() > 1) {
					System.out.println("PERFORMING POST-PROCESSING OPERATIONS...");

					try {
						Collections.sort(results);
					} catch (UnsupportedOperationException e) {
						// do nothing; this list was already sorted
					}

					if (KPMParameters.NUM_SOLUTIONS >= 0) {
						List<Result> filter = new ArrayList<Result>();
						int toAdd = 1;
						for (Result result : results) {
							if (toAdd > KPMParameters.NUM_SOLUTIONS) {
								break;
							} else {
								filter.add(result);
								toAdd++;
							}
						}
						results = filter;
					}

					System.out.println("BEST FITNESS: " + results.get(0).getFitness());
				} else {
					System.out.println("NO PATHWAYS EXIST FOR THE GIVEN PARAMETERS !");
				}
			 */



			/*System.out.println("\n********* GENERATING STATS AND WRITING SOLUTIONS TO FILES *******\n");

				File theDir = new File(Globals.RESULTS_FOLDER);

				// if the directory does not exist, create it
				if (!theDir.exists()) {
					System.out.println("creating directory: " + Globals.RESULTS_FOLDER);
					boolean created = theDir.mkdir();
					if (created) {
						System.out.println("Results folder created");
					}
				}

				StatisticsUtility.writeResultsStats(results, g);
				if (KPMParameters.GENERATE_SUMMARY_FILE) {
					System.out.println("Generating summary file... ");
					StatisticsUtility.writeSummaryFile(Globals.RESULTS_FOLDER + 
							Globals.FILE_SEPARATOR + KPMParameters.SUMMARY_FILE + "-"
							+ Globals.SUFFIX + Globals.FILE_EXTENSION, results);
				}
				if (KPMParameters.GENERATE_PATHWAYS_FILE) {
					System.out.println("Generating pathway file(s)... ");
					if (!Globals.PATHWAYS_IN_SINGLE_FILE) {
						StatisticsUtility.writeIndividualPathwayFiles(Globals.RESULTS_FOLDER,
								results, g);
					} else {
						StatisticsUtility.writePathwaysFile(Globals.RESULTS_FOLDER + 
								Globals.FILE_SEPARATOR + KPMParameters.PATHWAYS_FILE + "-"
								+ Globals.SUFFIX + Globals.FILE_EXTENSION, results, g);
					}
				}
				if (KPMParameters.GENERATE_PATHWAYS_STATS_FILE) {
					System.out.println("Generating pathways stats file... ");
					StatisticsUtility.writePathwaysStatsFile(Globals.RESULTS_FOLDER + 
							Globals.FILE_SEPARATOR + KPMParameters.PATHWAYS_STATS_FILE + "-"
							+ Globals.SUFFIX + Globals.FILE_EXTENSION, results, g);
				}
				if (KPMParameters.GENERATE_GENE_STATS_FILE) {
					System.out.println("Generating gene stats file... ");
					StatisticsUtility.writeGeneStatsFile(Globals.RESULTS_FOLDER + 
							Globals.FILE_SEPARATOR + KPMParameters.GENE_STATS_FILE + "-"
							+ Globals.SUFFIX + Globals.FILE_EXTENSION, results, g);
				}
				if (KPMParameters.GENERATE_DATASETS_STATS_FILE) {
					System.out.println("Generating datasets stats file... ");
					StatisticsUtility.writeDatasetsStats(Globals.RESULTS_FOLDER + 
							Globals.FILE_SEPARATOR + KPMParameters.DATASETS_STATS_FILE + "-"
							+ Globals.SUFFIX + Globals.FILE_EXTENSION, p, g);
				}            

				if (Globals.PROGRAM == Program.KPM_SP) {
					shortestPathways(results, g);
				}*/
		}
	}
}

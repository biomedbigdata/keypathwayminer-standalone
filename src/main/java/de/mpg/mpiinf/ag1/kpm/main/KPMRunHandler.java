package de.mpg.mpiinf.ag1.kpm.main;

import de.mpg.mpiinf.ag1.kpm.parsers.PriorityFileParser;
import de.mpg.mpiinf.ag1.kpm.utils.KPMStandaloneTaskMonitor;
import de.mpg.mpiinf.ag1.kpm.utils.Parser;
import de.mpg.mpiinf.ag1.kpm.utils.StatisticsUtility;
import de.mpg.mpiinf.ag1.kpm.utils.OutputSettings;

import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.charts.IChart;
import dk.sdu.kpm.results.IKPMResultSet;
import dk.sdu.kpm.results.IKPMRunListener;
import dk.sdu.kpm.runners.BatchRunWithPerturbationParameters;
import dk.sdu.kpm.runners.BatchRunWithPerturbationRunner;
import dk.sdu.kpm.runners.BatchRunner;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class that handles the actual KPM runs, and ensures we get the results.
 * Assumes all other parameters has already been set.
 *
 * @author Martin
 */
public class KPMRunHandler implements IKPMRunListener {

    private volatile KPMSettings kpmSettings;

    private volatile Parameters params;

    public KPMRunHandler(KPMSettings settings) {
        this.kpmSettings = settings;
    }

    /**
     * Starts the runs, ensures the settings are correctly set up.
     * Assumes all other parameters has already been set.
     */
    public void runBatch(Parameters params) {
        this.params = params;
        // Register starting time
        long start, end, time;
        kpmSettings.STARTING_TIME = System.nanoTime();
        start = System.currentTimeMillis();

        System.out.println("\n********** CREATING GRAPH **********");
        // Parse the graph and matrix files
        params.PARSER = new Parser(params.GRAPH_FILE, kpmSettings.MATRIX_FILES_MAP, Parser.TAB);

        // Create KPM graph
        params.INPUT_GRAPH = params.PARSER.createGraph(kpmSettings);

        // Print graph statistics
        if (params.PRINT_GRAPH_STATS) {
            System.out.print("--------------\nNETWORK STATS|\n--------------\n");
            StatisticsUtility.printGraphStatistics(params.INPUT_GRAPH);
        }
        // Print dataset statistics
        if (params.PRINT_DATASETS_STATS) {
            System.out.println("--------------\nDATASETS STATS|\n--------------");
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
        try {
            if (params.VALIDATION_FILE != null) {
                params = setValidationFile(params);
            }
        } catch (FileNotFoundException e) {
            System.out.println(String.format("\nValidation file could not be found: '%s'.\n\nKPM will now terminate.", params.VALIDATION_FILE));
            return;
        } catch (IOException e) {
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
                params.STRATEGY = KPMStrategy.INES;
                params.ALGORITHM = KPMAlgorithm.GREEDY;
                break;

            case LCG:
                params.STRATEGY = KPMStrategy.INES;
                params.ALGORITHM = KPMAlgorithm.ACO;
                break;

            case OPTIMAL:
                params.STRATEGY = KPMStrategy.INES;
                params.ALGORITHM = KPMAlgorithm.OPTIMAL;
                break;

            case EXCEPTIONSUMGREEDY:
                params.STRATEGY = KPMStrategy.GLONE;
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
        if (!configurator.areSettingsValid()) {
            System.out.println("Settings are invalid. Stopping.");
            return;
        }

        System.out.println("\n*************** RUNNING " + params.STRATEGY + " (" + params.ALGORITHM + ") " + " ***************");
        if (params.IS_PERTURBATION_RUN) {
            runBatchWithPerturbation(params);
        } else {
            runStandard(params);
        }

        end = System.currentTimeMillis();
        time = (end - start) / 1000;
        // Do stats and post-processing stuff after obtaining the results
        System.out.println("\n Finished with a total running time of: " + "\033[0;1m" + time + " seconds");

        kpmSettings.TOTAL_RUNNING_TIME = time;
    }

    private void runBatchWithPerturbation(Parameters params) {
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
        } catch (Exception e) {
            Logger.getLogger(KPMRunHandler.class.getName()).log(Level.SEVERE, null, e);
            System.out.println("An error occurred during the KPM run with perturbed graphs. Exiting.");
        }
    }

    private void runStandard(Parameters params) {
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

    @Override
    public void runFinished(IKPMResultSet results) {
        System.out.println("\n********** SAVING RESULTS **********");
        //Check if path for the RESULT_FOLDER was provided
        if (params.RESULTS_FOLDER.isEmpty()) {
            System.out.println("-> Name for the results folder is empty.\nSaving results process is aborted");
            return;
        }
        //Create results folder if folder does not already exist
        try {
            Files.createDirectories(Paths.get(params.RESULTS_FOLDER));
        } catch (IOException e) {
            e.printStackTrace();
        }

        String time = new SimpleDateFormat("yyyy-MM-dd_HH.mm.ss").format(Calendar.getInstance().getTime());

        Path currentResultsFolder = Paths.get(params.RESULTS_FOLDER + File.separator + time + "_" + params.STRATEGY + "_" + results.getKpmSettings().ALGO);

        //Create folder to save the tables and charts for the current run
        try {
            Files.createDirectories(currentResultsFolder);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        System.out.println("Results will be saved in: " + currentResultsFolder.toString());
        /*
             Saving the charts of the current run if charts were created
         */
        if (!results.getCharts().isEmpty()) {
            System.out.print("--------------\nSAVING CHARTS|\n--------------\n");
            Path resultsChartsFolder = Paths.get(currentResultsFolder + File.separator + "charts");
            //Create folder to save the charts for the current run
            try {
                Files.createDirectories(resultsChartsFolder);
            } catch (IOException ioe) {
                ioe.printStackTrace();
            }
            //Saving all charts
            int chartNumber = 0;
            for (IChart chart : results.getCharts().values()) {
                String chartFileName = resultsChartsFolder.toString() + File.separator + String.format("chart%d.png", chartNumber);
                //System.out.println(chartFileName);
                File chartsFile = new File(chartFileName);
                try {
                    chart.saveAsPng(chartsFile, 1600, 1200);
                } catch (IOException e) {
                    e.printStackTrace();
                }
                chartNumber++;
            }
            if (chartNumber > 0) {
                System.out.println("-> Saved " + chartNumber + " charts to: '" + resultsChartsFolder.toString() + "'");
            }
        }

        /*
            Saving result tables to files
         */
        System.out.print("--------------\nSAVING TABLES|\n--------------\n");

        Path resultsTableFolder = Paths.get(currentResultsFolder + File.separator + "tables");

        //Create folder to store the tables of the current run
        try {
            Files.createDirectories(resultsTableFolder);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        // Start writing tables
        if (OutputSettings.GENERATE_SUMMARY_FILE) {
            String resultSummary = resultsTableFolder.toString() + File.separator + OutputSettings.SUMMARY_FILE + params.FILE_EXTENSION;
            StatisticsUtility.writeSummaryFile(resultSummary, results, kpmSettings, params);
        }
        // Start writing dataset statistics files for every dataset
        if (OutputSettings.GENERATE_DATASETS_STATS_FILE) {
            String resultDatasetStats = resultsTableFolder.toString() + File.separator + OutputSettings.DATASETS_STATS_FILE + params.FILE_EXTENSION;
            StatisticsUtility.writeDatasetsStats(resultDatasetStats, params.PARSER, kpmSettings.MAIN_GRAPH);
        }

        // Start writing gene statistics files for every dataset
        if (OutputSettings.GENERATE_GENE_STATS_FILE) {
            String resultGeneStats = resultsTableFolder.toString() + File.separator + OutputSettings.GENE_STATS_FILE;
            StatisticsUtility.writeGeneStatsFile(resultGeneStats, results, kpmSettings.MAIN_GRAPH, kpmSettings);
        }

        // Start writing pathways file which contains node and interactions in one file
        if (OutputSettings.GENERATE_PATHWAYS_FILE) {
            String resultIndividualPathwaysFile = resultsTableFolder.toString();
            StatisticsUtility.writeIndividualPathwayFiles(resultIndividualPathwaysFile, results, kpmSettings.MAIN_GRAPH, params);
            String resultPathwaysFile = resultsTableFolder.toString() + File.separator + OutputSettings.PATHWAYS_FILE;
            StatisticsUtility.writePathwaysFile(resultPathwaysFile, results, kpmSettings.MAIN_GRAPH, kpmSettings);
        }

        //Start writing pathway statistics file
        if (OutputSettings.GENERATE_PATHWAYS_STATS_FILE) {
            String resultPathwayStatsFile = resultsTableFolder.toString() + File.separator + OutputSettings.PATHWAYS_STATS_FILE + params.FILE_EXTENSION;
            StatisticsUtility.writePathwaysStatsFile(resultPathwayStatsFile, results, kpmSettings.MAIN_GRAPH);
        }
        //Start writing general statistics file
        if (OutputSettings.GENERATE_GENERAL_STATS_FILE) {
            StatisticsUtility.writeResultsStats(resultsTableFolder.toString(), results, kpmSettings.MAIN_GRAPH, params);
        }

        System.out.println("-> Saved tables to: " + resultsTableFolder.toString());

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


}

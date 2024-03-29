package de.mpg.mpiinf.ag1.kpm.parsers;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.mpg.mpiinf.ag1.kpm.main.Parameters;
import de.mpg.mpiinf.ag1.kpm.utils.Program;
import de.mpg.mpiinf.ag1.kpm.main.Main;
import de.mpg.mpiinf.ag1.kpm.output.OutputSettings;
import de.mpg.mpiinf.ag1.kpm.utils.Separator;
import dk.sdu.kpm.Algo;
import dk.sdu.kpm.Combine;
import dk.sdu.kpm.Heuristic;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.algo.glone.LocalSearch;
import dk.sdu.kpm.algo.glone.RhoDecay;
import dk.sdu.kpm.perturbation.PerturbationService;
import dk.sdu.kpm.perturbation.IPerturbation.PerturbationTags;

public class ArgsParametersParser {
    private Parameters parameters;
    private volatile KPMSettings kpmSettings;

    //List of all algorithms to test if the right algorithm was provided
    private final List<String> algoList = Arrays.asList("GREEDY", "ACO", "OPTIMAL");
    private String algorithm;

    //List of all strategies to test if the right strategy was provided
    private final List<String> strategyList = Arrays.asList("INES", "GLONE");
    private String strategy;

    private HashMap<String, String> id2path = new HashMap<String, String>();
    private HashMap<String, Integer> id2param = new HashMap<String, Integer>();


    public ArgsParametersParser(KPMSettings settings, Parameters params) {
        this.kpmSettings = settings;
        this.parameters = params;
        //Set default algorithm and strategy before parsing the arguments
        algorithm = parameters.ALGORITHM.toString();
        strategy = parameters.STRATEGY.toString();
    }

    /**
     * Parses command line arguments and
     * saves them to the corresponding objects.
     *
     * @param args Arguments passed from the commandline.
     * @return Parameters object.
     * @throws Exception when input is not valid
     */
    public Parameters parse(String[] args) throws Exception {
        //Read command line arguments
        parseParameters(args, parameters);
        // Number of processors used
        System.out.println("Processors used: " + kpmSettings.NUMBER_OF_PROCESSORS+"/"+Runtime.getRuntime().availableProcessors());

        if (parameters.PROGRAM == Program.SP) {
            // do nothing
        } else if (!id2path.isEmpty()) {
            //It is important to print the following information here so that the user knows if he named a parameter wrong
            for (String id : id2path.keySet()) {
                if (!(new File(id2path.get(id))).isFile()) {
                    System.err.println(id2path.get(id) + " not found !");
                    System.exit(-1);
                }
            }
            kpmSettings.MATRIX_FILES_MAP = id2path;

            //Is this if loop necessary????
            if (kpmSettings.IS_BATCH_RUN) {
                List<String> invalids = new ArrayList<String>();
                for (String id : id2param.keySet()) {
                    if (id.startsWith("L")) {
                        invalids.add(id);
                    }
                }
                for (String invalidID : invalids) {
                    id2param.remove(invalidID);
                }
            }

            kpmSettings.CASE_EXCEPTIONS_MAP = id2param;
        } else {
            //If a dataset file exists read input from dataset file
            if (parameters.DATASETS_FILE == null) {
                System.out.println("No datasets file was specified.");
                System.exit(-1);
            }

            if (!new File(parameters.DATASETS_FILE).isFile()) {
                System.out.println("Datasets file " + parameters.DATASETS_FILE + " does not exist.");
                System.exit(-1);
            }
            //If a dataset file was provide take input arguments from there
            DatasetsFileParser dfp = new DatasetsFileParser(kpmSettings);
            parameters = dfp.parse(parameters.DATASETS_FILE_SEPARATOR.charValue(), parameters.DATASETS_FILE_HAS_HEADER, parameters);
        }
            /*
                In Batch run:
                Check if some parameters were entered without a batch flag.
                If this is the case set MIN=MAX and INC=0.
            */
        if (kpmSettings.IS_BATCH_RUN) {
            System.out.println("BATCH run will be performed with following parameters:");
            //Case exceptions
            for (String identifier : id2path.keySet()) {
                if (kpmSettings.MAX_L.get(identifier) == null || kpmSettings.INC_L.get(identifier) == null) {
                    kpmSettings.MAX_L.put(identifier, kpmSettings.MIN_L.get(identifier));
                    kpmSettings.INC_L.put(identifier, 0);
                }
            }
        }

        // Print parameter and dataset information
        printInformation();
        //Run checks on identifiers of the data sets and the case-exception parameters
        runChecks();

        //Set the right algorithm depending on the strategy
        setAlgorithm();


        return parameters;
    }


    /*
        Method which goes through the command line arguments
     */
    private void parseParameters(String[] args, Parameters params) throws Exception {
        for (String arg : args) {
            // Split argument on equals symbol
            String[] options = arg.split("=");

            // Input files
            if (options[0].equals("-graphFile")) {
                params.GRAPH_FILE = options[1];
            } else if (options[0].startsWith("-isDirected")) {
                params.IS_GRAPH_DIRECTED = Boolean.parseBoolean(options[1]);
            } else if (options[0].startsWith("-matrix")) {
                String id = "L" + options[0].substring(7);
                String path = options[1];
                String internalID = kpmSettings.externalToInternalIDManager.getOrCreateInternalIdentifier(id);

                id2path.put(internalID, path);
            } else if (options[0].equals("-gfHeader")) {
                params.GRAPH_FILE_HAS_HEADER = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gfSep")) {
                params.GRAPH_FILE_SEPARATOR = Separator.valueOf(options[1]);
            } else if (options[0].equals("-datasetsFile")) {
                params.DATASETS_FILE = options[1];
            } else if (options[0].equals("-dfHeader")) {
                params.DATASETS_FILE_HAS_HEADER = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-dfSep")) {
                params.DATASETS_FILE_SEPARATOR = Separator.valueOf(options[1]);
            } else if (options[0].equals("-positiveFile")) {
                params.POSITIVE_FILE = options[1];
            } else if (options[0].equals("-negativeFile")) {
                params.NEGATIVE_FILE = options[1];
            } else if (options[0].equals("-mfHeader")) {
                params.MATRIX_FILES_HAVE_HEADER = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-mfSep")) {
                params.MATRIX_FILES_SEPARATOR = Separator.valueOf(options[1]);
                // OUTPUT FILES
            } else if (options[0].equals("-fileExt")) {
                params.FILE_EXTENSION = options[1];
            } else if (options[0].equals("-summaryFile")) {
                OutputSettings.SUMMARY_FILE = options[1];
            } else if (options[0].equals("-pathwaysFile")) {
                OutputSettings.PATHWAYS_FILE = options[1];
            } else if (options[0].equals("-pathwaysStatsFile")) {
                OutputSettings.PATHWAYS_STATS_FILE = options[1];
            } else if (options[0].equals("-geneStatsFile")) {
                OutputSettings.GENE_STATS_FILE = options[1];
            } else if (options[0].equals("statsFile")) {
                OutputSettings.GENERAL_STATS_FILE = options[1];
            } else if (options[0].equals("-dataStatsFile")) {
                OutputSettings.DATASETS_STATS_FILE = options[1];
            } else if (options[0].equals("-resultsDir")) {
                params.RESULTS_FOLDER = options[1];
            } else if (options[0].equals("-spStatsFile")) {
                params.SHORTEST_PATHS_STATS_FILE = options[1];
            } else if (options[0].equals("-spFile")) {
                params.SHORTEST_PATH_FILE = options[1];
            } else if (options[0].equals("-spNodeStatsFile")) {
                params.SHORTEST_PATHS_NODE_STATS_FILE = options[1];
            } else if (options[0].equals("-spEdgeStatsFile")) {
                params.SHORTEST_PATHS_EDGE_STATS_FILE = options[1];
            } else if (options[0].equals("-pSingleFile")) {
                params.PATHWAYS_IN_SINGLE_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSummary")) {
                OutputSettings.GENERATE_SUMMARY_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gPathways")) {
                OutputSettings.GENERATE_PATHWAYS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gPathwayStats")) {
                OutputSettings.GENERATE_PATHWAYS_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gGeneStats")) {
                OutputSettings.GENERATE_GENE_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gDataStats")) {
                OutputSettings.GENERATE_DATASETS_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPStats")) {
                params.GENERATE_SHORTEST_PATHS_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPFiles")) {
                params.GENERATE_SHORTEST_PATH_FILES = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPNodes")) {
                params.GENERATE_SHORTEST_PATHS_NODE_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-gSPEdges")) {
                params.GENERATE_SHORTEST_PATHS_EDGE_STATS_FILE = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("runID")) {
                params.RUN_ID = options[1];
                // OUTPUT TO TERMINAL
            } else if (options[0].equals("-pGraphStats")) {
                params.PRINT_GRAPH_STATS = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-pDataStats")) {
                params.PRINT_DATASETS_STATS = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-suffix")) {
                params.SUFFIX = options[1];
                // BASIC params
            } else if (options[0].equals("-program")) {
                params.PROGRAM = Program.valueOf(options[1]);
            } else if (options[0].equals("-sp")) {
                params.SHORTEST_PATHS_REPORTED = Integer.valueOf(options[1]);
            } else if (options[0].equals("-spPathways")) {
                params.PATHWAYS_SHORTEST_PATHS = Integer.valueOf(options[1]);
            } else if (options[0].equals("-spMinLength")) {
                params.MINIMUM_LENGTH_SHORTEST_PATHS = Integer.valueOf(options[1]);
            } else if (options[0].equals("-K")) {
                kpmSettings.GENE_EXCEPTIONS = Integer.parseInt(options[1]);
                // In case the K parameter without a batch flag was used in a batch run
                kpmSettings.MIN_K = kpmSettings.GENE_EXCEPTIONS;
                kpmSettings.INC_K = 0;
                kpmSettings.MAX_K = kpmSettings.GENE_EXCEPTIONS;
            } else if (options[0].equals("-K_batch")) {
                String[] values = options[1].split(",");

                if (values.length != 3) {
                    throw new Exception("Invalid settings for " + options[0]);
                }

                kpmSettings.MIN_K = Integer.parseInt(values[0]);
                kpmSettings.INC_K = Integer.parseInt(values[1]);
                kpmSettings.MAX_K = Integer.parseInt(values[2]);

            } else if (options[0].startsWith("-L") && !options[0].contains("batch")) {
                String id = "L" + options[0].substring(2);
                int l = Integer.parseInt(options[1]);

                String internalID = kpmSettings.externalToInternalIDManager.getOrCreateInternalIdentifier(id);

                kpmSettings.MIN_L.put(internalID, l);
                kpmSettings.VARYING_L_ID.add(internalID);

                id2param.put(internalID, l);
            } else if (options[0].matches("-L[1-9][0-9]*[_]batch")) {
                // If batch option is set assign ranged L value to n th matrix
                String id = options[0].substring(1, options[0].indexOf('_'));

                String internalID = kpmSettings.externalToInternalIDManager.getOrCreateInternalIdentifier(id);

                String[] values = options[1].split(",");

                if (values.length != 3) {
                    throw new Exception("Invalid settings for " + options[0]);
                }
                kpmSettings.MIN_L.put(internalID, Integer.parseInt(values[0]));
                kpmSettings.INC_L.put(internalID, Integer.parseInt(values[1]));
                kpmSettings.MAX_L.put(internalID, Integer.parseInt(values[2]));
                kpmSettings.VARYING_L_ID.add(internalID);
            } else if (options[0].equals("-algo")) {
                if (!algoList.contains(options[1])) {
                    System.err.println(options[1]
                            + " is not a valid algorithm. " + ("Valid options are: GREEDY, ACO or OPTIMAL"));
                    System.exit(-1);
                } else {
                    algorithm = options[1];
                }

            } else if (options[0].equals("-strategy")) {
                if (!strategyList.contains(options[1])) {
                    System.err.println(options[1]
                            + " is not a valid strategy. " + ("Valid options are: INES or GLONE"));
                    System.exit(-1);

                } else {
                    strategy = options[1];
                }

                // ADVANCED PARAMETERS (GENERAL)

            } else if (options[0].equals("-numProc")) {
                int input = Integer.parseInt(options[1]);
                int available = Runtime.getRuntime().availableProcessors();
                if (input < 1) {
                    kpmSettings.NUMBER_OF_PROCESSORS = 1;
                } else if (input > available) {
                    kpmSettings.NUMBER_OF_PROCESSORS = available;
                } else {
                    kpmSettings.NUMBER_OF_PROCESSORS = input;
                }
            } else if (options[0].equals("-nodeHeuristic")) {
                kpmSettings.NODE_HEURISTIC_VALUE = Heuristic.valueOf(options[1]);
            } else if (options[0].equals("-combineOp")) {
                kpmSettings.COMBINE_OPERATOR = Combine.valueOf(options[1]);
            } else if (options[0].equals("-combineFormula")) {
                kpmSettings.COMBINE_FORMULA = options[1].replace('~', '!');
            } else if (options[0].equals("-eval")) {
                kpmSettings.EVAL = Boolean.parseBoolean(options[1]);
            } else if (options[0].equals("-maxSolutions")) {
                kpmSettings.NUM_SOLUTIONS = Integer.parseInt(options[1]);
            } else if (options[0].equals("-doubleSolutions")) {
                kpmSettings.DOUBLE_SOLUTIONS_ALLOWED = Boolean.parseBoolean(options[1]);

                // ADVANCED PARAMETERS (ACO ONLY)
            } else if (options[0].equals("-alpha")) {
                kpmSettings.ALPHA = Double.parseDouble(options[1]);
            } else if (options[0].equals("-beta")) {
                kpmSettings.BETA = Double.parseDouble(options[1]);
            } else if (options[0].equals("-rho")) {
                kpmSettings.RHO = Double.parseDouble(options[1]);
            } else if (options[0].equals("-iterations")) {
                if (Integer.parseInt(options[1]) == 0) {
                    kpmSettings.MAX_ITERATIONS = Integer.MAX_VALUE;
                } else {
                    kpmSettings.MAX_ITERATIONS = Integer.parseInt(options[1]);
                }
            } else if (options[0].equals("-seed")) {
                kpmSettings.SEED = Long.parseLong(options[1]);
                kpmSettings.R = new Random(kpmSettings.SEED);
            } else if (options[0].equals("-tradeoff")
                    && options[1].equals("multiplicative")) {
                kpmSettings.MULTIPLICATIVE_TRADEOFF = true;
            } else if (options[0].equals("-tradeoff")
                    && options[1].equals("additive")) {
                kpmSettings.MULTIPLICATIVE_TRADEOFF = false;
            } else if (options[0].equals("-solutions")) {
                kpmSettings.NUMBER_OF_SOLUTIONS_PER_ITERATION = Integer.parseInt(options[1]);
            } else if (options[0].equals("-startNodes")) {
                kpmSettings.NUM_STARTNODES = Integer.parseInt(options[1]);
            } else if (options[0].equals("-rhoDecay")) {
                try {
                    kpmSettings.RHO_DECAY = RhoDecay.valueOf(options[1]);
                } catch (IllegalArgumentException e) {
                    System.err.println(options[1]
                            + " is not a valid rho decay function. "
                            + Algo.values());
                    System.exit(-1);
                }
            } else if (options[0].equals("-iterationBased")) {
                try {
                    kpmSettings.ITERATION_BASED = Boolean.parseBoolean(options[1]);
                } catch (IllegalArgumentException e) {
                    System.err.println(options[1]
                            + " is not a valid local search method. "
                            + LocalSearch.values());
                    System.exit(-1);
                }
            } else if (options[0].equals("-localSearch")) {
                kpmSettings.L_SEARCH = LocalSearch.valueOf(options[1]);
            } else if (options[0].equals("-maxRunsWithoutChange")) {
                if (Integer.parseInt(options[1]) == 0) {
                    kpmSettings.MAX_RUNS_WITHOUT_CHANGE = Integer.MAX_VALUE;
                } else {
                    kpmSettings.MAX_RUNS_WITHOUT_CHANGE = Integer.parseInt(options[1]);
                }
            } else if (options[0].equals("-tauMin")) {
                kpmSettings.TAU_MIN = Double.parseDouble(options[1]);
            } else if (options[0].equals("-batch")) {
                kpmSettings.IS_BATCH_RUN = true;
            } else if (options[0].equals("-removeBens")) {
                kpmSettings.REMOVE_BENs = true;
            } else if (options[0].equals("-validationFile")) {
                params.VALIDATION_FILE = options[1];
            } else if (options[0].equals("-perturbation")) {
                params.IS_PERTURBATION_RUN = true;
                String[] values = options[1].split(",");

                if (values.length != 4) {
                    throw new Exception("Invalid settings for " + options[0]);
                }
                params.MIN_PERCENTAGE = Integer.parseInt(values[0]);
                params.STEP_PERCENTAGE = Integer.parseInt(values[1]);
                params.MAX_PERCENTAGE = Integer.parseInt(values[2]);
                params.GRAPHS_PER_STEP = Integer.parseInt(values[3]);
            } else if (options[0].equals("-perturbationTechnique")) {
                if (options[1].equals("edgeremove")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.EdgeRemoval);

                } else if (options[1].equals("edgerewire")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.EdgeRewire);

                } else if (options[1].equals("nodeswap")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.NodeSwap);

                } else if (options[1].equals("noderemove")) {
                    params.PERTURBATION = PerturbationService.getPerturbation(PerturbationTags.NodeRemoval);

                }

                if (params.PERTURBATION == null) {
                    throw new Exception("Invalid perturbation technique.");
                }

                System.out.println("Perturbation technique: " + params.PERTURBATION.getName());
            } else if (options[0].equals("-processors")) {
                // Number of the processors to be used in the current run
                // Maximal number of available processors
                int max = Runtime.getRuntime().availableProcessors();
                if (options[1].equals("MAX")) {
                    kpmSettings.NUMBER_OF_PROCESSORS = max;
                } else {
                    int val = Integer.parseInt(options[1]);
                    if (val > max) {
                        kpmSettings.NUMBER_OF_PROCESSORS = max;
                    } else kpmSettings.NUMBER_OF_PROCESSORS = Math.max(val, 1);
                }
            } else if (options[0].equals("-help")) {
                printHelp();
                System.exit(0);
            } else {
                System.out.println("Unknown command line argument: "
                        + options[0]);
                System.out.println("Please type \"-help\" for a list of commands and usage.");
                System.exit(-1);
            }
        }
    }

    private void printInformation() {
        System.out.println("Algorithm:\t" + algorithm);
        System.out.println("Strategy:\t" + strategy);
        System.out.println("Graph file:\t" + parameters.GRAPH_FILE);
        //Print gene exceptions
        if (strategy.equals("INES")) {
            if (kpmSettings.IS_BATCH_RUN) {
                if (kpmSettings.MIN_K == kpmSettings.MAX_K && kpmSettings.INC_K == 0) {
                    System.out.println("K:\t" + kpmSettings.MIN_K);
                } else {
                    System.out.println("K:\t" + "[Min: " + kpmSettings.MIN_K + ", Max: "
                            + kpmSettings.MAX_K + ", Step: " + kpmSettings.INC_K + "]");
                }
            } else {
                System.out.println("K:\t" + kpmSettings.GENE_EXCEPTIONS);
            }
        }

        //Print dataset paths corresponding to the identifiers of the id2path map
        for (String id : id2path.keySet()) {
            System.out.println(kpmSettings.externalToInternalIDManager.getExternalIdentifier(id) + ": " + id2path.get(id));
        }

        // Print all case exception(L) parameters and the corresponding identifier
        for (String identifier : kpmSettings.MIN_L.keySet()) {
            if (kpmSettings.IS_BATCH_RUN) {
                if (kpmSettings.MIN_L.get(identifier).equals(kpmSettings.MAX_L.get(identifier)) && kpmSettings.INC_L.get(identifier).equals(0)) {
                    System.out.println(kpmSettings.externalToInternalIDManager.getExternalIdentifier(identifier) + ":\t" + kpmSettings.MIN_L.get(identifier));
                } else {
                    System.out.println(kpmSettings.externalToInternalIDManager.getExternalIdentifier(identifier) + ":\t" + "[Min: " + kpmSettings.MIN_L.get(identifier) + ", Max: "
                            + kpmSettings.MAX_L.get(identifier) + ", Step: " + kpmSettings.INC_L.get(identifier) + "]");
                }
            } else {
                System.out.println(kpmSettings.externalToInternalIDManager.getExternalIdentifier(identifier) + ":\t" + kpmSettings.MIN_L.get(identifier));
            }
        }

        // Print information on how the datasets will be combined if more than one dataset is used
        if (kpmSettings.MATRIX_FILES_MAP.size() > 1) {
            if (kpmSettings.COMBINE_OPERATOR == Combine.CUSTOM) {
                //TODO  checks for custom formula if. Or just remove try catch in KPM core
                //If everything goes well print out custom formula
                System.out.println("Combine Formula: " + kpmSettings.COMBINE_FORMULA);
            } else {
                //In case OR or AND operator was selected
                StringBuilder formula = new StringBuilder("");

                // Getting a Set of Key-value pairs
                Set<String> keySet = kpmSettings.MATRIX_FILES_MAP.keySet();
                int numOfVars = keySet.size();

                // Determine operator
                String operator = ((kpmSettings.COMBINE_OPERATOR == Combine.OR) ? "||" : "&&");
                for (String indentifier : keySet) {
                    if (numOfVars == 1) {
                        //last identifier
                        formula.append(kpmSettings.externalToInternalIDManager.getExternalIdentifier(indentifier));
                        break;
                    }
                    formula.append(kpmSettings.externalToInternalIDManager.getExternalIdentifier(indentifier)).append(operator);
                    numOfVars--;
                }
                kpmSettings.COMBINE_FORMULA = formula.toString();
                System.out.println("Combine Formula:\t" + formula);
            }
        }
    }

    /*
        Set one of 6 algorithm options in kpmSettings
     */
    private void setAlgorithm() {
        if (strategy.equals("INES")) {
            switch (algorithm) {
                case "GREEDY":
                    kpmSettings.ALGO = Algo.GREEDY;
                    break;
                case "ACO":
                    kpmSettings.ALGO = Algo.LCG;
                    break;
                case "OPTIMAL":
                    kpmSettings.ALGO = Algo.OPTIMAL;
                    break;
            }
        } else if (strategy.equals("GLONE")) {
            switch (algorithm) {
                case "GREEDY":
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMGREEDY;
                    break;
                case "ACO":
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMACO;
                    break;
                case "OPTIMAL":
                    kpmSettings.ALGO = Algo.EXCEPTIONSUMOPTIMAL;
                    break;
            }
        }
    }

    /*
    Runs checks between the identifiers of the identifiers of the matrices
    and the identifiers of the L parameters
    */
    private void runChecks() {
        //Compares the number of matrices against the number of given L parameters
        if (kpmSettings.MIN_L.keySet().size() != kpmSettings.MATRIX_FILES_MAP.size()) {
            System.err.println(String.format(
                    "\nThe were found setup for %d L-parameters, this amount does not match the %d found for the" +
                            " matrix files map.\n->\tKPM will now terminate.",
                    kpmSettings.MIN_L.keySet().size(),
                    kpmSettings.MATRIX_FILES_MAP.size()
            ));
            System.exit(-1);
        }
        //Check if the identifiers for the matrices is the same as the identifier for the Case-Exceptions(L)
        for (String parameterIdentifier : kpmSettings.MIN_L.keySet()) {
            if (!kpmSettings.MATRIX_FILES_MAP.keySet().contains(parameterIdentifier)) {
                System.err.println("No matrix was specified for identifier:\t" + parameterIdentifier + "\nMake sure numerical identifiers for matrices and" +
                        " case-exception parameters are the same\n->\tKPM will now terminate.");
                System.exit(-1);
            }
        }
    }

    //Prints Readme on command line
    private static void printHelp() {
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader("README.md"));
            String line = "";
            while ((line = br.readLine()) != null) {
                System.out.println(line);
            }
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ioe) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ioe);
        } finally {
            try {
                br.close();
            } catch (IOException ex) {
                Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }


}

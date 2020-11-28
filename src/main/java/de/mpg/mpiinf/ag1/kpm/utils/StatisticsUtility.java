package de.mpg.mpiinf.ag1.kpm.utils;

import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import de.mpg.mpiinf.ag1.kpm.main.Parameters;
import de.mpg.mpiinf.ag1.kpm.main.Program;
import dk.sdu.kpm.Algo;
import dk.sdu.kpm.Combine;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.graph.GeneEdge;
import dk.sdu.kpm.graph.GeneNode;
import dk.sdu.kpm.graph.KPMGraph;
import dk.sdu.kpm.graph.Result;
import dk.sdu.kpm.results.IKPMResultItem;
import dk.sdu.kpm.results.IKPMResultSet;

public class StatisticsUtility {

    public static final boolean PRINT_GENE_LIST = false;
    public static boolean GENERATE_STATISTICS = false;
    public static boolean GENERATE_HIST_LVALUES = false;


    //Write summary file with general information on the current run
    public static void writeSummaryFile(String file, IKPMResultSet results, KPMSettings kpmSettings, Parameters parameters) {
        System.out.println(">Writing summary file");
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));

            bw.write("Gene exceptions (K):\t" + kpmSettings.GENE_EXCEPTIONS + "\n");

            for (String expId : kpmSettings.CASE_EXCEPTIONS_MAP.keySet()) {
                bw.write(expId + ":\t[Case exceptions (L): " + kpmSettings.CASE_EXCEPTIONS_MAP.get(expId)
                        + ", Matrix: " + kpmSettings.MATRIX_FILES_MAP.get(expId) + "]" + "\n");
            }
            //In case several data sets were combined print information on how the data sets were combined
            if (kpmSettings.MATRIX_FILES_MAP.size() > 1) {
                bw.write("Combine Formula:\t" + kpmSettings.COMBINE_FORMULA + "\n");
            }
            //Write which strategy and algorithm was used
            if (kpmSettings.ALGO == Algo.GREEDY) {
                bw.write("Strategy:\t" + "INES" + "\n");
                bw.write("Algorithm:\t" + "GREEDY" + "\n");
            } else if (kpmSettings.ALGO == Algo.LCG) {
                bw.write("Strategy:\t" + "INES" + "\n");
                bw.write("Algorithm:\t" + "ACO" + "\n");
            } else if (kpmSettings.ALGO == Algo.OPTIMAL) {
                bw.write("Strategy:\t" + "INES" + "\n");
                bw.write("Algorithm:\t" + "OPTIMAL" + "\n");
            } else if (kpmSettings.ALGO == Algo.EXCEPTIONSUMGREEDY) {
                bw.write("Strategy:\t" + "GLONE" + "\n");
                bw.write("Algorithm:\t" + "GREEDY" + "\n");
            } else if (kpmSettings.ALGO == Algo.EXCEPTIONSUMACO) {
                bw.write("Strategy:\t" + "GLONE" + "\n");
                bw.write("Algorithm:\t" + "ACO" + "\n");
            } else if (kpmSettings.ALGO == Algo.EXCEPTIONSUMOPTIMAL) {
                bw.write("Strategy:\t" + "GLONE" + "\n");
                bw.write("Algorithm:\t" + "OPTIMAL" + "\n");
            }
            //Program that was used
            if (parameters.PROGRAM == Program.KPM) {
                bw.write("Program:\t" + "KeyPathwayMiner" + "\n");
            } else if (parameters.PROGRAM == Program.SP) {
                bw.write("Program:\t" + "Shortest paths" + "\n");
            } else if (parameters.PROGRAM == Program.KPM_SP) {
                bw.write("Program:\t" + "KeyPathwayMiner followed by SP of the resulting pathways" + "\n");
            }

            //Perturbation information
            if (parameters.IS_PERTURBATION_RUN) {
                bw.write("Perturbation technique:\t" + parameters.PERTURBATION.getName() + "\n");
                bw.write("Perturbation parameters:\t" +
                        "[Start percent:" + parameters.MIN_PERCENTAGE +
                        ", Step percent:" + parameters.STEP_PERCENTAGE +
                        ", Max percent:" + parameters.MAX_PERCENTAGE +
                        ", Graphs per step:" + parameters.GRAPHS_PER_STEP + "]\n");
            }

            //TODO: check whether this is actually the best result aka. if resultList is ordered
            IKPMResultItem best = results.getResults().get(0);
            bw.write("Size of the best result:\t" + best.getResultsInfoTable()[0][2] + "\n");
            bw.write("Total running time:\t" + kpmSettings.TOTAL_RUNNING_TIME + "\n");
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void writeDatasetsStats(String file, Parser p, KPMGraph g) {
        System.out.println(">Writing dataset statistics files");

        BufferedWriter bw = null;
        int totalNodes = g.getVertexCount();
        try {
            bw = new BufferedWriter(new FileWriter(file));
            Map<String, Integer> numCasesMap = p.getNumCasesMap();
            Map<String, Integer> numGenesMap = p.getNumGenesMap();
            Map<String, Set<String>> backNodesMap = p.getBackNodesByExpMap();
            Map<String, Set<String>> backGenesMap = p.getBackGenesMap();
            Map<String, Double> avgCasesMap = p.getAvgExpressedCasesMap();
            Map<String, Double> avgGenesMap = p.getAvgExpressedGenesMap();

            bw.write("MATRIX_ID" + "\t");
            bw.write("CASES" + "\t");
            bw.write("GENES" + "\t");
            bw.write("AVG_EXP_CASES_p_GENE" + "\t");
            bw.write("AVG_EXP_GENES_p_CASE" + "\t");
            bw.write("NODES_WITH_EXPRESSION_MAPPINGS" + "\n");
//          bw.write("GENES_NOT_MAPPED" + "\n");

            DecimalFormat df = new DecimalFormat("#.00");

            for (String expId : numCasesMap.keySet()) {
                bw.write(expId + "\t");
                bw.write(String.valueOf(numCasesMap.get(expId)) + "\t");
                bw.write(String.valueOf(numGenesMap.get(expId)) + "\t");
                bw.write(String.valueOf(avgCasesMap.get(expId)) + "\t");
                bw.write(String.valueOf(avgGenesMap.get(expId)) + "\t");
                int backNodes = backNodesMap.get(expId).size();
                int mappedNodes = totalNodes - backNodes;
                double mappedNodesP = ((double) mappedNodes / (double) totalNodes) * 100.0;
                String mappedNodesPS = df.format(mappedNodesP);
                bw.write(String.valueOf(mappedNodes) + "(" + mappedNodesPS + "%)" + "\n");

                //bw.write(String.valueOf(backGenesMap.get(expId).size()) + "\n");
                System.out.println("\t>Dataset statistics file was written for dataset:\t" + expId);
            }
            bw.close();
        } catch (IOException ex) {
            System.out.println("ERROR IN WRITING DATASETS FILE");
            Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public static void writeGeneStatsFile(String file, IKPMResultSet results, KPMGraph g, KPMSettings kpmSettings) {
        System.out.println(">Writing gene statistics files");

        DecimalFormat df = new DecimalFormat("#.00");
        boolean isInes = kpmSettings.ALGO == Algo.GREEDY || kpmSettings.ALGO == Algo.LCG
                || kpmSettings.ALGO == Algo.OPTIMAL;
        BufferedWriter bw = null;

        try {
            //legacy code
            /*Map<String, Integer> nodeCount = new HashMap<String, Integer>();
            for (IKPMResultItem result : results.getResults()) {
                //TODO:  unionNode set
                for (String nodeId : result.getUnionNodeSet().keySet()) {
                    if (!nodeCount.containsKey(nodeId)) {
                        nodeCount.put(nodeId, 1);
                    } else {
                        nodeCount.put(nodeId, nodeCount.get(nodeId) + 1);
                    }
                }
            }*/

            for (IKPMResultItem result : results.getResults()) {
                // TODO variable file extension
                System.out.println("\t>Writing gene statistics file for: " + "Pathway-k-" + result.getK() + "-l-" + result.getL());
                bw = new BufferedWriter(new FileWriter(file + "_" + "PATHWAY-k-" + result.getK() + "-l-" + result.getL() + ".txt"));
                bw.write("ID" + "\t");
                bw.write("# PATHWAYS CONTAINING" + "\t");
                bw.write("AVG. DIFF. EXP. CASES" + "\t");
                bw.write("DEGREE");
                if (isInes) {
                    bw.write("\t" + "IS EXCEPTION?");
                }
                bw.write("\n");
                for (GeneNode node : g.getVertices()) {
                    String nodeId = node.getNodeId();
                    bw.write(nodeId + "\t");
                    if (result.getUnionNodeSetCounts().containsKey(nodeId)) {
                        bw.write(result.getUnionNodeSetCounts().get(nodeId) + "\t");
                    } else {
                        bw.write("0" + "\t");
                    }
                    if (node.getAverageExpressedCases() != 0) {
                        bw.write(df.format(node.getAverageExpressedCases()) + "\t");
                    } else {
                        bw.write(0 + "\t");
                    }
                    bw.write(String.valueOf(g.degree(node)));
                    if (isInes) {
                        bw.write("\t" + String.valueOf(!node.isValid()));
                    }
                    bw.write("\n");
                }
            }
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    //TODO unite writeIndividualPathwayFile and writePathwaysFile?
    public static void writeIndividualPathwayFiles(String dir, IKPMResultSet results, KPMGraph g, Parameters params) {
        System.out.println(">Writing individual pathways files");
        for (IKPMResultItem result : results.getResults()) { // list: order consistent
            System.out.println("\t>Writing individual pathways file for: " + "Pathway-k-" + result.getK() + "-l-" + result.getL());
            String nodeFile = dir + params.FILE_SEPARATOR + "pathway-k-" + result.getK() + "-l-" + result.getL()
                    + "-nodes-" + params.SUFFIX + params.FILE_EXTENSION;
            nodeFile = checkFileName(nodeFile);
            String edgeFile = dir + params.FILE_SEPARATOR + "pathway-k-" + result.getK() + "-l-" + result.getL()
                    + "-interactions-" + params.SUFFIX + params.FILE_EXTENSION;
            edgeFile = checkFileName(edgeFile);
            int i = 0;
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(nodeFile))) {
                try (BufferedWriter br = new BufferedWriter(new FileWriter(edgeFile))) {
                    for (Map<String, GeneNode> nodes : result.getAllComputedNodeSets()) { // list: order consistent
                        List<String[]> edges = g.getEdgesConnecting(nodes.values());
                        for (String nodeId : nodes.keySet()) { // set: order not consistent, but it does not matter
                            bw.write(i + "\t");
                            bw.write(nodeId);
                            bw.newLine();

                        }
                        for (String[] edge : edges) {
                            br.write(i + "\t");
                            br.write(edge[0] + "\t" + edge[1]);
                            br.newLine();
                        }
                        i = i + 1;
                    }
                }
            } catch (IOException ioe) {
                Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ioe);
            }

        }
    }

    public static void writePathwaysFile(String file, IKPMResultSet results, KPMGraph g, KPMSettings kpmSettings) {
        System.out.println(">Writing pathways files");
        int size = results.getResults().size();
        int zeros = (int) Math.floor(Math.log10(size)) + 1;
        String f = "";
        for (int i = 1; i <= zeros; i++) {
            f += "0";
        }
        DecimalFormat df = new DecimalFormat(f);

        boolean isInes = kpmSettings.ALGO == Algo.GREEDY || kpmSettings.ALGO == Algo.LCG
                || kpmSettings.ALGO == Algo.OPTIMAL;


        for (IKPMResultItem result : results.getResults()) {
            int id = 1;
            try (BufferedWriter bw = new BufferedWriter(new FileWriter(file + "-k-" + result.getK() + "-l-" + result.getL() + ".txt"))) {
                System.out.println("\t>Writing pathways file for: " + "Pathway-k-" + result.getK() + "-l-" + result.getL());
                for (Map<String, GeneNode> nodes : result.getAllComputedNodeSets()) {
                    bw.write(id + "\n");
                    bw.write("NODES" + "\t" + "IS EXCEPTION?" + "\n");
                    if (isInes) {
                        for (String node : nodes.keySet()) {
                            bw.write(node + "\t" + nodes.get(node).isValid() + "\n");
                        }
                    } else {
                        for (String node : nodes.keySet()) {
                            bw.write(node + "\n");
                        }
                    }
                    List<String[]> edges = g.getEdgesConnecting(nodes.values());
                    bw.write("EDGES" + "\n");
                    for (String[] edge : edges) {
                        bw.write(edge[0] + "\t" + "pp" + "\t" + edge[1] + "\n");
                    }
                    id++;
                }
            } catch (IOException ex) {
                Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }

    public static void writeResultsStats(String tableFolder, IKPMResultSet results, KPMGraph g, Parameters params) {
        System.out.println(">Writing node and edge pathway hits files");

        Map<String, Integer> nodeHits = new HashMap<String, Integer>();
        Map<String, Integer> edgeHits = new HashMap<String, Integer>();
        String nodeFileName = tableFolder + params.FILE_SEPARATOR + "node_pathway_hits-" + params.SUFFIX + params.FILE_EXTENSION;
        String edgeFileName = tableFolder + params.FILE_SEPARATOR + "edge_pathway_hits-" + params.SUFFIX + params.FILE_EXTENSION;
        nodeFileName = checkFileName(nodeFileName);
        edgeFileName = checkFileName(edgeFileName);

        int maxNode = 0;
        int maxEdge = 0;
        for (GeneNode node : g.getVertices()) {
            String nodeId = node.getNodeId();
            nodeHits.put(nodeId, 0);
        }

        for (GeneEdge edge : g.getEdges()) {
            String edgeId = edge.getEdgeId();
            edgeId = edgeId.replace("-", "(pp)");
            edgeHits.put(edgeId, 0);
        }

        for (IKPMResultItem result : results.getResults()) {
            //TODO: unionNodeset
            for (String nodeId : result.getUnionNodeSet().keySet()) {
                nodeHits.put(nodeId, nodeHits.get(nodeId) + 1);
            }
            // TODO: check whether change from "getVisitedNodes" to getUnionNodeset is equal
            List<GeneEdge> edges = g.getConnectingEdges(result.getUnionNodeSet().values());
            for (GeneEdge edge : edges) {
                String edgeId = edge.getEdgeId();
                edgeId = edgeId.replace("-", "(pp)");
                edgeHits.put(edgeId, edgeHits.get(edgeId) + 1);
            }
        }

        for (String nodeId : nodeHits.keySet()) {
            int hits = nodeHits.get(nodeId);
            if (hits > maxNode) {
                maxNode = hits;
            }
        }

        for (String edgeId : edgeHits.keySet()) {
            int hits = edgeHits.get(edgeId);
            if (hits > maxEdge) {
                maxEdge = hits;
            }
        }

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(nodeFileName));
            bw.write("id" + "\t");
            bw.write("hits" + "\t");
            bw.write("hits.normalized" + "\t");
            bw.write("diff.expr.cases" + "\t");
            bw.write("is.exception" + "\n");

            for (GeneNode node : g.getVertices()) {
                String nodeId = node.getNodeId();
                bw.write(nodeId + "\t");
                int hits = nodeHits.get(nodeId);
                double avg = (double) hits / (double) maxNode;
                bw.write(hits + "\t");
                bw.write(avg + "\t");
                bw.write(node.getTotalDiffExpressedCases() + "\t");
                bw.write(node.isValid() + "\n");
            }
            bw.close();

            BufferedWriter bw2 = new BufferedWriter(new FileWriter(edgeFileName));
            bw2.write("id" + "\t");
            bw2.write("hits" + "\n");
            for (GeneEdge edge : g.getEdges()) {
                String edgeId = edge.getEdgeId();
                String[] fields = edgeId.split(" ");
                String from = fields[0].trim();
                String to = fields[2].trim();
                edgeId = edgeId.replace("-", "(pp)");
                bw2.write(edgeId + "\t");
                int hits = edgeHits.get(edgeId);
                double avg = (double) hits / (double) maxEdge;
                bw2.write(hits + "\t");
                bw2.write(avg + "\n");
                String edgeId2 = to + " (pp) " + from;
                bw2.write(edgeId2 + "\t");
                bw2.write(hits + "\n");
            }

            bw2.close();

        } catch (IOException ex) {
            Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
        }

//        System.out.println("MAX NODE HITS: " + maxNode);
//        System.out.println("MAX EDGE HITS: " + maxEdge);

    }

    public static void generateHistogramLValues(KPMGraph g) {
        if (!GENERATE_HIST_LVALUES) {
            return;
        }

        System.out.println("Generating Histogram L Values...");
        try {
            FileWriter fw = new FileWriter(new File(
                    "results/histogram-l-values.txt"));
            BufferedWriter out = new BufferedWriter(fw);

            int[] hist = new int[g.getVertices().iterator().next().getTotalDiffExpressedCases() + 1];
            for (int i = 0; i < hist.length; i++) {
                hist[i] = 0;
            }

            for (GeneNode node : g.getVertices()) {
                hist[node.getTotalDiffExpressedCases()] += 1;
            }
            // out.write(node.getNumDiffExpressedCases() + "\n");

            for (int i = 0; i < hist.length; i++) {
                out.write("l = " + i + " => " + hist[i] + " nodes.\n");
            }

            out.close();
            fw.close();
            System.out.println("Wrote histogram-l-values.txt...");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Generates a List that contains all genes as labeled in the node file.
     * Writes it to geneList.txt in the results folder.
     */
    @SuppressWarnings("unused")
    public static void generateGeneList(Map<String, String> nodeIdToGeneId) {
        if (!PRINT_GENE_LIST) {
            return;
        }
        try {
            FileWriter fw = new FileWriter(new File("results/geneList.txt"));
            BufferedWriter out = new BufferedWriter(fw);

            Collection<String> vertices = nodeIdToGeneId.keySet();
            for (String v : vertices) {
                out.write(nodeIdToGeneId.get(v));
                out.newLine();
            }
            out.close();
            fw.close();
            System.out.println("Wrote geneList.txt...");
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
            return;
        }
    }

    /*
     * Prints the results of the ACO algorithm into a file.
     */
    public static void writeResultsToFile2(String fileout, List<Result> results, KPMGraph g, KPMSettings kpmSettings) {

        try {
            FileOutputStream fos = new FileOutputStream(fileout, true);
            DataOutputStream dos = new DataOutputStream(fos);

            //k
            dos.writeInt(kpmSettings.GENE_EXCEPTIONS);
            // l
            //dos.writeInt(KPMParameters.CASE_EXCEPTIONS);
            // local search
            dos.writeUTF(kpmSettings.L_SEARCH.toString());
            // rho decay
            dos.writeUTF(kpmSettings.RHO_DECAY.toString());
            dos.writeBoolean(kpmSettings.ITERATION_BASED);
            dos.writeUTF(kpmSettings.ALGO.toString());

            // best fitness
            dos.writeInt(results.get(0).getFitness());
            // runtime
            dos.writeLong((System.nanoTime() - kpmSettings.STARTING_TIME) / 1000000000);

            dos.close();

        } catch (IOException ex) {
            ex.printStackTrace();
        }

    }

    public static void printGraphStatistics(KPMGraph g) {
        // Used for statistics...
        System.out.println("NUMBER OF NODES: " + g.getVertexCount());
        System.out.println("NUMBER OF INTERACTIONS: " + g.getEdgeCount());
        System.out.println("AVERAGE NODE DEGREE: " + g.getAverageDegree());
        //       KPMParameters.NUM_CASES = g.getVertices().iterator().next().getTotalCases();
        //       System.out.println("NUMBER OF CASES: " + KPMParameters.NUM_CASES);
        //       System.out.println("SEED: " + KPMParameters.SEED);
        generateHistogramLValues(g);
    }

    public static void printMatrixStatistics(Parser p, KPMGraph g) {
        Map<String, Integer> numCasesMap = p.getNumCasesMap();
        Map<String, Integer> numGenesMap = p.getNumGenesMap();
        Map<String, Set<String>> backNodesMap = p.getBackNodesByExpMap();
        Map<String, Set<String>> backGenesMap = p.getBackGenesMap();
        Map<String, Double> avgCasesMap = p.getAvgExpressedCasesMap();
        Map<String, Double> avgGenesMap = p.getAvgExpressedGenesMap();

        System.out.print("MATRIX_ID" + "\t");
        System.out.print("CASES" + "\t");
        System.out.print("GENES" + "\t");
        System.out.print("AVG_EXP_CASES_p_GENE" + "\t");
        System.out.print("AVG_EXP_GENES_p_CASE" + "\t");
        System.out.print("NODES_WITH_EXPRESSION_MAPPINGS" + "\n");
//      System.out.print("EXPRESSION_MAPPED" + "\n");

        DecimalFormat df = new DecimalFormat("#.00");
        int totalNodes = g.getVertexCount();
        for (String expId : numCasesMap.keySet()) {
            int totalGenes = numGenesMap.get(expId);
            System.out.print(String.format("%-" + "MATRIX_ID".length() + "s", expId) + "\t");
            System.out.print(String.format("%-" + "CASES".length() + "s", numCasesMap.get(expId)) + "\t");
            System.out.print(String.format("%-" + "GENES".length() + "s", totalGenes) + "\t");
            System.out.print(String.format("%-" + "AVG_EXP_CASES_p_GENE".length() + "s", avgCasesMap.get(expId)) + "\t");
            System.out.print(String.format("%-" + "AVG_EXP_GENES_p_CASE".length() + "s", avgGenesMap.get(expId)) + "\t");
            int backNodes = backNodesMap.get(expId).size();
            int mappedNodes = totalNodes - backNodes;
            double mappedNodesP = ((double) mappedNodes / (double) totalNodes) * 100.0;
            String mappedNodesPS = df.format(mappedNodesP);
            System.out.print(mappedNodes + " (" + mappedNodesPS + "%)" + "\n");

//            int backGenes = backGenesMap.get(expId).size();
//            int mappedGenes = totalGenes - backGenes;
//            double mappedGenesP = (double)mappedGenes / (double)totalGenes;
//            String mappedGenesPS = df.format(mappedGenesP);
//            System.out.print(String.valueOf(mappedGenes) + " (" + mappedGenesPS + "%)" + "\n");
        }


    }

    public static void writePathwaysStatsFile(String file, IKPMResultSet results, KPMGraph g) {
        System.out.println(">Writing pathways statistics file");
        DecimalFormat df = new DecimalFormat("#.###");
        int size = results.getResults().size();
        int zeros = (int) Math.floor(Math.log10(size)) + 1;
        String f = "";
        for (int i = 1; i <= zeros; i++) {
            f += "0";
        }
        DecimalFormat df2 = new DecimalFormat(f);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(file))) {
            bw.write("PATHWAY_ID" + "\t");
            bw.write("k" + "\t");
            bw.write("l" + "\t");
            bw.write("# NODES" + "\t");
            bw.write("# EDGES" + "\t");
            bw.write("AVG. DIFF. EXP. CASES" + "\t");
            bw.write("AVG. INFO. CONTENT" + "\n");
            int j = 0; // Variable encoding Pathway ID
            for (IKPMResultItem result : results.getResults()) { // List: order consistent
                int i = 0;
                for (Map<String, GeneNode> nodeSet : result.getAllComputedNodeSets()) { // List: order consistent
                    String id = df2.format(j + 1);
                    List<String[]> edges = g.getEdgesConnecting(nodeSet.values());
                    bw.write(id + "\t");
                    bw.write(result.getK() + "\t");
                    bw.write(result.getL() + "\t");
                    bw.write(nodeSet.size() + "\t");
                    bw.write(edges.size() + "\t");
                    bw.write(df.format(result.getResultsInfoTable()[i][3]) + "\t");
                    bw.write(df.format(result.getResultsInfoTable()[i][4]) + "\n");
                    i = i + 1;
                    j = j + 1;
                }
            }
        } catch (IOException ex) {
            Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    private static String checkFileName(String file) {
        String ret = file;
        int count = 0;
        while ((new File(ret)).exists()) {
            count++;
            String name = ret;
            if (count == 1) {
                int index = name.lastIndexOf(".");
                ret = name.substring(0, index) + "(" + String.valueOf(count) + ")"
                        + name.substring(index, name.length());
            } else {
                int aux = count - 1;
                ret = name.replace("(" + aux + ")", "(" + count + ")");
            }
        }
        return ret;
    }


    public static DecimalFormat getIntegerFormat(int size) {
        int zeros = (int) Math.floor(Math.log10(size)) + 1;
        String f = "";
        for (int i = 1; i <= zeros; i++) {
            f += "0";
        }
        DecimalFormat df = new DecimalFormat(f);
        return df;
    }

    public static DecimalFormat getDoubleFormat(int zeros) {
        String format = "#.";
        for (int i = 1; i <= zeros; i++) {
            format += "0";
        }

        return new DecimalFormat(format);
    }
}

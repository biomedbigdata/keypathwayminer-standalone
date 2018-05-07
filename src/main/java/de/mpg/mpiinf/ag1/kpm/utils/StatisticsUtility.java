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

import de.mpg.mpiinf.ag1.kpm.Parameters;
import dk.sdu.kpm.Algo;
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
    public static boolean GENERATE_MORE_STATISTICS = false;
    public static boolean GENERATE_HIST_LVALUES = false;

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

    /**
     * Prints the results of the ACO algorithm into a file.
     */
//    public static void writeResultsToFile(List<Result> results, KPMGraph g) {
//        java.util.Date today = new java.util.Date();
//        String stamp = new java.sql.Timestamp(today.getTime()).toString().replaceAll(" ", "_");
//        stamp = stamp.replaceAll(":", "_");
//
//        File directory = new File("results");
//        directory.mkdir();
//
//        String filename = stamp + "-summary.txt";
//        File f = new File("results", filename);
//        FileWriter fstream;
//        HTTStats stats = new HTTStats();
//
//        try {
//            fstream = new FileWriter(f);
//            BufferedWriter out = new BufferedWriter(fstream);
//            out.write("********** PARAMETERS ***************\n");
//            out.write("GENE EXCEPTIONS: " + Globals.GENE_EXCEPTIONS + "\n");
//            out.write("CASE EXCEPTIONS: " + Globals.CASE_EXCEPTIONS + "\n");
//            //out.write("Data set: " + Globals.DIRECTORY + "\n");
//            out.write("Number of Solutions per Iteration: "
//                    + Globals.NUMBER_OF_SOLUTIONS_PER_ITERATION + "\n");
//            out.write("alpha: " + Globals.ALPHA + "\n");
//            out.write("beta: " + Globals.BETA + "\n");
//            out.write("rho: " + Globals.RHO + "\n");
//            out.write("randomness seed: " + Globals.SEED + "\n");
//            out.write("iterations: " + Globals.MAX_ITERATIONS + "\n");
//            out.write("algorithm: " + Globals.ALGO + "\n");
//            out.write("runtime (in nanoseconds): "
//                    + (System.nanoTime() - Globals.STARTING_TIME) + "\n");
//            out.write("*************************************\n");
//
//            int rank = 1;
//            for (Result result : results) {
//                if (rank > Globals.NUM_SOLUTIONS) {
//                    break;
//                }
//
//                HashMap<String, LinkedList<String>> lists = stats.getStats(result.getVisitedNodes().keySet());
//                out.write("Rank: " + rank + " Nodes: "
//                        + result.getVisitedNodes().size()
//                        + " Case Exceptions: "
//                        + result.getNonDifferentiallyExpressedCases()
//                        + " Avg. Expressed: "
//                        + result.getAverageDiffExpressedCases()
//                        + " Info gain: " + result.getInformationGainExpressed()
//                        + "\n");
//                out.write("Genes in KEGG: " + lists.get(HTTStats.KEGGID).size()
//                        + ", CALCIUM: " + lists.get(HTTStats.CALCIUMID).size()
//                        + ", INTERACTION: "
//                        + lists.get(HTTStats.INTERACTID).size() + "\n");
//                out.write("\n");
//                result.flagExceptionNodes();
//                for (GeneNode node : result.getVisitedNodes().values()) {
//                    String isValid = "";
//                    if (node.isValid()) {
//                        isValid = "\tYES";
//                    } else {
//                        isValid = "\tNO";
//                    }
//                    out.write(node.getSymbol() + " " + isValid + "\t"
//                            + node.getTotalNoDiffCases() + "\n");
//                }
//                out.write("\n");
//                rank++;
//            }
//            out.close();
//            rank = 0;
//            for (Result result : results) {
//                rank++;
//                if (rank > Globals.NUM_SOLUTIONS) {
//                    break;
//                }
//
//                BufferedWriter nodeFile = new BufferedWriter(new FileWriter(
//                        "results/" + stamp + "-rank" + rank + "-nodes.txt"));
//                BufferedWriter edgeFile = new BufferedWriter(new FileWriter(
//                        "results/" + stamp + "-rank-" + rank + "-edges.txt"));
//
//                for (GeneNode node : result.getVisitedNodes().values()) {
//                    nodeFile.write(node.getNodeId() + " " + node.getSymbol()
//                            + "\n");
//                }
//
//                for (GeneNode node1 : result.getVisitedNodes().values()) {
//                    for (GeneNode node2 : result.getVisitedNodes().values()) {
//                        if (g.findEdge(node1, node2) != null) {
//                            String fromNodeId = node1.getNodeId();
//                            String toNodeId = node2.getNodeId();
//                            edgeFile.write(fromNodeId + " " + toNodeId + "\n");
//                        }
//                    }
//                }
//                nodeFile.close();
//                edgeFile.close();
//            }
//        } catch (IOException ex) {
//            ex.printStackTrace();
//        }
//
//    }
    /**
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
//        System.out.print("EXPRESSION_MAPPED" + "\n");

        DecimalFormat df = new DecimalFormat("#.00");
        int totalNodes = g.getVertexCount();
        for (String expId : numCasesMap.keySet()) {
            int totalGenes = numGenesMap.get(expId);
            System.out.print(expId + "\t");
            System.out.print(String.valueOf(numCasesMap.get(expId)) + "\t");
            System.out.print(String.valueOf(totalGenes) + "\t");
            System.out.print(String.valueOf(avgCasesMap.get(expId)) + "\t");
            System.out.print(String.valueOf(avgGenesMap.get(expId)) + "\t");
            int backNodes = backNodesMap.get(expId).size();
            int mappedNodes = totalNodes - backNodes;
            double mappedNodesP = ((double)mappedNodes / (double)totalNodes) * 100.0;
            String mappedNodesPS = df.format(mappedNodesP);
            System.out.print(String.valueOf(mappedNodes) + " (" + mappedNodesPS + "%)" + "\n");
//            int backGenes = backGenesMap.get(expId).size();
//            int mappedGenes = totalGenes - backGenes;
//            double mappedGenesP = (double)mappedGenes / (double)totalGenes;
//            String mappedGenesPS = df.format(mappedGenesP);
//            System.out.print(String.valueOf(mappedGenes) + " (" + mappedGenesPS + "%)" + "\n");
        }



    }

    public static void writeDatasetsStats(String file, Parser p, KPMGraph g) {
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
//            bw.write("GENES_NOT_MAPPED" + "\n");

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
                bw.write(String.valueOf(mappedNodes) + " (" + mappedNodesPS + "%)" + "\n");

                //                bw.write(String.valueOf(backGenesMap.get(expId).size()) + "\n");
                System.out.println("WROTE for dataset: " + expId);
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

    public static void writePathwaysFile(String file, List<Result> results, KPMGraph g, KPMSettings kpmSettings) {
        BufferedWriter bw = null;
        try {
            int size = results.size();
            int zeros = (int) Math.floor(Math.log10(size)) + 1;
            String f = "";
            for (int i = 1; i <= zeros; i++) {
                f += "0";
            }
            DecimalFormat df = new DecimalFormat(f);
            bw = new BufferedWriter(new FileWriter(file));

            boolean isInes = kpmSettings.ALGO == Algo.GREEDY || kpmSettings.ALGO == Algo.LCG
                    || kpmSettings.ALGO == Algo.OPTIMAL;

            for (int i = 0; i < results.size(); i++) {
                Result result = results.get(i);
                Map<String, GeneNode> nodes = result.getVisitedNodes();
                String id = df.format(i + 1);
                bw.write(id + "\n");
                bw.write("NODES" + "\t" + "IS EXCEPTION?" + "\n");
                if (isInes) {
                    for (GeneNode node : nodes.values()) {
                        bw.write(node.getNodeId() + "\t" + !node.isValid() + "\n");
                    }
                } else {
                    for (GeneNode node : nodes.values()) {
                        bw.write(node.getNodeId() + "\n");
                    }
                }
                List<String[]> edges = g.getEdgesConnecting(nodes.values());
                bw.write("EDGES" + "\n");
                for (String[] edge : edges) {
                    bw.write(edge[0] + "\t" + "pp" + "\t" + edge[1] + "\n");
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
    
    public static void writeIndividualPathwayFiles(String dir, List<Result> results, KPMGraph g, Parameters params) {
        BufferedWriter bw1 = null;
        BufferedWriter bw2 = null;
        
        //try {
            int size = results.size();
            int zeros = (int) Math.floor(Math.log10(size)) + 1;
            String f = "";
            for (int i = 1; i <= zeros; i++) {
                f += "0";
            }
            DecimalFormat df = new DecimalFormat(f);
            
            for (int i = 0; i < results.size(); i++) {
                Result result = results.get(i);
                String id = df.format(i+1);
                String nodeFile = dir + params.FILE_SEPARATOR + "Pathway-" + id 
                        + "-NODES-" + params.SUFFIX + params.FILE_EXTENSION;
                nodeFile = checkFileName(nodeFile);
                String edgeFile = dir + params.FILE_SEPARATOR + "Pathway-" + id 
                        + "-INTERACTIONS-" + params.SUFFIX + params.FILE_EXTENSION;
                edgeFile = checkFileName(edgeFile);
                
                
                Set<String> nodes = result.getVisitedNodes().keySet();
                List<String[]> edges = g.getEdgesConnecting(result.getVisitedNodes().values());
                try {
                    bw1 = new BufferedWriter(new FileWriter(nodeFile));
                    bw2 = new BufferedWriter(new FileWriter(edgeFile));
                    for (String nodeId: nodes) {
                        bw1.write(nodeId);
                        bw1.newLine();
                    }
                    bw1.close(); 
                    for (String[] edge: edges) {
                        bw2.write(edge[0] + "\t" + edge[1]);
                        bw2.newLine();
                    }
                    bw2.close();
                } catch (IOException ex) {
                    Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
                } finally {
                    try {
                        bw1.close();
                        bw2.close();
                    } catch (IOException ex) {
                        Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
                    }                    
                }
                   
                
                
            }
    }

    public static void writePathwaysStatsFile(String file, List<Result> results, KPMGraph g) {
        DecimalFormat df = new DecimalFormat("#.###");
        int size = results.size();
        int zeros = (int) Math.floor(Math.log10(size)) + 1;
        String f = "";
        for (int i = 1; i <= zeros; i++) {
            f += "0";
        }
        DecimalFormat df2 = new DecimalFormat(f);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));
            bw.write("PATHWAY ID" + "\t");
            bw.write("# NODES" + "\t");
            bw.write("# EDGES" + "\t");
            bw.write("AVG. DIFF. EXP. CASES" + "\t");
            bw.write("AVG. INFO. CONTENT" + "\n");
            for (int i = 0; i < size; i++) {
                String id = df2.format(i + 1);
                Result result = results.get(i);
                Map<String, GeneNode> nodes = result.getVisitedNodes();
                List<String[]> edges = g.getEdgesConnecting(nodes.values());
                bw.write(id + "\t");
                bw.write(nodes.size() + "\t");
                bw.write(edges.size() + "\t");
                bw.write(df.format(result.getAverageDiffExpressedCases()) + "\t");
                bw.write(df.format(result.getInformationGainExpressed()) + "\n");
            }
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(StatisticsUtility.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    public static void writeGeneStatsFile(String file, List<Result> results, KPMGraph g, KPMSettings kpmSettings) {
        DecimalFormat df = new DecimalFormat("#.00");
        boolean isInes = kpmSettings.ALGO == Algo.GREEDY || kpmSettings.ALGO == Algo.LCG
                || kpmSettings.ALGO == Algo.OPTIMAL;
        BufferedWriter bw = null;

        try {
            Map<String, Integer> nodeCount = new HashMap<String, Integer>();
            for (Result result : results) {
                for (String nodeId : result.getVisitedNodes().keySet()) {
                    if (!nodeCount.containsKey(nodeId)) {
                        nodeCount.put(nodeId, 1);
                    } else {
                        nodeCount.put(nodeId, nodeCount.get(nodeId) + 1);
                    }
                }
            }
            bw = new BufferedWriter(new FileWriter(file));
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
                if (nodeCount.containsKey(nodeId)) {
                    bw.write(nodeCount.get(nodeId) + "\t");
                } else {
                    bw.write("0" + "\t");
                }
                bw.write(df.format(node.getAverageExpressedCases()) + "\t");
                bw.write(String.valueOf(g.degree(node)));
                if (isInes) {
                    bw.write("\t" + String.valueOf(!node.isValid()));
                }
                bw.write("\n");
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

    public static void writeSummaryFile(String file, IKPMResultSet results, KPMSettings kpmSettings) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));

            bw.write("GENE EXCEPTIONS: " + kpmSettings.GENE_EXCEPTIONS + "\n");
            bw.write("CASE EXCEPTIONS" + "\t" + "MATRIX" + "\n");

            for (String expId : kpmSettings.CASE_EXCEPTIONS_MAP.keySet()) {
                bw.write(kpmSettings.CASE_EXCEPTIONS_MAP.get(expId)
                        + kpmSettings.MATRIX_FILES_MAP.get(expId) + "\n");
            }
            if (kpmSettings.ALGO == Algo.GREEDY) {
                bw.write("STRATEGY: " + "\t" + "INES" + "\n");
                bw.write("ALGORITHM: " + "\t" + "GREEDY" + "\n");
            } else if (kpmSettings.ALGO == Algo.LCG) {
                bw.write("STRATEGY: " + "\t" + "INES" + "\n");
                bw.write("ALGORITHM: " + "\t" + "ACO" + "\n");
            } else if (kpmSettings.ALGO == Algo.OPTIMAL) {
                bw.write("STRATEGY: " + "\t" + "INES" + "\n");
                bw.write("ALGORITHM: " + "\t" + "OPTIMAL" + "\n");
            } else if (kpmSettings.ALGO == Algo.EXCEPTIONSUMGREEDY) {
                bw.write("STRATEGY: " + "\t" + "GLONE" + "\n");
                bw.write("ALGORITHM: " + "\t" + "GREEDY" + "\n");
            } else if (kpmSettings.ALGO == Algo.EXCEPTIONSUMACO) {
                bw.write("STRATEGY: " + "\t" + "GLONE" + "\n");
                bw.write("ALGORITHM: " + "\t" + "ACO" + "\n");
            } else if (kpmSettings.ALGO == Algo.EXCEPTIONSUMOPTIMAL) {
                bw.write("STRATEGY: " + "\t" + "GLONE" + "\n");
                bw.write("ALGORITHM: " + "\t" + "OPTIMAL" + "\n");
            }
            //TODO: check whether this is actually the best result aka. if resultList is ordered
            IKPMResultItem best = results.getResults().get(0);
            bw.write("SIZE OF BEST RESULT: " + "\t" + best.getResultsInfoTable()[0][2] + "\n");
            bw.write("TOTAL RUNNING TIME: " + "\t" + kpmSettings.TOTAL_RUNNING_TIME + "\n");
            bw.close();
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
    
    public static void writeResultsStats(List<Result> results, KPMGraph g, Parameters params) {
        Map<String, Integer> nodeHits = new HashMap<String, Integer>();
        Map<String, Integer> edgeHits = new HashMap<String, Integer>();
        String nodeFileName = params.RESULTS_FOLDER + params.FILE_SEPARATOR 
                + "node_pathway_hits-" + params.SUFFIX + params.FILE_EXTENSION;
        String edgeFileName = params.RESULTS_FOLDER + params.FILE_SEPARATOR 
                + "edge_pathway_hits-" + params.SUFFIX + params.FILE_EXTENSION;
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

        for (Result result : results) {
            for (String nodeId : result.getVisitedNodes().keySet()) {
                nodeHits.put(nodeId, nodeHits.get(nodeId) + 1);
            }
            List<GeneEdge> edges = g.getConnectingEdges(result.getVisitedNodes().values());
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

        System.out.println("MAX NODE HITS: " + maxNode + "\n");
        System.out.println("MAX EDGE HITS: " + maxEdge + "\n");

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

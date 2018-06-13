package de.mpg.mpiinf.ag1.kpm.utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;


import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.graph.KPMGraph;

/**
 *
 * @author nalcaraz
 */
public class Parser {

    String graphFile;
    HashMap<String, String> expressionFiles;
    HashMap<String, Set<String>> backNodesMap;
    HashMap<String, Set<String>> backNodesByExpMap;
    HashMap<String, Set<String>> backGenesMap;
    String sep;
    HashMap<String, Integer> numCasesMap;
    HashMap<String, Integer> numGenesMap;
    HashMap<String, Double> avgExpressedCasesMap;
    HashMap<String, Double> avgExpressedGenesMap;
    HashMap<String, Integer> totalExpressedMap;
    Map<String, String> pvalueFilesMap;


    Map<String, Double> threshold;
    public boolean is_binary_matrix;
    public Comparator comparator;
    public boolean use_double_values;

    public static final String TAB = "\t";
    public static final String COMMA = ",";
    public static final String SPACE = " ";

    public Parser() {
        graphFile = "";
        expressionFiles = new HashMap<String, String>();
    }

    public Parser(String graphFile, Map<String, String> expressionFiles, String sep,
                  Map<String, Double> pvalueCutoffs, boolean is_binary_matrix, Comparator comparator,
                  Map<String, String> pValueFiles, boolean use_double_values) {
        this.graphFile = graphFile;
        this.expressionFiles = new HashMap<String, String>(expressionFiles);
        backNodesMap = new HashMap<String, Set<String>>();
        backNodesByExpMap = new HashMap<String, Set<String>>();
        backGenesMap = new HashMap<String, Set<String>>();        
        numCasesMap = new HashMap<String, Integer>();
        numGenesMap = new HashMap<String, Integer>();        
        avgExpressedCasesMap = new HashMap<String, Double>();
        avgExpressedGenesMap = new HashMap<String, Double>();
        totalExpressedMap = new HashMap<String, Integer>();
        this.threshold = pvalueCutoffs;
        this.is_binary_matrix = is_binary_matrix;
        this.comparator = comparator;
        this.pvalueFilesMap = pValueFiles;
        this.use_double_values = use_double_values;
        
        for (String expId : expressionFiles.keySet()) {
            backGenesMap.put(expId, new HashSet());
            numCasesMap.put(expId, 0);
        }

        if (!(sep.equals(TAB) || sep.equals(COMMA) || sep.equals(COMMA))) {
            this.sep = TAB;
        } else {
            this.sep = sep;
        }
    }

    public Parser(String graphFile, Map<String, String> expressionFiles) {
        this.graphFile = graphFile;
        this.expressionFiles = new HashMap<String, String>(expressionFiles);
        backNodesMap = new HashMap<String, Set<String>>();
        backNodesByExpMap = new HashMap<String, Set<String>>();
        backGenesMap = new HashMap<String, Set<String>>();
        numCasesMap = new HashMap<String, Integer>();
        numGenesMap = new HashMap<String, Integer>();
        avgExpressedCasesMap = new HashMap<String, Double>();
        avgExpressedGenesMap = new HashMap<String, Double>();
        totalExpressedMap = new HashMap<String, Integer>();
 
        for (String expId : expressionFiles.keySet()) {
            backGenesMap.put(expId, new HashSet());
            numCasesMap.put(expId, 0);
        }

        sep = TAB;

    }

    public void setGraphFile(String graphFile) {
        this.graphFile = graphFile;
    }

    public void setExpressionFiles(Map<String, String> expressionFiles) {
        this.expressionFiles = new HashMap<String, String>(expressionFiles);
    }

    public void addExpressionFile(String id, String path) {
        expressionFiles.put(id, path);
    }

    public String getGraphFile() {
        return graphFile;
    }

    public Map<String, String> getExpressionFiles() {
        return expressionFiles;
    }

    public Map<String, Set<String>> getBackNodesMap() {
        return backNodesMap;
    }

    public Map<String, Set<String>> getBackGenesMap() {
        return backGenesMap;
    }

    public Map<String, Integer> getNumCasesMap() {
        return numCasesMap;
    }

    public KPMGraph createGraph(KPMSettings kpmSettings) {
        HashMap<String, String> nodeId2Symbol = new HashMap<String, String>();
        Map<String, Map<String, double[]>> expressionMap = new HashMap<>();
        LinkedList<String[]> edgeList = new LinkedList<String[]>();
        HashMap<String, Integer> without_exp = new HashMap<String, Integer>();
        HashSet<String> inNetwork = new HashSet<String>();
        HashMap<String, Map<String,Double>> pvals = new HashMap<String, Map<String,Double>>();
        for (String fileId : expressionFiles.keySet()) {
            numCasesMap.put(fileId, 0);
            without_exp.put(fileId, 0);
        }
        try {

            String line = "";
            BufferedReader graphReader = new BufferedReader(new FileReader(graphFile));
            int cont = 0;
            while ((line = graphReader.readLine()) != null) {
                String[] fields = line.split("\t");
                String id1 = fields[0].trim();
                nodeId2Symbol.put(id1, id1);
                cont++;
                if (fields.length < 3) {
                    System.out.println("LINE NUMBER:" + cont);
                    System.out.print(line);
                }
                String id2 = fields[2].trim();
                nodeId2Symbol.put(id2, id2);

                String[] edge = {id1, id2};
                edgeList.add(edge);
                inNetwork.add(id1);
                inNetwork.add(id2);
            }
            for (String fileId : expressionFiles.keySet()) {
                int totalExp = 0;
                int numCases = 0;
                int numGenes = 0;
                HashMap<String, int[]> nodeId2Expression = new HashMap<String, int[]>();
                Set<String> inExp = new HashSet<String>();

                BufferedReader expressionReader =
                        new BufferedReader(new FileReader(expressionFiles.get(fileId)));


                if(kpmSettings.MATRIX_FILES_HAVE_HEADER){
                    expressionReader.readLine();
                }
                while ((line = expressionReader.readLine()) != null) {
                    numGenes++;
                    String[] fields = line.split("\t");
                    String nodeId = fields[0].trim();
                    inExp.add(nodeId);

                    double[] exp = new double[fields.length - 1];

                    // Currently any nonsensical value will be parsed to 0
                    if (is_binary_matrix) {
                        for (int i = 1; i < fields.length; i++) {
                            String val = fields[i].trim();
                            if (val.equals("1")) {
                                exp[i - 1] = 1;
                                totalExp++;
                            } else if (val.equals("-1")) {
                                exp[i - 1] = -1;
                                totalExp++;
                            } else if (val.equals("0")) {
                                exp[i - 1] = 0;
                            } else {
                                exp[i - 1] = 0;
                            }
                        }
                    }
                    // do not convert matrix into binary. Differential expression based on p-value
                    // is still counted
                    else if (use_double_values) {
                        for (int i = 1; i < fields.length; i++) {
                            String val = fields[i].trim();
                            exp[i - 1] = Double.parseDouble(val);
                            if (Comparison.evalate(Double.parseDouble(val), threshold.get(fileId), comparator)) {
                                totalExp++;
                            }
                        }
                    }
                    // TODO: make it work for different thresholds and different comparators. PRIO: Low
                    else {
                        for (int i = 1; i < fields.length; i++) {
                            String val = fields[i].trim();
                            if (Comparison.evalate(Double.parseDouble(val),  threshold.get(fileId), comparator)) {
                                exp[i - 1] = 1;
                                totalExp++;
                            } else {
                                exp[i - 1] = 0;
                            }
                        }
                    }

                    if (expressionMap.containsKey(nodeId)) {
                        expressionMap.get(nodeId).put(fileId, exp);
                    } else {
                        Map<String, double[]> aux = new HashMap<>();
                        aux.put(fileId, exp);
                        expressionMap.put(nodeId, aux);
                    }
                    numCases = exp.length;
                    numCasesMap.put(fileId, numCases);
                }
                totalExpressedMap.put(fileId, totalExp);
                double avgExpCases = 0;
                double avgExpGenes = 0;
                if (totalExp > 0) {
                    avgExpCases = (double) numCases / (double) totalExp;
                    avgExpGenes = (double) numGenes / (double) totalExp;
                }
                numGenesMap.put(fileId, inExp.size());
                avgExpressedCasesMap.put(fileId, avgExpCases);
                avgExpressedGenesMap.put(fileId, avgExpGenes);
                Set<String> bckN = new HashSet(inNetwork);
                Set<String> bckG = new HashSet(inExp);
                for (String id : inNetwork) {
                    if (inExp.contains(id)) {
                        bckN.remove(id);
                    }
                }
                for (String id : inExp) {
                    if (inNetwork.contains(id)) {
                        bckG.remove(id);
                    }
                }

                backNodesByExpMap.put(fileId, bckN);
                backGenesMap.put(fileId, bckG);
                expressionReader.close();

                if (!is_binary_matrix && use_double_values) {
                    try (BufferedReader pvalReader = new BufferedReader(new FileReader(kpmSettings.PVALUE_FILES_MAP.get(fileId)))) {
                        String ll = "";
                        while ((ll = pvalReader.readLine()) != null) {
                            String[] l = ll.split("\t");
                            pvals.putIfAbsent(l[0], new HashMap<String, Double>());
                            pvals.get(l[0]).put(fileId, Double.parseDouble(l[1]));
                        }

                    } catch (IOException ioe) {
                        Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, "Your p-value file produced an excpetion" +
                                "- check,if you have added a valid file", ioe);
                    } catch (NumberFormatException nfe) {
                        Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, "Your p-value file seems to" +
                                "have an incorrect format", nfe);
                    }
                }
            }
 


            graphReader.close();


        } catch (FileNotFoundException ex) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ioe) {
            Logger.getLogger(Parser.class.getName()).log(Level.SEVERE, null, ioe);
        }

        for (String nodeId : inNetwork) {
            if (expressionMap.containsKey(nodeId)) {
                Map<String, double[]> expMap = expressionMap.get(nodeId);
                for (String expId : expressionFiles.keySet()) {
                    if (!expMap.containsKey(expId)) {
                        if (backNodesMap.containsKey(nodeId)) {
                            backNodesMap.get(nodeId).add(expId);
                        } else {
                            HashSet<String> aux = new HashSet<String>();
                            aux.add(expId);
                            backNodesMap.put(nodeId, aux);
                        }

                    }
                }
            } else {
                if (backNodesMap.containsKey(nodeId)) {
                    backNodesMap.get(nodeId).addAll(expressionFiles.keySet());
                } else {
                    HashSet<String> aux = new HashSet<String>();
                    aux.addAll(expressionFiles.keySet());
                    backNodesMap.put(nodeId, aux);
                }
            }
        }

//        for (String expId : expressionFiles.keySet()) {
//            for (String nodeId : expressionMap.keySet()) {
//                if (!expressionMap.get(nodeId).containsKey(expId)) {
//                    backGenesMap.get(expId).add(nodeId);
//                }
//            }
//        }

 //       System.out.println("Nodes NOT mapped to AT LEAST 1 expression vector: " + backNodesMap.size());
 //       System.out.println("Expression vectors not mapped to a node in the network: ");
//        for (String fileId : expressionFiles.keySet()) {
//            Set<String> aux = backGenesMap.get(fileId);
//            System.out.println("  " + fileId + ": " + aux.size());
//        }

        kpmSettings.NUM_CASES_MAP = numCasesMap;
        kpmSettings.NUM_STUDIES = numCasesMap.size();
        return new KPMGraph(expressionMap, edgeList, nodeId2Symbol, backNodesMap, backGenesMap, kpmSettings.NUM_CASES_MAP, pvals, use_double_values);
    }

    public boolean isNumber(String input) {
        return input.matches("((-|\\+)?[0-9]+(\\.[0-9]+)?)+");
    }

    public Map<String, Double> getAvgExpressedCasesMap() {
        return avgExpressedCasesMap;
    }

    public Map<String, Double> getAvgExpressedGenesMap() {
        return avgExpressedGenesMap;
    }

    public Map<String, Integer> getTotalExpressedMap() {
        return totalExpressedMap;
    }

    public Map<String, Set<String>> getBackNodesByExpMap() {
        return backNodesByExpMap;
    }

    public Map<String, Integer> getNumGenesMap() {
        return numGenesMap;
    }
    
    
}

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.main;

import de.mpg.mpiinf.ag1.kpm.KPMRunHandler;
import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.Program;
import de.mpg.mpiinf.ag1.kpm.parsers.ArgsParametersParser;
import de.mpg.mpiinf.ag1.kpm.parsers.InputFileParser;
import de.mpg.mpiinf.ag1.kpm.parsers.PropertiesFileParser;
import de.mpg.mpiinf.ag1.kpm.shortestpath.ShortestPathAlgorithms;
import dk.sdu.kpm.KPMSettings;

/**
 * @author nalcaraz
 */

public class Main {

    private volatile KPMSettings kpmSettings;

    public static void main(String[] args) throws Exception {
        Main main = new Main();
        main.run(args);
    }

    public void run(String[] args) {
        //Used memory in mb
        System.out.println("Used Memory:" + (Runtime.getRuntime().freeMemory()) / (1024 * 1024));

        //Path to folder with input files, output files .
        String datasetFolder = "resources/";
        System.out.println("Properties file:" + "kpm.properties");

        //KPM default settings
        kpmSettings = new KPMSettings();

        //Parse parameters from the kpm.properties file, initialize default parameters
        PropertiesFileParser pfp = new PropertiesFileParser(kpmSettings);
        Parameters params = pfp.parse(datasetFolder);

        try {
            //Read commandline arguments
            ArgsParametersParser argsParser = new ArgsParametersParser(kpmSettings);
            params = argsParser.parse(args, params);

            //Check if all input files are correct
            InputFileParser ifp = new InputFileParser(kpmSettings);
            params = ifp.parse(params);

            //Start KeyPathwayMiner
            start(params);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void start(Parameters params) {
        if (params.PROGRAM == Program.SP) {
            System.out.println("PROGRAM SELECTED: Shortest Paths");
            ShortestPathAlgorithms.shortestPathways(params.GRAPH_FILE, params);
        } else {
            KPMRunHandler kpmHandler = new KPMRunHandler(kpmSettings);
            kpmHandler.runBatch(params);
        }
    }
}

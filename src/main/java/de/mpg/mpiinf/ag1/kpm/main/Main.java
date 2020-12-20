/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.main;

import de.mpg.mpiinf.ag1.kpm.output.OutputSettings;
import de.mpg.mpiinf.ag1.kpm.parsers.ArgsParametersParser;
import de.mpg.mpiinf.ag1.kpm.parsers.InputFileParser;
import de.mpg.mpiinf.ag1.kpm.parsers.PropertiesFileParser;
import de.mpg.mpiinf.ag1.kpm.shortestpath.ShortestPathAlgorithms;
import de.mpg.mpiinf.ag1.kpm.utils.Program;
import dk.sdu.kpm.KPMSettings;

/**
 * @author nalcaraz
 */

public class Main {

    public static void main(String[] args) {
        Main main = new Main();
        main.runR(args, "resources/","kpm.properties");
    }

    //Runner method for calling KeyPathwayMiner
    public void run(String[] args) {
        String datasetFolder = "resources/";
        String propertiesFile = "kpm.properties";

        System.out.println("********** KEYPATHWAYMINER STANDALONE VERSION 5 **********");
        System.out.println("Properties file: " + propertiesFile);
        System.out.println("Default datasets from: " + datasetFolder);
        kpm(args, datasetFolder, propertiesFile);
    }

    //Runner method for calling KeyPathwayMiner from R
    public String runR(String[] args, String datasetFolder, String propertiesFile) {
        System.out.println("********** KEYPATHWAYMINER STANDALONE VERSION 5 VIA R **********");
        System.out.println("Properties file: " + propertiesFile);
        System.out.println("Default datasets from: " + datasetFolder);
        kpm(args, datasetFolder, propertiesFile);
        return OutputSettings.CURRENT_OUTPUT_FOLDER;
    }

    // Prepares objects, parses arguments and runs KeyPathwayMiner
    private void kpm(String[] args, String datasetFolder, String propertiesFile) {
        //KPM default settings
        KPMSettings kpmSettings = new KPMSettings();

        //Parse parameters from the kpm.properties file
        PropertiesFileParser pfp = new PropertiesFileParser(kpmSettings, propertiesFile);
        Parameters params = pfp.parse(datasetFolder);

        try {
            //Read commandline arguments
            ArgsParametersParser argsParser = new ArgsParametersParser(kpmSettings, params);
            params = argsParser.parse(args);

            //Check files and set output file names
            InputFileParser ifp = new InputFileParser(kpmSettings);
            params = ifp.parse(params);

            //Start KeyPathwayMiner
            if (params.PROGRAM == Program.SP) {
                System.out.println("Program selected: Shortest Paths");
                ShortestPathAlgorithms.shortestPathways(params.GRAPH_FILE, params);
            } else {
                KPMRunHandler kpmHandler = new KPMRunHandler(kpmSettings);
                kpmHandler.runBatch(params);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


}

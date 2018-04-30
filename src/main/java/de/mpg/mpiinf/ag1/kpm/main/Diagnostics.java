package de.mpg.mpiinf.ag1.kpm.main;
import de.mpg.mpiinf.ag1.kpm.main.Main;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.stream.Stream;

/** This class assesses different parameter settings for KMP standalone
 */
public class Diagnostics {

    /* Diagnostics:
     * loop through different KMP versions with settings from input file args[0]
     * */

    public static void main(String[] args) {
        String filename = args[0];
        String logfile = args[1];
        Diagnostics d = new Diagnostics();
        ArrayList<String[]> commandLineEmulator = d.fileReader(filename);
        for (String[] commandLine : commandLineEmulator) {
            try {
                Main main = new Main();
                main.run(commandLine);
            } catch (Exception e) {
                errorLogger(logfile, commandLine, e);
            }

        }
    }

    public ArrayList<String[]> fileReader(String filename) {
        ArrayList<String[]> paramEmulator = new ArrayList<String[]>();
        try (BufferedReader br = new BufferedReader(new FileReader(filename))) {
            String line = br.readLine();
            while (line != null){
                String[] split = line.split(" ");
                paramEmulator.add(split);
                line = br.readLine();
            }
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return(paramEmulator);
    }

    public static void errorLogger(String logfile, String [] commandLine, Exception e){
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(logfile));
            //e.printStackTrace();
            bw.write(e.getStackTrace().toString());
            for (String s : commandLine) {
                bw.write(s + " ");
            }

        }
        catch (IOException ioe){
            ioe.printStackTrace();
        }
    }
    public static void generateDiagnosticsParameterScript(){

    }
}
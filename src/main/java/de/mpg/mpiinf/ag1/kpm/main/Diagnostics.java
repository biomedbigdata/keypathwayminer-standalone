package de.mpg.mpiinf.ag1.kpm.main;
import cern.colt.matrix.doublealgo.Stencil;
import de.mpg.mpiinf.ag1.kpm.main.Main;

import java.io.*;
import java.nio.Buffer;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
//import java.util.stream.Stream;

/** This class assesses different parameter settings for KMP standalone
 */
public class Diagnostics {

    /* Diagnostics:
     * loop through different KMP versions with settings
     * */

    public static void main(String[] args) {
        String possibleParamFile =  args[0];
        String outputFile = args[1];
        String logfile = args[2];
        Diagnostics d = new Diagnostics();
        //outputfile for the Parameterscript is same as input file for commandLineEmulator
        generateDiagnosticsParameterScript(possibleParamFile, outputFile);
        ArrayList<String[]> commandLineEmulator = d.fileReader(outputFile);
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
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(logfile))) {
            bw.write(e.getStackTrace().toString());
            for (String s : commandLine) {
                bw.write(s + " ");
            }

        }
        catch (IOException ioe){
            ioe.printStackTrace();
        }
    }

    /**
     *
     * @param possibleParameterfile line wise format: -paramName options,no,space,comma,separated
     * @returns Filename of parameter file.
     */
    public static void generateDiagnosticsParameterScript(String possibleParameterfile, String outputFile){
        // read the file that contains all the parameters for the algorithms and their respective possible values
        // if they need to be modified. Examples for params that need to be modified:
        // -algo: ACO -> GREEDY
        // examples for params that need not to be modified
        // input files, output files (can be done in an automated fashion)
        HashMap<String,String[]> possibleParams = new HashMap<String, String[]>();
        ArrayList<String> flags = new ArrayList<String>();
        HashMap<String, String> noChoiceParams = new HashMap<String, String>();
        try{

            BufferedReader br = new BufferedReader(new FileReader(possibleParameterfile));
            String line = br.readLine();
            while(line != null){
                String[] paramsAndPossibleValues = line.split(" ");
                if(paramsAndPossibleValues.length==2 && paramsAndPossibleValues[1].split(",").length>1) {
                    possibleParams.put(paramsAndPossibleValues[0], paramsAndPossibleValues[1].split(","));
                }
                else if (paramsAndPossibleValues.length==2){
                    noChoiceParams.put(paramsAndPossibleValues[0], paramsAndPossibleValues[1]);
                }
                else {
                    flags.add(paramsAndPossibleValues[0]);
                }
                line = br.readLine();
            }
        }
        catch (IOException ioe){
            ioe.printStackTrace();
        }
        int total = 1;
        for(String a: possibleParams.keySet()){
            total*=possibleParams.get(a).length;
        }

        // Create array with variying calls
        String[] ar = possibleParams.keySet().toArray(new String[3]);
        ArrayList<String> re =createCall(ar , possibleParams, "",  new ArrayList<String>());

        // put everything together, create final calls and print it into a file
        StringBuilder noChoice = new StringBuilder();
        for (String s : noChoiceParams.keySet()) {
            noChoice.append(s);
            noChoice.append("=");
            noChoice.append(noChoiceParams.get(s));
            noChoice.append(" ");
        }
        System.out.print(noChoice.toString());
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {
            for(String s: re){
                //System.out.println(noChoice.toString()+s);
                //System.out.println(noChoice.insert(0, s).toString());
                bw.write(noChoice.toString()+s);
            }
        }
        catch (IOException ioe){
            ioe.printStackTrace();
        }
    }


    public static ArrayList<String> createCall(String[] params, HashMap<String, String[]> p, String result, ArrayList<String> re) {
        if(params.length == 1){
            for(String pa: p.get(params[0])) {
                re.add(result +params[0]+"="+pa+"\n");
            }
            //System.out.print(re);
            return(re);
        }
        if(params.length > 1) {
            for(String pa: p.get(params[0])) {
                createCall(Arrays.copyOfRange(params, 1, params.length), p, result +params[0]+"="+pa+ " ", re);
            }
            return(re);
        }

        return(re);
    }
}
package com.github.tDBN.cli;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import com.github.tDBN.dbn.DynamicBayesNet;
import com.github.tDBN.dbn.LLScoringFunction;
import com.github.tDBN.dbn.MDLScoringFunction;
import com.github.tDBN.dbn.Observations;
import com.github.tDBN.dbn.ObservationsStatic;
import com.github.tDBN.dbn.Scores;
import com.github.tDBN.utils.Utils;
import com.github.tDBN.dbn.ObservationsToInference;

/**
 * Class with main method, to perform learning and inference of a tDBN using both dynamic and static attributes.
 * 
 * This class improves LearnFromFile, to learn a tDBN only with dynamic attributes, adding to it MainWithStatic, 
 * to learn a tDBN also with static attributes, and also allowing the user to perform inference on the learned tDBN.
 * 
 * @author Tiago Leao
 * 
 */
public class Inference {

	@SuppressWarnings({ "static-access" })
	public static void main(String[] args) {

		// create Options object
		Options options = new Options();

		Option inputFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Input CSV file to be used for network learning.").withLongOpt("inputFile")
				.create("i");

		Option numParents = OptionBuilder.withArgName("int").hasArg()
				.withDescription("Maximum number of parents from preceding time-slice(s).").withLongOpt("numParents")
				.create("p");

		Option outputFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Writes output to <file>. If not supplied, output is written to terminal.")
				.withLongOpt("outputFile").create("o");

		Option rootNode = OptionBuilder.withArgName("int").hasArg()
				.withDescription("Root node of the intra-slice tree. By default, root is arbitrary.")
				.withLongOpt("root").create("r");

		Option scoringFunction = OptionBuilder.hasArg()
				.withDescription("Scoring function to be used, either MDL or LL. MDL is used by default.")
				.withLongOpt("scoringFunction").create("s");

		Option dotFormat = OptionBuilder
				.withDescription(
						"Outputs network in dot format, allowing direct redirection into Graphviz to visualize the graph.")
				.withLongOpt("dotFormat").create("d");

		Option compact = OptionBuilder
				.withDescription(
						"Outputs network in compact format, omitting intra-slice edges. Only works if specified together with -d and with --markovLag 1.")
				.withLongOpt("compact").create("c");

		Option maxMarkovLag = OptionBuilder
				.withArgName("int")
				.hasArg()
				.withDescription(
						"Maximum Markov lag to be considered, which is the longest distance between connected time-slices. Default is 1, allowing edges from one preceding slice.")
				.withLongOpt("markovLag").create("m");

		Option spanningTree = OptionBuilder
				.withDescription(
						"Forces intra-slice connectivity to be a tree instead of a forest, eventually producing a structure with a lower score.")
				.withLongOpt("spanning").create("sp");

		Option nonStationary = OptionBuilder
				.withDescription(
						"Learns a non-stationary network (one transition network per time transition). By default, a stationary DBN is learnt.")
				.withLongOpt("nonStationary").create("ns");

		Option parameters = OptionBuilder.withDescription("Learns and outputs the network parameters.")
				.withLongOpt("parameters").create("pm");
		
		// New parameters introduced by Tiago Leao come next:
		
		Option newObservation = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with the observations where inference should be done.")
				.withLongOpt("obsFile").create("obs");

		Option inferenceFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with variables to perform inference on.")
				.withLongOpt("inferenceFile").create("inf");
		
		Option inferenceFormat = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Format to present inference. Can be distrSampl, to give only a value sampled according to the distribution; mostProb, to give only the most probable value; or distrib, to give the full distribution, for each attribute specified (where distrSampl is applied to intermmediate nodes). Default is distrSampl.")
				.withLongOpt("inferenceFormat").create("infFmt");
		
		Option outputInferenceFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Writes inference output to <file>. If not supplied, inference output is written to terminal.")
				.withLongOpt("outputInferenceFile").create("outInf");
		
		Option outputTrajectoryFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Writes predicted trajectories to <file>. If not supplied, output is written to terminal.")
				.withLongOpt("outputTrajectoryFile").create("tf");

		Option trajMaxTimestep = OptionBuilder.withArgName("int").hasArg()
				.withDescription("Timestep until which trajectory is to be determined.")
				.withLongOpt("trajectory").create("t");
		
		Option inputStatic = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Input CSV file with static features to be used for network learning.").withLongOpt("inputStaticFile")
				.create("is");
		
		Option numStaticParents = OptionBuilder.withArgName("int").hasArg()
				.withDescription("Maximum number of static parents of a certain node (default = 2).").withLongOpt("numStaticParents")
				.create("b");
		
		Option newStaticObservation = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with the static observations to make inference.")
				.withLongOpt("obsStaticFile").create("obsStatic");

		Option fromObjectFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with the serialized object of the sdtDBN.")
				.withLongOpt("fromObjFile").create("fromFile");
		
		Option toObjectFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File in which the serialized object with the sdtDBN should be stored.")
				.withLongOpt("toObjFile").create("toFile");
				
		Option mustNotAppear_dynPast = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File that, for each node Xi[t], contains the dynamic nodes from t'<t that cannot be parents of each Xi[t].")
				.withLongOpt("mustNotAppear_dynPast").create("mNotA_dynPast");
		
		Option mustAppear_dynPast = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File that, for each node Xi[t], contains the dynamic nodes from t'<t that must be parents of each Xi[t].")
				.withLongOpt("mustAppear_dynPast").create("mA_dynPast");
		
		Option mustNotAppear_static = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File that, for each node Xi[t], contains the static nodes that cannot be parents of each Xi[t].")
				.withLongOpt("mustNotAppear_static").create("mNotA_static");
		
		Option mustAppear_static = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File that, for each node Xi[t], contains the static nodes that must be parents of each Xi[t].")
				.withLongOpt("mustAppear_static").create("mA_static");
				
		Option mustNotAppear_dynSameTimestep = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File that, for each node Xi[t], contains the dynamic nodes from t that cannot be parents of each Xi[t].")
				.withLongOpt("mustNotAppear_dynSameTimestep").create("mNotA_dynSame");
		
		Option mustAppear_dynSameTimestep = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File that, for each node Xi[t], contains the dynamic nodes from t that must be parents of each Xi[t].")
				.withLongOpt("mustAppear_dynSameTimestep").create("mA_dynSame");
				
		options.addOption(inputFile);
		options.addOption(numParents);
		options.addOption(outputFile);
		options.addOption(rootNode);
		options.addOption(scoringFunction);
		options.addOption(dotFormat);
		options.addOption(compact);
		options.addOption(maxMarkovLag);
		options.addOption(spanningTree);
		options.addOption(nonStationary);
		options.addOption(parameters);
		
		options.addOption(inferenceFile);
		options.addOption(newObservation);
		options.addOption(outputTrajectoryFile);
		options.addOption(trajMaxTimestep);
		options.addOption(outputInferenceFile);
		options.addOption(inferenceFormat);
		options.addOption(inputStatic);
		options.addOption(numStaticParents);
		options.addOption(newStaticObservation);

		options.addOption(fromObjectFile);
		options.addOption(toObjectFile);

		options.addOption(mustNotAppear_dynPast);
		options.addOption(mustAppear_dynPast);
		options.addOption(mustNotAppear_static);
		options.addOption(mustAppear_static);
		options.addOption(mustNotAppear_dynSameTimestep);
		options.addOption(mustAppear_dynSameTimestep);
		
		CommandLineParser parser = new GnuParser();
		try {

			CommandLine cmd = parser.parse(options, args);

			boolean verbose = !cmd.hasOption("d");
			boolean stationary = !cmd.hasOption("nonStationary");
			boolean spanning = cmd.hasOption("spanning");
			boolean printParameters = cmd.hasOption("parameters");
			
			boolean makeInference = cmd.hasOption("obsFile");
			boolean makeInferenceAtt = cmd.hasOption("inferenceFile");
			boolean getMostProbTraj = cmd.hasOption("trajectory");
			boolean definedTrajOutput = cmd.hasOption("outputTrajectoryFile");
			boolean definedOutputInferenceFile = cmd.hasOption("outputInferenceFile");
			boolean definedinferenceFormat = cmd.hasOption("inferenceFormat");
			boolean hasStatic = cmd.hasOption("inputStaticFile");
			boolean hasStaticObservFile = cmd.hasOption("obsStatic");

			boolean learnDBNfromObjFile  = cmd.hasOption("fromObjFile");
			boolean storeObjInFile = cmd.hasOption("toObjFile");
			
			if( learnDBNfromObjFile == false && ( cmd.hasOption("inputFile") == false || cmd.hasOption("numParents") == false ) ) {
				System.out.println("No file with DBN object was given and either inputFile or numParents (or both) not specified!!");
				System.out.println("Check the following usage:\n");
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp("sdtDBN", options);
				System.exit(1);
			}
			
			
			DynamicBayesNet dbn;
			int markovLag;
			Observations o = null;
			ObservationsStatic staticObservations = null;
			
			// Parse all the files with restrictions
			String mustNotAppear_dynPast_filename = cmd.hasOption("mNotA_dynPast") == true ? cmd.getOptionValue("mustNotAppear_dynPast") : null;
			String mustAppear_dynPast_filename = cmd.hasOption("mA_dynPast") == true ? cmd.getOptionValue("mustAppear_dynPast") : null;
			String mustNotAppear_static_filename = cmd.hasOption("mNotA_static") == true ? cmd.getOptionValue("mustNotAppear_static") : null;
			String mustAppear_static_filename = cmd.hasOption("mA_static") == true ? cmd.getOptionValue("mustAppear_static") : null;
			String mustNotAppear_dynSameTimestep_filename = cmd.hasOption("mNotA_dynSame") == true ? cmd.getOptionValue("mustNotAppear_dynSameTimestep") : null;
			String mustAppear_dynSameTimestep_filename = cmd.hasOption("mA_dynSame") == true ? cmd.getOptionValue("mustAppear_dynSameTimestep") : null;
			
			// If restricting the intra-slice connectivity, force it to be a spanning tree (instead of forest) so that all the restrictions are applied
			if(mustNotAppear_dynSameTimestep_filename != null || mustAppear_dynSameTimestep_filename != null) {
				spanning = true;
			}
			
			if(learnDBNfromObjFile == false) {

				// TODO: check sanity
				markovLag = Integer.parseInt(cmd.getOptionValue("m", "1"));
				int root = Integer.parseInt(cmd.getOptionValue("r", "-1"));
	
				o = new Observations(cmd.getOptionValue("i"), markovLag);
				
				// Fill staticObservations only if they are provided
				if(hasStatic == true)
					staticObservations = new ObservationsStatic(cmd.getOptionValue("is"), o.getSubjLinePerMtrx(), o.numTransitions(), o.getNumbSubjects());			
	
				Scores s = new Scores(o, Integer.parseInt(cmd.getOptionValue("p")), stationary, verbose, staticObservations, Integer.parseInt(cmd.getOptionValue("b", "2")), mustNotAppear_dynPast_filename, mustAppear_dynPast_filename, mustNotAppear_static_filename, mustAppear_static_filename, mustNotAppear_dynSameTimestep_filename, mustAppear_dynSameTimestep_filename );
				if (cmd.hasOption("s") && cmd.getOptionValue("s").equalsIgnoreCase("ll")) {
					if (verbose)
						System.out.println("Evaluating network with LL score.");
					s.evaluate(new LLScoringFunction());
				} else {
					if (verbose)
						System.out.println("Evaluating network with MDL score.");
					s.evaluate(new MDLScoringFunction());
				}
	
				// if (verbose)
				// System.out.println(s);
	
				if (verbose) {
					if (cmd.hasOption("r"))
						System.out.println("Root node specified: " + root);
					if (spanning)
						System.out.println("Finding a maximum spanning tree.");
					else
						System.out.println("Finding a maximum branching.");
				}
	
				dbn = s.toDBN(root, spanning);
	
				if (printParameters || storeObjInFile)
					dbn.learnParameters(o, stationary, staticObservations);
			
			} else { // Case where DBN must be learned from a given file with the object
				FileInputStream fi = new FileInputStream(new File(cmd.getOptionValue("fromObjFile")));
				ObjectInputStream oi = new ObjectInputStream(fi);
				
				// Read objects
				dbn = (DynamicBayesNet) oi.readObject();

				oi.close();
				fi.close();
				
				// Get parameters from the dbn learned
				markovLag = dbn.getMarkovLag();
				hasStatic = dbn.hasStaticAtts();
				stationary = dbn.isStationary();
				printParameters = true;
			}
			
			if(storeObjInFile == true) {
				FileOutputStream file = new FileOutputStream(new File(cmd.getOptionValue("toObjFile")));
				ObjectOutputStream object = new ObjectOutputStream(file);

				// Write objects to file
				object.writeObject(dbn);

				object.close();
				file.close();
			}

			String output;
			if (cmd.hasOption("d")) {
				if (cmd.hasOption("c") && markovLag == 1)
					output = dbn.toDot(true);
				else
					output = dbn.toDot(false);
			} else
				output = dbn.toString(printParameters);

			if (cmd.hasOption("o")) {
				try {
					Utils.writeToFile(cmd.getOptionValue("o"), output);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			} else {
				if (verbose) {
					System.out.println();
					System.out.println("-----------------");
					System.out.println();
				}
				System.out.println(output);
			}
			
			// Starting here is the inference part of the program:
			
			if(makeInference == true) {
				
				if(makeInferenceAtt == false && getMostProbTraj==false) {
					System.err.println("Desired inference not specified (neither a specific set of attributes nor a trajectory were defined!)");
					System.exit(1);
				}
				
				// Parameters were not learned before so they should be now to perform inference
				// If sdtDBN was given in file, it comes with parameters for sure
				if(printParameters == false && storeObjInFile == false)
					dbn.learnParameters(o, stationary, staticObservations);
				
				ObservationsToInference observToInference; // Check whether the observations to inference should be filled only with dynamic observations or also with sattic observations
				if(hasStatic == true && hasStaticObservFile){
					observToInference = new ObservationsToInference(cmd.getOptionValue("obs"), markovLag, dbn.getDynAttributes(), cmd.getOptionValue("obsStatic"), dbn.getStaticAttributes());
				} else {
					observToInference = new ObservationsToInference(cmd.getOptionValue("obs"), markovLag, dbn.getDynAttributes(), null, null);
				}
				
				int inferenceFormatAux = 0;
				if(definedinferenceFormat == true) {
					if(cmd.getOptionValue("inferenceFormat").equalsIgnoreCase("mostProb")) {
						inferenceFormatAux = 1;
					} else if(cmd.getOptionValue("inferenceFormat").equalsIgnoreCase("distrib")){
						inferenceFormatAux = 2;
					}	
				}
				
				if (getMostProbTraj == true) { // To perform the option -t
					
					int timestepMax = Integer.parseInt(cmd.getOptionValue("t", "-1"));
					if(timestepMax == -1) {
						System.err.println("Error parsing parameter -t");
						System.exit(1);
					}
					
					observToInference.getMostProbableTrajectory(timestepMax, dbn, stationary, inferenceFormatAux);
					observToInference.printMostProbableTrajectory(definedTrajOutput, cmd.getOptionValue("tf"));
				}
				
				
				if (makeInferenceAtt == true) { // To perform inference on specific attributes specified by user
					observToInference.parseAttributes(cmd.getOptionValue("inf"));
					
					if(definedOutputInferenceFile== true) {
						observToInference.makeInference(stationary, dbn, cmd.getOptionValue("outputInferenceFile"), inferenceFormatAux);
					} else {
						observToInference.makeInference(stationary, dbn, null, inferenceFormatAux);
					}
						
				}
			}
			
		} catch (FileNotFoundException e) {
			System.out.println("File not found");
			e.printStackTrace();
		} catch (IOException e) {
			System.out.println("Error initializing stream");
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("sdtDBN", options);
		}

	}
}

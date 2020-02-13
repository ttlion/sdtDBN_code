package com.github.tDBN.cli;

import java.io.FileNotFoundException;

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

public class Inference {

	@SuppressWarnings({ "static-access" })
	public static void main(String[] args) {

		// create Options object
		Options options = new Options();

		Option inputFile = OptionBuilder.withArgName("file").hasArg().isRequired()
				.withDescription("Input CSV file to be used for network learning.").withLongOpt("inputFile")
				.create("i");

		Option numParents = OptionBuilder.withArgName("int").hasArg().isRequired()
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
		
		Option newObservation = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with the observations where inference should be done")
				.withLongOpt("obsFile").create("obs");

		Option inferenceFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with variables to perform inference on")
				.withLongOpt("inferenceFile").create("inf");
		
		Option inferenceFormat = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Format to present inference. Can be mostProb or distrib, to give only the most probable value or the full distribution, for each attribute specified. Default is mostProb.")
				.withLongOpt("inferenceFormat").create("infFmt");
		
		Option outputInferenceFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File to write the inference performed")
				.withLongOpt("outputInferenceFile").create("outInf");
		
		Option outputTrajectoryFile = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Writes predicted trajectories to <file>. If not supplied, output is written to terminal.")
				.withLongOpt("outputTrajectoryFile").create("tf");

		Option trajMaxTimestep = OptionBuilder.withArgName("int").hasArg()
				.withDescription("Timestep until which trajectory is to be determined")
				.withLongOpt("trajectory").create("t");
		
		Option inputStatic = OptionBuilder.withArgName("file").hasArg()
				.withDescription("Input CSV file with static features to be used for network learning.").withLongOpt("inputStaticFile")
				.create("is");
		
		Option numStaticParents = OptionBuilder.withArgName("int").hasArg()
				.withDescription("Maximum number of static parents of a certain node (default = 2)").withLongOpt("numStaticParents")
				.create("b");
		
		Option newStaticObservation = OptionBuilder.withArgName("file").hasArg()
				.withDescription("File with the static observations to make inference")
				.withLongOpt("obsStaticFile").create("obsStatic");
				
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

			// TODO: check sanity
			int markovLag = Integer.parseInt(cmd.getOptionValue("m", "1"));
			int root = Integer.parseInt(cmd.getOptionValue("r", "-1"));

			Observations o = new Observations(cmd.getOptionValue("i"), markovLag);

			ObservationsStatic staticObservations = null;

			if(hasStatic == true)
				staticObservations = new ObservationsStatic(cmd.getOptionValue("is"), o.getSubjLinePerMtrx(), o.numTransitions(), o.getNumbSubjects());			

			Scores s = new Scores(o, Integer.parseInt(cmd.getOptionValue("p")), stationary, verbose, staticObservations, Integer.parseInt(cmd.getOptionValue("b", "2")));
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

			DynamicBayesNet dbn;

			if (verbose) {
				if (cmd.hasOption("r"))
					System.out.println("Root node specified: " + root);
				if (spanning)
					System.out.println("Finding a maximum spanning tree.");
				else
					System.out.println("Finding a maximum branching.");
			}

			dbn = s.toDBN(root, spanning);

			if (printParameters)
				dbn.learnParameters(o, stationary, staticObservations);

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

			if(makeInference == true) {
				
				if(makeInferenceAtt == false && getMostProbTraj==false) {
					System.err.println("Desired inference not specified (neither a specific set of attributes nor a trajectory were defined!)");
					System.exit(1);
				}
				
				if(printParameters == false) // Parameters were not learned before so they should be now to perform inference
					dbn.learnParameters(o, stationary, staticObservations);
				
				ObservationsToInference observToInference;
				if(staticObservations != null && hasStaticObservFile){
					observToInference = new ObservationsToInference(cmd.getOptionValue("obs"), markovLag, o.getAttributes(), cmd.getOptionValue("obsStatic"), staticObservations.getAttributes());
				} else {
					observToInference = new ObservationsToInference(cmd.getOptionValue("obs"), markovLag, o.getAttributes(), null, null);
				}
				
				if (getMostProbTraj == true) {
					
					int timestepMax = Integer.parseInt(cmd.getOptionValue("t", "-1"));
					if(timestepMax == -1) {
						System.err.println("Error parsing parameter -t");
						System.exit(1);
					}
					
					observToInference.getMostProbableTrajectory(timestepMax, dbn, stationary);
					observToInference.printMostProbableTrajectory(definedTrajOutput, cmd.getOptionValue("tf"));
				}
				
				boolean presentDistr = false;
				
				if (makeInferenceAtt == true) {
					observToInference.parseAttributes(cmd.getOptionValue("inf"));
					
					if(definedinferenceFormat == true && cmd.getOptionValue("inferenceFormat").equalsIgnoreCase("distrib") )
						presentDistr = true;
					
					if(definedOutputInferenceFile== true) {
						observToInference.makeInference(stationary, dbn, cmd.getOptionValue("outputInferenceFile"), presentDistr);
					} else {
						observToInference.makeInference(stationary, dbn, null, presentDistr);
					}
						
				}
			}
			
		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("tDBN", options);
		}

	}
}

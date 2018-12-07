package rnaSeqGUI;

import java.util.*;

import java.util.function.UnaryOperator;
import java.io.*;
import java.text.*;

import javafx.application.*;
import javafx.event.*;
import javafx.geometry.*;
import javafx.scene.*;
import javafx.scene.control.*;
import javafx.scene.control.Alert.*;
import javafx.scene.control.ScrollPane.*;
import javafx.scene.control.TextFormatter.*;
import javafx.scene.layout.*;
import javafx.scene.text.*;
import javafx.stage.*;

public class RNA_Seq_GUI extends Application {

	// GUI Controls
	private Text welcomeText, enterPWDText, currentPWDText, sampleNamesInputInstructions, fastqDumpInstructions,
			kallistoInstructions, singleEndInstructions, deseq2Instructions;
	private ScrollPane mainScrollPane;
	private HBox sampleNamesInputHBox, singleEndInputHBox;
	private VBox mainBox, fastqDumpOptionsVBox, kallistoInstructionsVBox, deseq2InstructionsVBox;
	private Scene createFileStructureStage;
	private TextField pWDTextFieldInput, fragmentLengthInput, standardDevInput;
	private TextArea executionLog, sampleNamesInputTextArea;
	private Button createFileStructureButton, clearPWDButton, enterSampleNamesButtonKallisto, generateFastQDumpButton,
			generateKallistoInstructionsButton, generateDeseq2InstructionsButton;
	private RadioButton downloadNCBI, noDownloadNCBI, singleEnd, pairedEnd;

	// Files
	private File pWD, pipelineFile, samplesAccession, kallisto, samplesKallisto, shellScriptsKallisto, deseq2,
			samplesDeseq2, shellScriptsDeseq2, deseq2_kallisto, deseq2_results;

	// Instance variables
	private final int WIDTH_OF_GUI = 800, HEIGHT_OF_GUI = 800;
	private final int WIDTH_OF_CONTROLS = WIDTH_OF_GUI - 20;
	private final int DEFAULT_SPACING = 10;
	private static int executionCounter = 1;
	private List<String> sampleListSorted;
	private boolean singleTPairedF;
	private String appendToRNASeqPipelineFileError;
	private double fragmentLength, standardDev;
	private final static LengthLexicoComparator LENGTH_LEXICO_COMPARATOR = new LengthLexicoComparator();

	// For text file input of instructions
	ClassLoader classLoader = RNA_Seq_GUI.class.getClassLoader();
	private InputStream deseq2InstructionsFile = classLoader.getResourceAsStream("deseq2Instructions.md"),
			fastqDumpInstructionsBeginFile = classLoader.getResourceAsStream("fastqDumpInstructionsBegin.md"),
			fastqDumpInstructionsEndFile = classLoader.getResourceAsStream("fastqDumpInstructionsEnd.md"),
			fastqDumpInstructionsScriptPairedFile = classLoader
					.getResourceAsStream("fastqDumpInstructionsScriptPaired.md"),
			fastqDumpInstructionsScriptSingleFile = classLoader
					.getResourceAsStream("fastqDumpInstructionsScriptSingle.md"),
			fragmentLengthStdDevInstructionsFile = classLoader
					.getResourceAsStream("fragmentLengthStdDevInstructions.md"),
			kallistoInstructionsBeginFile = classLoader.getResourceAsStream("kallistoInstructionsBegin.md"),
			kallistoInstructionsEndFile = classLoader.getResourceAsStream("kallistoInstructionsEnd.md"),
			kallistoInstructionsScriptPairedFile = classLoader
					.getResourceAsStream("kallistoInstructionsScriptPaired.md"),
			kallistoInstructionsScriptSingleFile = classLoader
					.getResourceAsStream("kallistoInstructionsScriptSingle.md"),
			moveFastqFilesInstructionsFile = classLoader.getResourceAsStream("moveFastqFilesInstructions.md");

	// Flags are final by nature of the FlagForCodeSubstitution being immutable
	private FlagForCodeSubstitution PWD_ABSOLUTEPATH_FLAG, KALLISTO_ABSOLUTEPATH_FLAG,
			KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH_FLAG, KALLISTO_SAMPLES_ABSOLUTEPATH_FLAG,
			ACCESSION_SAMPLES_ABSOLUTEPATH_FLAG, FRAGMENT_LENGTH_FLAG, STD_DEV_FLAG, DESEQ2_ABSOLUTEPATH_FLAG,
			DESEQ2_SHELLSCRIPTS_ABSOLUTEPATH_FLAG, DESEQ2_KALLISTO_SUBDIR_ABSOLUTEPATH_FLAG,
			DESEQ2_RESULTS_ABSOLUTEPATH_FLAG;
	private List<FlagForCodeSubstitution> flagList = new ArrayList<>();

	/*
	 * Comparator to order sample names first by length, then lexicographically
	 * ignoring case
	 */
	private static class LengthLexicoComparator implements Comparator<String> {

		@Override
		public int compare(String o1, String o2) {
			if (o1.length() == o2.length()) {
				return o1.compareToIgnoreCase(o2);
			} else {
				return Integer.compare(o1.length(), o2.length());
			}
		}

	}

	/*
	 * Set up the GUI
	 */
	public void start(Stage primaryStage) {
		/*
		 * Welcome text
		 */
		welcomeText = new Text("Welcome to the RNA-Sequencing Pipeline Generator. This program will "
				+ "assist in setting up and executing sequence alignment via Kallisto, which "
				+ "leads into differential expression analysis via DESeq2.");
		TextFlow textFlowWelcome = new TextFlow(welcomeText);
		setTextFlowLayout(textFlowWelcome);
		textFlowWelcome.setPadding(new Insets(DEFAULT_SPACING, DEFAULT_SPACING, 0, DEFAULT_SPACING));

		/*
		 * Creating the PWD and necessary file structure/files
		 */
		enterPWDText = new Text(
				"Enter the path of your present working directory. If it does not exist, it will be created. "
						+ "Your present working directory will be used to set up the necessary file structure. "
						+ "Pressing the Clear button will allow you to re-enter your present working directory, "
						+ "but you will have to perform all subsequent steps in the same order again.\n"
						+ "\nExample: /Users/userName/Desktop/RNA_Seq_Analysis\n"
						+ "Note: if using either Windows or Mac, enter the path with forward slashes as shown.");
		TextFlow textFlowPWD = new TextFlow(enterPWDText);
		setTextFlowLayout(textFlowPWD);

		pWDTextFieldInput = new TextField("Enter present working directory here");
		pWDTextFieldInput.setMinWidth(WIDTH_OF_CONTROLS - 220);
		createFileStructureButton = new Button("Create File Structure");
		createFileStructureButton.setOnAction(this::handleFileStructureCreation);
		clearPWDButton = new Button("Clear");
		clearPWDButton.setOnAction(this::handleFileStructureCreation);
		HBox pWDInputBox = new HBox(pWDTextFieldInput, createFileStructureButton, clearPWDButton);
		pWDInputBox.setSpacing(DEFAULT_SPACING);

		currentPWDText = new Text("");
		currentPWDText.setVisible(false);

		/*
		 * Enter the sample names and create directories for samples
		 */
		// Text Area for input
		sampleNamesInputTextArea = new TextArea();
		sampleNamesInputTextArea.setWrapText(false);
		sampleNamesInputTextArea.setEditable(true);
		sampleNamesInputTextArea.setMaxSize(WIDTH_OF_GUI / 3, (HEIGHT_OF_GUI / 4));
		// So that the user can paste from excel - replaces \r with \n for a
		// newline
		UnaryOperator<Change> filter = c -> {
			c.setText(c.getText().replaceAll("\r", "\n"));
			return c;
		};
		sampleNamesInputTextArea.setTextFormatter(new TextFormatter<>(filter));

		// Instructions
		sampleNamesInputInstructions = new Text();
		sampleNamesInputInstructions
				.setText("Check the execution log to make sure the file structure was set up properly.\n\n"
						+ "Please enter the accession number (if performing fastq dump, otherwise put \"none\"), "
						+ "sample names, and the condition for that sample, with a comma "
						+ "separating each component, one <accession#>,<sample name>,<condition> per line. "
						+ "Example: \"SRR1027605,C9_52i_run1,exp\". "
						+ "Sample names should not include any spaces or illegal "
						+ "characters for filenames. Press the button to continue. Directories for each sample will be "
						+ "made if necessary.");
		TextFlow textFlowSampleNamesInputInstructions = new TextFlow(sampleNamesInputInstructions);
		textFlowSampleNamesInputInstructions.setScaleShape(true);
		textFlowSampleNamesInputInstructions.setMaxSize(WIDTH_OF_GUI / 2, HEIGHT_OF_GUI / 8);

		// Enter Button
		enterSampleNamesButtonKallisto = new Button("Enter Sample Names");
		enterSampleNamesButtonKallisto.setOnAction(this::writeSamplesFiles);

		// Final HBox to add to mainBox
		VBox sampleNamesInputVBox = new VBox(textFlowSampleNamesInputInstructions, enterSampleNamesButtonKallisto);
		sampleNamesInputVBox.setAlignment(Pos.CENTER);
		sampleNamesInputVBox.setSpacing(DEFAULT_SPACING);
		sampleNamesInputHBox = new HBox(sampleNamesInputTextArea, sampleNamesInputVBox);
		setGUIVBoxHBoxLayout(sampleNamesInputHBox);

		/*
		 * Choose whether downloading from NCBI or not, and single-end or paired-end -
		 * For fastq dump
		 */
		// Instructions
		fastqDumpInstructions = new Text();
		fastqDumpInstructions.setText("Check the execution log to make sure the directories for each sample were "
				+ "created properly and that your samples_kallisto.txt file in the kallisto directory was filled with the "
				+ "correct sample names. Then make the proper choices below.\n\nWhen you click \"Generate fastq "
				+ "dump instructions\", instructions will be provided for performing a fastq dump (or moving "
				+ "your fastq files to the proper place if you don't need to download from NCBI). These "
				+ "instructions will be appended to your RNA_Seq_Pipeline.md file.");
		TextFlow textFlowFastqDumpInstructions = new TextFlow(fastqDumpInstructions);
		setTextFlowLayout(textFlowFastqDumpInstructions);

		// Downloading from NCBI Toggle Group
		ToggleGroup downloadNCBIToggle = new ToggleGroup();
		downloadNCBI = new RadioButton("I need to download my .fastq files from NCBI");
		downloadNCBI.setToggleGroup(downloadNCBIToggle);
		noDownloadNCBI = new RadioButton("I have my own .fastq files already");
		noDownloadNCBI.setToggleGroup(downloadNCBIToggle);
		VBox downloadNCBIVBox = new VBox(downloadNCBI, noDownloadNCBI);
		// Single-end or paired-end Toggle Group
		ToggleGroup singleOrPairedToggle = new ToggleGroup();
		singleEnd = new RadioButton("Single-End");
		singleEnd.setToggleGroup(singleOrPairedToggle);
		pairedEnd = new RadioButton("Paired-End");
		pairedEnd.setToggleGroup(singleOrPairedToggle);
		VBox singleOrPairedVBox = new VBox(singleEnd, pairedEnd);

		// HBox for RadioButtons
		HBox fastqDumpOptionsHBox = new HBox(downloadNCBIVBox, new Separator(Orientation.VERTICAL), singleOrPairedVBox);
		fastqDumpOptionsHBox.setAlignment(Pos.CENTER);
		fastqDumpOptionsHBox.setSpacing(DEFAULT_SPACING);

		// Button to generate fastq instructions
		generateFastQDumpButton = new Button("Generate fastq dump instructions");
		generateFastQDumpButton.setOnAction(this::writeFastqDump);

		// Final VBox to be added to mainBox
		fastqDumpOptionsVBox = new VBox(textFlowFastqDumpInstructions, fastqDumpOptionsHBox, generateFastQDumpButton);
		setGUIVBoxHBoxLayout(fastqDumpOptionsVBox);

		/*
		 * Write to the pipeline file the instructions for Kallisto alignment
		 */

		// Instructions
		kallistoInstructions = new Text(
				"After following instructions for the fastq dump, close the RNA_Seq_Pipeline.md "
						+ "file, then click the button below to generate instructions for Kallisto alignment, "
						+ "which will be appended to the RNA_Seq_Pipeline.md file.");
		TextFlow textFlowKallistoInstructions = new TextFlow(kallistoInstructions);
		setTextFlowLayout(textFlowKallistoInstructions);

		// Button to generate Kallisto instructions
		generateKallistoInstructionsButton = new Button("Generate Kallisto alignment instructions");
		generateKallistoInstructionsButton.setOnAction(this::writeKallistoInstructions);

		// Single end instructions
		singleEndInstructions = new Text();
		TextFlow singleEndInstructionsTextFlow = new TextFlow(singleEndInstructions);
		setTextFlowLayout(singleEndInstructionsTextFlow);

		// Single end Textfields for fragment length and std. dev. entry
		Text fragmentLengthText = new Text("Enter fragment length: ");
		fragmentLengthInput = new TextField("Ex: 50.0");
		Text standardDevText = new Text("Enter standard deviation: ");
		standardDevInput = new TextField("Ex: 0.0");
		singleEndInputHBox = new HBox(fragmentLengthText, fragmentLengthInput, standardDevText, standardDevInput);
		setGUIVBoxHBoxLayout(singleEndInputHBox);

		// Final VBox to be added to mainBox
		kallistoInstructionsVBox = new VBox(textFlowKallistoInstructions, singleEndInstructionsTextFlow,
				singleEndInputHBox, generateKallistoInstructionsButton);
		setGUIVBoxHBoxLayout(kallistoInstructionsVBox);

		/*
		 * Write to the pipeline file the instructions for DESeq2
		 */
		deseq2Instructions = new Text(
				"After following instructions for kallisto alignment, close the RNA_Seq_Pipeline.md "
						+ "file, then click the button below to generate instructions for differential expression analysis, "
						+ "which will be appended to the RNA_Seq_Pipeline.md file.");
		TextFlow textFlowDeseq2Instructions = new TextFlow(deseq2Instructions);
		setTextFlowLayout(textFlowDeseq2Instructions);

		generateDeseq2InstructionsButton = new Button("Generate DESeq2 instructions");
		generateDeseq2InstructionsButton.setOnAction(this::writeDESeq2Instructions);

		// Final VBox to be added to mainBox
		deseq2InstructionsVBox = new VBox(textFlowDeseq2Instructions, generateDeseq2InstructionsButton);
		setGUIVBoxHBoxLayout(deseq2InstructionsVBox);

		/*
		 * set the scene and stage
		 */
		mainBox = new VBox(textFlowWelcome, new Separator(), textFlowPWD, pWDInputBox, currentPWDText, new Separator(),
				sampleNamesInputHBox, new Separator(), fastqDumpOptionsVBox, new Separator(), kallistoInstructionsVBox,
				new Separator(), deseq2InstructionsVBox, new Separator());
		mainBox.setAlignment(Pos.TOP_CENTER);
		mainBox.setSpacing(10);

		mainScrollPane = new ScrollPane();
		mainScrollPane.setContent(mainBox);
		mainScrollPane.setHbarPolicy(ScrollBarPolicy.NEVER);
		mainScrollPane.setVbarPolicy(ScrollBarPolicy.ALWAYS);
		executionLog = new TextArea("Execution Log\n");
		executionLog.setEditable(false);
		executionLog.setWrapText(true);
		BorderPane scenePane = new BorderPane();
		scenePane.setCenter(mainScrollPane);
		scenePane.setBottom(executionLog);
		createFileStructureStage = new Scene(scenePane, WIDTH_OF_GUI, WIDTH_OF_GUI);
		primaryStage.setScene(createFileStructureStage);
		primaryStage.setResizable(false);
		primaryStage.setTitle("RNA Sequencing with Kallisto and DESeq2");
		primaryStage.show();
	}

	/*
	 * Start the application
	 */
	public static void main(String[] args) {
		launch(args);
	}

	/*
	 * Appends a string to a new line in the execution log along with the execution
	 * counter
	 */
	private void appendToExecutionLogNewLine(String str) {
		executionLog.appendText("\n" + executionCounter + ") " + str);
		executionCounter++;
	}

	/*
	 * Sets the layout of TextFlow objects to fit the GUI dimensions
	 */
	private void setTextFlowLayout(TextFlow textFlow) {
		textFlow.setScaleShape(true);
		textFlow.setPadding(new Insets(0, DEFAULT_SPACING, 0, DEFAULT_SPACING));
		textFlow.setMaxWidth(WIDTH_OF_CONTROLS);
	}

	/*
	 * Sets the layout of VBox's for positioning, spacing, and visibility
	 */
	private void setGUIVBoxHBoxLayout(VBox box) {
		box.setAlignment(Pos.CENTER);
		box.setSpacing(DEFAULT_SPACING);
		box.setVisible(false);
	}

	/*
	 * Sets the layout of HBox's for positioning, spacing, and visibility
	 */
	private void setGUIVBoxHBoxLayout(HBox box) {
		box.setAlignment(Pos.CENTER);
		box.setSpacing(DEFAULT_SPACING);
		box.setVisible(false);
	}

	/*
	 * Takes a String input as the file name and applies all flag substitutions to
	 * the text of that file, returning a String of the text with flag substitutions
	 */
	private String getInstructionsTextWithFlagSubstitutions(InputStream instructionsFile) {
		Alert supportingFileNotFoundError = new Alert(AlertType.ERROR);
		supportingFileNotFoundError.setContentText(
				"Supporting file for generating pipeline instructions not found. JAR file may be corrupt.");

		try (Scanner scanFile = new Scanner(instructionsFile)) {
			// \Z delimiter is for the end of the input, but for the final terminator -
			// reads the entire file
			scanFile.useDelimiter("\\Z");
			String flagsAppliedString = scanFile.next();

			for (FlagForCodeSubstitution flag : flagList) {
				flagsAppliedString = flagsAppliedString.replace(flag.getFlag(), flag.getValueToSubstitute());
			}

			return flagsAppliedString;
		}
	}

	/*
	 * Creates a new file or directory only if it does not already exist. Will
	 * recursively create directories if need be. Asks the user to create a file
	 * that cannot be created due to Security Exception. Returns a File object
	 * referencing new file created
	 */
	private File createNonexistentFile(String pathToFile, String description, boolean isDir) {

		File newFile = new File(pathToFile);

		// if the file does not exist, recursively create parent directories then create
		// the file
		// if the file cannot be created, prints a message about the error and how to
		// fix it
		if (!newFile.exists()) {
			boolean fileCreated = false;
			String errorMessage = null;
			// if newFile is a directory, it will recursively make the directory.
			// Otherwise newFile is a File and it will recursively create directories and
			// then the file
			try {
				if (isDir)
					fileCreated = newFile.mkdirs();
				else {
					newFile.getParentFile().mkdirs();
					fileCreated = newFile.createNewFile();
				}
			} catch (SecurityException se) {
				errorMessage = "Cannot create directories and/or file due to security exception - not accessible. "
						+ "You will need to change the working directory or create them yourself as this path: "
						+ newFile.getAbsolutePath();
			} catch (IOException ex) {
				errorMessage = "Cannot create directories and/or file due to Input/Output error. This may "
						+ "be due to a problem with your path input. Re-enter your path to try again.";
			} finally {
				if (errorMessage != null) {
					appendToExecutionLogNewLine(errorMessage);
				}
			}

			if (fileCreated) {
				appendToExecutionLogNewLine(description + " created:\n\t" + newFile.getAbsolutePath());
			}
		} else {
			appendToExecutionLogNewLine(description + " already exists:\n\t" + newFile.getAbsolutePath());
		}
		// returns a file that exists or will be created by user - will crash the
		// program if file is not created
		return newFile;
	}

	/*
	 * Creates the necessary file structure for analysis. Includes creation of
	 * necessary directories (based on the present working directory path provided)
	 * and text files
	 */
	private void handleFileStructureCreation(ActionEvent event) {

		// Clear Button
		if (event.getSource() == clearPWDButton) {
			pWDTextFieldInput.clear();

			// Enable Text Input and Create File Structures button
			pWDTextFieldInput.setDisable(false);
			createFileStructureButton.setDisable(false);
			appendToExecutionLogNewLine("\n\nClear button pressed: if entering a new working directory, restart "
					+ "from the beginning of the program.\n");
			executionCounter = 1;
		}

		// Create File Structure Button
		else if (event.getSource() == createFileStructureButton) {

			// Print beginning of creating file structure to execution log
			appendToExecutionLogNewLine("BEGIN CREATING FILE STRUCTURE");

			// Disable Text Input and Create File Structures button
			pWDTextFieldInput.setDisable(true);
			createFileStructureButton.setDisable(true);

			// Set next section visible
			sampleNamesInputHBox.setVisible(true);

			// Create PWD
			pWD = createNonexistentFile(pWDTextFieldInput.getText(), "Present working directory", true);
			currentPWDText.setText("Your present working directory is: " + pWD.getAbsolutePath());
			currentPWDText.setVisible(true);

			// Error message for later
			appendToRNASeqPipelineFileError = "Attempted to append text to RNA_Seq_Pipeline.md file, but file "
					+ "does not exist. Please create a text file at: " + pWD.getAbsolutePath()
					+ "/RNA_Seq_Pipeline.md\nYou will have to click Clear and then Create File Structure "
					+ "once more.";

			// Create pipelineFile
			pipelineFile = createNonexistentFile(pWD.getAbsolutePath() + "/RNA_Seq_Pipeline.md",
					"pipeline markdown file", false);
			try (PrintWriter printToPipelineFile = new PrintWriter(new FileWriter(pipelineFile))) {
				String timeStamp = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss").format(Calendar.getInstance().getTime());
				printToPipelineFile.println("Pipeline created: " + timeStamp + "\n");
			} catch (IOException ex) {
				appendToExecutionLogNewLine(appendToRNASeqPipelineFileError);
			}

			// Create kallisto directories and files
			kallisto = createNonexistentFile(pWD.getAbsolutePath() + "/kallisto", "kallisto directory", true);
			samplesKallisto = createNonexistentFile(kallisto.getAbsolutePath() + "/samples_kallisto.txt",
					"samples file for kallisto alignment", false);
			shellScriptsKallisto = createNonexistentFile(kallisto.getAbsolutePath() + "/shell_scripts",
					"kallisto shell scripts directory", true);
			samplesAccession = createNonexistentFile(kallisto.getAbsolutePath() + "/samples_Accession.txt",
					"samples file with accession numbers for fastq dump", false);

			// create DESeq2 directories and files
			deseq2 = createNonexistentFile(pWD.getAbsolutePath() + "/DESeq2", "DESeq2 directory", true);
			samplesDeseq2 = createNonexistentFile(deseq2.getAbsolutePath() + "/samples_deseq2.txt",
					"samples file for DESeq2 analysis", false);
			shellScriptsDeseq2 = createNonexistentFile(deseq2.getAbsolutePath() + "/shell_scripts",
					"DESeq2 shell scripts directory", true);
			deseq2_kallisto = createNonexistentFile(deseq2.getAbsolutePath() + "/kallisto",
					"kallisto subdirectory within DESeq2 directory", true);
			deseq2_results = createNonexistentFile(deseq2.getAbsolutePath() + "/results",
					"results subdirectory within DESeq2 directory", true);

			// Print end of creating file structure to execution log
			appendToExecutionLogNewLine("END CREATING FILE STRUCTURE\n");

			// Set up Flags for text substitution
			// Omit FRAGMENT_LENGTH_FLAG and STD_DEV_FLAG until they are set up in the
			// writeKallistoInstructions method
			PWD_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution("@@@PWD_ABSOLUTEPATH@@@", pWD.getAbsolutePath());
			KALLISTO_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution("@@@KALLISTO_ABSOLUTEPATH@@@",
					kallisto.getAbsolutePath());
			KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution(
					"@@@KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH@@@", shellScriptsKallisto.getAbsolutePath());
			KALLISTO_SAMPLES_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution("@@@KALLISTO_SAMPLES_ABSOLUTEPATH@@@",
					samplesKallisto.getAbsolutePath());
			ACCESSION_SAMPLES_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution("@@@ACCESSION_SAMPLES_ABSOLUTEPATH@@@",
					samplesAccession.getAbsolutePath());
			DESEQ2_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution("@@@DESEQ2_ABSOLUTEPATH@@@",
					deseq2.getAbsolutePath());
			DESEQ2_SHELLSCRIPTS_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution(
					"@@@DESEQ2_SHELLSCRIPTS_ABSOLUTEPATH@@@", shellScriptsDeseq2.getAbsolutePath());
			DESEQ2_KALLISTO_SUBDIR_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution(
					"@@@DESEQ2_KALLISTO_SUBDIR_ABSOLUTEPATH@@@", deseq2_kallisto.getAbsolutePath());
			DESEQ2_RESULTS_ABSOLUTEPATH_FLAG = new FlagForCodeSubstitution("@@@DESEQ2_RESULTS_ABSOLUTEPATH@@@",
					deseq2_results.getAbsolutePath());

			// Omit FRAGMENT_LENGTH_FLAG and STD_DEV_FLAG until they are set up in the
			// writeKallistoInstructions method
			Collections.addAll(flagList, PWD_ABSOLUTEPATH_FLAG, KALLISTO_ABSOLUTEPATH_FLAG,
					KALLISTO_SHELLSCRIPTS_ABSOLUTEPATH_FLAG, KALLISTO_SAMPLES_ABSOLUTEPATH_FLAG,
					ACCESSION_SAMPLES_ABSOLUTEPATH_FLAG, DESEQ2_ABSOLUTEPATH_FLAG,
					DESEQ2_SHELLSCRIPTS_ABSOLUTEPATH_FLAG, DESEQ2_KALLISTO_SUBDIR_ABSOLUTEPATH_FLAG,
					DESEQ2_RESULTS_ABSOLUTEPATH_FLAG);
		}
	}

	/*
	 * Creates samples files for both kallisto and DESeq2. samples_kallisto.txt
	 * includes only sample names while samples_deseq2_.txt includes sample names
	 * with conditions separated by a tab and includes headers for each column. Text
	 * from the user is read in and parsed by a comma delimiter, then printed to the
	 * respective text file
	 */
	private void writeSamplesFiles(ActionEvent event) {

		// Print beginning of creating samples file to execution log
		appendToExecutionLogNewLine("BEGIN CREATING SAMPLES FILES AND SAMPLE DIRECTORIES");
		try {
			PrintWriter printToSamplesKallisto = new PrintWriter(new FileWriter(samplesKallisto), false);
			PrintWriter printToSamplesDESeq2 = new PrintWriter(new FileWriter(samplesDeseq2), false);
			PrintWriter printToSamplesAccession = new PrintWriter(new FileWriter(samplesAccession), false);
			printToSamplesDESeq2.println("run\tcondition");
			String allSamples = sampleNamesInputTextArea.getText();
			if (!allSamples.equals("")) { // if the text area is not empty
				String[] sampleArray = allSamples.split("\\n");
				sampleListSorted = new ArrayList<>(Arrays.asList(sampleArray));
				Collections.sort(sampleListSorted, LENGTH_LEXICO_COMPARATOR);
				for (String sample : sampleListSorted) {
					String[] sampleAttributes = sample.split(",");
					String sampleAccession = sampleAttributes[0];
					String sampleName = sampleAttributes[1];
					String sampleCondition = sampleAttributes[2];
					printToSamplesKallisto.println(sampleName);
					createNonexistentFile(kallisto.getAbsolutePath() + "/" + sampleName, "kallisto sample directory",
							true);
					printToSamplesDESeq2.println(sampleName + "\t" + sampleCondition);
					createNonexistentFile(deseq2.getAbsolutePath() + "/kallisto/" + sampleName,
							"DESeq2 sample directory", true);
					printToSamplesAccession.println(sample);
				}
			}
			printToSamplesKallisto.close();
			printToSamplesDESeq2.close();
			printToSamplesAccession.close();
		} catch (IOException ex) {
			// File should be created in previous step
			appendToExecutionLogNewLine("Attempted to append text to samples_kallisto.txt file, but file "
					+ "does not exist. Please create a text file at: " + kallisto.getAbsolutePath()
					+ "/samples_kallisto.txt\nYou will have to click Enter Sample Names " + "once more.");
		} catch (IndexOutOfBoundsException ex) {
			Alert errorAlert = new Alert(AlertType.ERROR);
			errorAlert.setContentText("Sample names entered improperly. Enter the accession number (or none) "
					+ "followed by a comma, then the sample name, followed by a comma, "
					+ "then the condition for that sample (with no spaces in between). "
					+ "Example: SRR292241,sampleName,positive");
			errorAlert.showAndWait();
			appendToExecutionLogNewLine("SAMPLE FILES AND SAMPLE DIRECTORIES NOT CREATED\n");
			return;
		}
		// Print end of creating samples file to execution log
		appendToExecutionLogNewLine("END CREATING SAMPLES FILES AND SAMPLE DIRECTORIES\n");

		// Make next section visible
		fastqDumpOptionsVBox.setVisible(true);
		// scroll to the bottom of the main scroll pane
		mainScrollPane.setVvalue(1.0);
	}

	/*
	 * Appends instruction for fastq-dump to RNA_Seq_Pipeline.md file. This is
	 * based on whether or not the user needs to download from NCBI or not, and if
	 * the fastq files are single-end or paired-end. Also appends instructions for
	 * determining fragment length and standard deviation of reads if single-end.
	 */
	private void writeFastqDump(ActionEvent event) {

		try (PrintWriter printToPipelineFile = new PrintWriter(new FileWriter(pipelineFile, true))) {

			Alert buttonUnpressedAlert = new Alert(AlertType.ERROR);
			buttonUnpressedAlert.setContentText("You must make a selection for both sets of choices.");

			boolean downloadFromNCBI;

			if (downloadNCBI.isSelected()) {
				downloadFromNCBI = true;
			} else if (noDownloadNCBI.isSelected()) {
				downloadFromNCBI = false;
			} else { // neither button pressed
				buttonUnpressedAlert.showAndWait();
				return;
			}

			if (singleEnd.isSelected()) {
				singleTPairedF = true;
				singleEndInstructions.setText("Since your fastq files are single end, you will need to specify your "
						+ "estimated average fragment length and estimated standard deviation of fragment length. "
						+ "Instructions for doing so have been provided in the RNA_Seq_Pipeline.md file. Please "
						+ "follow these instructions and enter the values below before proceeding with Kallisto "
						+ "alignment.");
				singleEndInputHBox.setVisible(true);
			} else if (pairedEnd.isSelected()) {
				singleEndInputHBox.setVisible(false);
				singleTPairedF = false;
			} else {
				buttonUnpressedAlert.showAndWait();
				return;
			}

			// Instructions to move fastq files to proper directory
			if (!downloadFromNCBI) {
				printToPipelineFile.println(getInstructionsTextWithFlagSubstitutions(moveFastqFilesInstructionsFile));
			}
			// Downloading from NCBI
			else {
				// For both single and paired
				printToPipelineFile.println(getInstructionsTextWithFlagSubstitutions(fastqDumpInstructionsBeginFile));
				// single end fastq dump discludes --split-files parameter
				if (singleTPairedF) {
					printToPipelineFile
							.println(getInstructionsTextWithFlagSubstitutions(fastqDumpInstructionsScriptSingleFile));
				}
				// paired-end fastq dump
				else {
					printToPipelineFile
							.println(getInstructionsTextWithFlagSubstitutions(fastqDumpInstructionsScriptPairedFile));
				}
				// For both single and paired
				printToPipelineFile.println(getInstructionsTextWithFlagSubstitutions(fastqDumpInstructionsEndFile));
			}

			if (singleTPairedF) {
				printToPipelineFile
						.println(getInstructionsTextWithFlagSubstitutions(fragmentLengthStdDevInstructionsFile));
			}

			appendToExecutionLogNewLine("Instructions for obtaining or moving fastq files to proper location "
					+ "is provided in: " + pipelineFile.getAbsoluteFile() + "\nPlease open this file and "
					+ "follow instructions, " + "then close the file when you are finished.\n");
		} catch (IOException ex) {
			appendToExecutionLogNewLine(appendToRNASeqPipelineFileError);
		}

		// make the next section visible
		kallistoInstructionsVBox.setVisible(true);
		// scroll to the bottom of the main scroll pane
		mainScrollPane.setVvalue(1.0);

	}

	/*
	 * Appends instructions for kallisto alignment to the RNA_Seq_Pipeline.md file
	 * depending on if the fastq files are single or paired end. Also reads in
	 * fragment length and standard deviation values from user.
	 */
	private void writeKallistoInstructions(ActionEvent event) {
		// For single end
		if (singleTPairedF) {
			Alert invalidInputError = new Alert(AlertType.ERROR);
			invalidInputError.setContentText("Invalid input for fragment length or standard deviation. "
					+ "These values should be positive numbers. Enter the values and press the "
					+ "\"Generate Kallisto alignment instructions\" button again.");
			try {
				DecimalFormat doubleFormatter2Dec = new DecimalFormat("0.00");
				DecimalFormat doubleFormatter6Dec = new DecimalFormat("0.000000");
				fragmentLength = Double.parseDouble(fragmentLengthInput.getText());
				standardDev = Double.parseDouble(standardDevInput.getText());
				if (fragmentLength < 0.0 || standardDev < 0.0) {
					throw new IllegalArgumentException();
				}
				if (standardDev == 0.0) {
					standardDev = 0.000001; // required or kallisto alignment will crash
				}
				// Set up flags for fragment length and standard deviation now that they are
				// entered by user
				FRAGMENT_LENGTH_FLAG = new FlagForCodeSubstitution("@@@FRAGMENT_LENGTH@@@",
						doubleFormatter2Dec.format(fragmentLength));
				flagList.add(FRAGMENT_LENGTH_FLAG);
				STD_DEV_FLAG = new FlagForCodeSubstitution("@@@STD_DEV@@@", doubleFormatter6Dec.format(standardDev));
				flagList.add(STD_DEV_FLAG);
			} catch (IllegalArgumentException ex) {
				invalidInputError.showAndWait();
				return;
			}
		}

		try (PrintWriter printToPipelineFile = new PrintWriter(new FileWriter(pipelineFile, true))) {
			// For both single and paired
			printToPipelineFile.println(getInstructionsTextWithFlagSubstitutions(kallistoInstructionsBeginFile));
			// For single end
			if (singleTPairedF) {
				printToPipelineFile
						.println(getInstructionsTextWithFlagSubstitutions(kallistoInstructionsScriptSingleFile));
			}
			// For paired end
			else {
				printToPipelineFile
						.println(getInstructionsTextWithFlagSubstitutions(kallistoInstructionsScriptPairedFile));
			}
			// For both single and paired
			printToPipelineFile.println(getInstructionsTextWithFlagSubstitutions(kallistoInstructionsEndFile));

			appendToExecutionLogNewLine("Instructions for performing kallisto alignment " + "is provided in: "
					+ pipelineFile.getAbsoluteFile() + "\nPlease open this file and "
					+ "follow instructions, then close the file when you are finished.\n");

			// make next section visible
			deseq2InstructionsVBox.setVisible(true);
		} catch (IOException ex) {
			appendToExecutionLogNewLine(appendToRNASeqPipelineFileError);
		}
		
		// scroll to the bottom of the main scroll pane
		mainScrollPane.setVvalue(1.0);
	}

	/*
	 * Appends instructions for DESeq2 to the RNA_Seq_Pipeline.md file. Includes a
	 * script for moving abundance.h5 file to the correct directory as well as the
	 * code for input into RStudio
	 */
	private void writeDESeq2Instructions(ActionEvent event) {
		try (PrintWriter printToPipelineFile = new PrintWriter(new FileWriter(pipelineFile, true))) {
			printToPipelineFile.println(getInstructionsTextWithFlagSubstitutions(deseq2InstructionsFile));

			appendToExecutionLogNewLine("Instructions for performing differential expression analysis "
					+ "is provided in: " + pipelineFile.getAbsoluteFile() + "\nPlease open this file and "
					+ "follow instructions.\n");
		} catch (IOException ex) {
			appendToExecutionLogNewLine(appendToRNASeqPipelineFileError);
		}
	}
}

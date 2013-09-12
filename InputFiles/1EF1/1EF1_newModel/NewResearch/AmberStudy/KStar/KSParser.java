///////////////////////////////////////////////////////////////////////////////////////////////
// KSParser.java
//
//  Version:           0.3
//
//
// authors:
//    initials    name                 organization                email
//   ---------   -----------------    ------------------------    ----------------------------
//     RHL        Ryan Lilien          Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		  Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Written by Ryan Lilien (2001-2004) and Ivelin Georgiev (2004-2007)
 * 
 * The main class that handles KStar computation and related functions
 * 
 * This class has some remnants of a more general molecular modeling /
 *   manipulation program. Only those functions that might be of use
 *   have been left in. These remnants appear in classes like
 *   ScopeObject. Most of the small functions (readSystem, writeSystem,
 *   translateRelative, dispCenterOfMass, and so on) are not necessary
 *   for computing KStar.
 * The KStar functions include
 *   computePairEnergyMatrix -- computes a pairwise minimium energy matrix
 *   singleKStar -- computes a single partition function (bound or unbound)
 *   KSMaster -- performs a mutation search
 *   adjustNodes -- adjusts the work nodes during a mutation search
 * 
 */

/*
	This library is free software; you can redistribute it and/or
	modify it under the terms of the GNU Lesser General Public
	License as published by the Free Software Foundation; either
	version 2.1 of the License, or (at your option) any later version.
	This library is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
	Lesser General Public License for more details.
	
	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, write to the Free Software
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
	USA
	
	Contact Info:
		Bruce Donald
		Duke University
		Department of Computer Science
		Levine Science Research Center (LSRC)
		Durham
		NC 27708-0129 
		USA
		brd@cs.duke.edu
	
	If you use or publish any results derived from the use of this
	program please cite:
	Georgiev I, Lilien R, Donald B. "Improved Pruning Algorithms and
	Divide-and-Conquer Strategies for Dead-End Elimination, with Application 
	to Protein Design" Bioinformatics, 22(14): e174-e183, 2006.
	
	Copyright (C) 2006 Ivelin Georgiev, Ryan H. Lilien, and Bruce R. Donald
		
	<signature of Bruce Donald>, 23 Aug, 2006
	Bruce Donald, Professor of Computer Science
*/

import java.io.*;
import java.nio.channels.*;
import java.lang.Runtime;
import java.util.*;
import java.lang.Integer;
import java.math.*;
// import com.neva.*;   // Not compatible with linux

import mpi.MPI;
import mpi.MPIException;
import mpi.Status;

public class KSParser
{
		
		final int hashSize = 61; //the size of the Hastables
		
		 // setup the scope stack
		Hashtable curScope = new Hashtable(hashSize);
		 // setup hashtable (scope) with hashSize spaces (a prime number)
		Hashtable rParams = new Hashtable(hashSize);
		 // setup hashtable for runtime parameters with hashSize spaces
		boolean printSegID = false;
		
		boolean hElect = true; // should hydrogens be used in electrostatic energy calculations
		boolean hVDW = true; // should hydrogens be used in vdw energy calculations
		boolean hSteric = false; // should hydrogens be used in steric checks
		
		final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)
		
		final float stericThresh = 1.5f; // allowed overlap between the vdW radii of two atoms, steric clash if larger overlap
 	   
	// Config file name
		String cfgName = "KStar.cfg";
	
    // For soft vdw potential
		double softvdwMultiplier = 1.0;
	// For electrostatics
		boolean distDepDielect = true;
		double dielectConst = 1.0;
		
		//Determine if dihedral and solvation energies should be computed
		boolean doDihedE = false;
		boolean doSolvationE = false;
		double solvScale = 1.0;
		
		//the rotamer library
		RotamerLibrary rl = null;		
	
		int numAAallowed = -1; //number of allowed AA types
		String resAllowed[] = null; //the type of allowed AA
		int rotamerIndexOffset[] = null; //the rotamer index offset for each allowed AA
		int totalNumRotamers = -1; //the total number of rotamers for the Lovell rotamer library
				
		
		final int regTag = 1; //regular tag for MPI messages
		final int killTag = 2; //this MPI message tag tells the slave node to return
		int numProc = 0; //number of processors for MPI
		
		boolean mpiRun = false; //determines if this is an MPI run
	    
	    //the algorithm options that define what pruning criteria will be applied
		//	NOTE!!! These must be the same as in MutationManager.java
	    final int optSimple = 1;
	    final int optBounds = 2;
	    final int optSplit = 3;
	    final int optPairs = 4;
	    
	    //the assigned protein, ligand, and cofactor strand numbers
	    final int sysStrNum = 0; //the protein strand number is always 0
	    int ligStrNum = -1;
	    int cofStrNum = -1;

	
	// Checks if this is an MPI run and calss the respective functions
	public void checkMPI(String[] args) {
		
		if ((args.length>0)&&(args[0].equalsIgnoreCase("mpi"))) { //MPI run
			
			mpiRun = true;
			
			String tmp[] = new String[args.length-1]; //remove the mpi argument
			System.arraycopy(args, 1, tmp, 0, tmp.length);
			args = tmp;
			
			if ((args.length>0)&&(args[0].equalsIgnoreCase("-c"))){
				cfgName = args[1];
				String temp []= new String[args.length-2];
				System.arraycopy(args,2,temp,0,args.length-2);
				args = temp;
			}
			
			try{ handleDoMPI(args);} catch (Exception e){};
		}
		else {
			mpiRun = false;
			
			if ((args.length>0)&&(args[0].equalsIgnoreCase("-c"))){
				cfgName = args[1];
				String temp []= new String[args.length-2];
				System.arraycopy(args,2,temp,0,args.length-2);
				args = temp;
			}
			
			outputProgInfo(); //output program information
			setConfigPars(); //set the parameters from the configuration file
			
			parse(args); //parse the arguments
		}
	}
	
	//The main function which handles all commands
	public void parse(String[] args) {
		
		boolean done = false;
		boolean commandLineScript = false;
		boolean firstCommandLine = false;
		byte bytebuff[];
		String s = new String("");  // line being parsed
		
		if (args.length > 0) {
			commandLineScript = true;
			firstCommandLine = true;
		}
		
		// MAIN LOOP
		while (!done) {
			
			bytebuff = new byte[150];
			if (!commandLineScript) {
				System.out.print("> ");
				try {
					System.in.read(bytebuff);
				}
				catch ( Exception e ){
					System.out.println("ERROR: An error occurred while reading input");
					System.exit(0);
				}
				s = new String(bytebuff);  // create a string from bytebuff
			}
			else if (commandLineScript && !firstCommandLine) {
				// If you were running a command line script and the file is over then quit
				s = new String("quit");
			}
			else if (firstCommandLine) {
				s = new String("");
				for(int i=0;i<args.length;i++)
					s = s.concat(args[i] + " ");
				firstCommandLine = false;
			}			
			 
			s = s.trim();  // remove whitespace from beginning and end of line
				
			StringTokenizer st = new StringTokenizer(s," ;\t\n\r\f");
			String firstToken = new String("");
			if (st.hasMoreTokens())
				firstToken = st.nextToken();  // snag a copy of the first token
			 
			if (firstToken.length() > 1)
				if (firstToken.substring(0,2).equals("//"))
					continue;
			else if (firstToken.equalsIgnoreCase("quit") || firstToken.equalsIgnoreCase("q")){
				if (mpiRun){
					CommucObj cObj[] = new CommucObj[1];
					cObj[0] = new CommucObj();
					for (int curProc=1; curProc<numProc; curProc++){
						try { MPI.COMM_WORLD.Send(cObj, 0, 1, MPI.OBJECT, curProc, killTag);}
						catch (Exception e){};
					}
				}
					
				done = true;
			}
			
			else if (firstToken.equalsIgnoreCase("computePairEnergyMatricesSampling"))
				if (numTokens(s) >= 1) handleComputeAllPairwiseRotamerEnergiesMaster(s,curScope,null);
					else errorTooFewParams("computePairEnergyMatricesSampling");
			else if (firstToken.equalsIgnoreCase("doSinglePairE"))
				if (numTokens(s) >= 1) doSinglePairE(s,curScope,null);
					else errorTooFewParams("doSinglePairE");
			else if (firstToken.equalsIgnoreCase("doResEntropy"))
				if (numTokens(s) >= 1) handleDoResEntropy(s,curScope,null);
					else errorTooFewParams("doResEntropy");
			else if (firstToken.equalsIgnoreCase("compEref"))
				if (numTokens(s) >= 1) handleCompEref(s,curScope,null);
					else errorTooFewParams("compEref");
			else if (firstToken.equalsIgnoreCase("selectResidues"))
				if (numTokens(s) >= 1) selectResidues(s,curScope,null);
					else errorTooFewParams("selectResidues");
			
			else if (firstToken.equalsIgnoreCase("doDEE"))
				if (numTokens(s) >= 3) handleDoDEE(s,curScope);
					else errorTooFewParams("doDEE");
			else if (firstToken.equalsIgnoreCase("minimizeSequences"))
				if (numTokens(s) >= 7) handleMinDEEApplyRot(s,curScope);
					else errorTooFewParams("minimizeSequences");
			else if (firstToken.equalsIgnoreCase("generateRandConfs"))
				if (numTokens(s) >= 3) generateRandConfs(s,curScope);
					else errorTooFewParams("generateRandConfs");
			else if (firstToken.equalsIgnoreCase("fitEparams"))
				if (numTokens(s) >= 7) fitEparams(s,curScope);
					else errorTooFewParams("fitEparams");
			
			else if (firstToken.equalsIgnoreCase("singleKStar"))
				if (numTokens(s) >= 2) handleKSTest(s,curScope);
					else errorTooFewParams("singleKStar");
			else if (firstToken.equalsIgnoreCase("computeEnergyMol"))
				if (numTokens(s) >= 2) handleComputeEnergyMol(s,curScope);
					else errorTooFewParams("computeEnergyMol");
			else if (firstToken.equalsIgnoreCase("KSMaster"))
				if (numTokens(s) >= 3) handleKSMaster(s,curScope);
					else errorTooFewParams("KSMaster");
			else if (firstToken.equalsIgnoreCase("calculateAAVolumes"))
				if (numTokens(s) >= 1) handleCalculateAAVolumes(s,curScope);
					else errorTooFewParams("calculateAAVolumes");
			
			else if (firstToken.equalsIgnoreCase("genBackbones"))
				if (numTokens(s) >= 3) generateBackbones(s,curScope);
					else errorTooFewParams("genBackbones");	 
		} // End main "while (!done)" loop

	} // End parse function	
	
	//Displays the program version and citations
	public void outputProgInfo() {
		
		System.out.println();
		System.out.println("KStar Version 0.3");
		System.out.println("Duke");
		System.out.println("");
		System.out.println("Lilien R, Stevens B, Anderson A, Donald B. \"A Novel Ensemble-Based Scoring");
		System.out.println(" and Search Algorithm for Protein Redesign, and its Application to Modify");
		System.out.println(" the Substrate Specificity of the Gramicidin Synthetase A Phenylalanine");
		System.out.println(" Adenylation Enzyme\" Proceedings of the Eighth Annual International");
		System.out.println(" Conference on Research in Computational Molecular Biology (RECOMB),");
		System.out.println(" San Diego (March 27-31, 2004) pp. 46-57.");
		System.out.println("");
		System.out.println("Georgiev I, Lilien R, Donald B. \"A Novel Minimized");
		System.out.println(" Dead-End Elimination Criterion and Its Application to Protein Redesign");
		System.out.println(" in a Hybrid Scoring and Search Algorithm for Computing Partition Functions");
		System.out.println(" over Molecular Ensembles\" Proceedings of the Tenth Annual International");
		System.out.println(" Conference on Research in Computational Molecular Biology (RECOMB),");
		System.out.println(" Venice, Italy (March, 2006) pp. 530-545.");
		System.out.println("");
		System.out.println("Georgiev I, Lilien R, Donald B. \"Improved Pruning Algorithms and");
		System.out.println(" Divide-and-Conquer Strategies for Dead-End Elimination, with Application");
		System.out.println(" to Protein Design\" Bioinformatics, 22(14): e174-e183, 2006.");
		System.out.println("");
	}
	
	//Sets the parameters from the configuration file
	public void setConfigPars() {
		
		readConfigFile(rParams,cfgName); //read configuration file
		
		hElect = (new Boolean((String)rParams.get("HELECT"))).booleanValue();
		hVDW = (new Boolean((String)rParams.get("HVDW"))).booleanValue();
		hSteric = (new Boolean((String)rParams.get("HSTERIC"))).booleanValue();
		distDepDielect = (new Boolean((String)rParams.get("DISTDEPDIELECT"))).booleanValue();
		dielectConst = (new Double((String)rParams.get("DIELECTCONST"))).doubleValue();
		doDihedE = (new Boolean((String)rParams.get("DODIHEDE"))).booleanValue();
		doSolvationE = (new Boolean((String)rParams.get("DOSOLVATIONE"))).booleanValue();
		solvScale = (new Double((String)rParams.get("SOLVSCALE"))).doubleValue();
		softvdwMultiplier = (new Double((String)rParams.get("VDWMULT"))).doubleValue();
		
		rl = new RotamerLibrary((String)rParams.get("ROTFILE"),(String)rParams.get("VOLFILE"));
		numAAallowed = rl.getNumAAallowed();
		resAllowed = rl.getAAtypesAllowed();
		rotamerIndexOffset = rl.getRotamerIndexOffset();
		totalNumRotamers = rl.getTotalNumRotamers();
	}
	
	
	// Reads the configuration file containing runtime options
	// rp is an initialized empty Hashtable
	private void readConfigFile(Hashtable rp, String fname){
	
		BufferedReader bufread = null;
		String curLine = null;
		String param = null;
		String value = null;
		boolean done = false;
		
		// First attempt to open and read the config file
		try{
			FileInputStream is = new FileInputStream(fname);
			bufread = new BufferedReader(new InputStreamReader(is));

			curLine = bufread.readLine();

			while (curLine != null){
				done = false;
				while (!done) {
					if (curLine.charAt(0) == '%'){
						curLine = bufread.readLine();
					} else {
						done = true;
					}
					if (curLine == null)
						done = true;
				}
				if (curLine != null){
					param = getToken(curLine,1);
					value = curLine.substring(param.length()+1);
					rp.put(param.toUpperCase(),value);
					curLine = bufread.readLine();
				}
			}
			bufread.close();
		}
		catch(Exception ex)
		{
			System.out.println("ERROR: An error occurred reading configuration file "+fname);
			System.exit(0);
		}

	}

	/******************************/
	// Prints an error message because too few parameters were in line
	private void errorTooFewParams(String s) {
	
		System.out.println("ERROR: Too few parameters for " + s);
	}


	/******************************/
	// This function true if the variable named s is in the current scope
	private boolean inScope(String s, Hashtable curScope) {
	
		if (curScope.get(s) != null)
			return true;
		return false;
	}
		
	/******************************/
	// This function parses s into a double, either it is a constant
	//  or it is a double variable
	private double parseToDouble(String s, Hashtable curScope) {
	
		Double tmpDouble;
		double d = 0.0;
		ScopeObject so = null;
		
		if (inScope(s,curScope)) {
			so = (ScopeObject) curScope.get(s);
			if (so == null) {
				System.out.println("ERROR: undefined variable " + getToken(s,1));
				return 0.0;
			}
			if (so.objType != 12) {
				System.out.println("ERROR: variable " + s + " is not of type double");
				return 0.0;
			}
			tmpDouble = (Double) so.data;	
			d = tmpDouble.doubleValue();
		}
		else
			d = new Double(s).doubleValue();
	
		return d;
	}

	
	/******************************/
	// This function parses s into a float, either it is a constant
	//  or it is a double variable (floats are stored as double variables)
	private float parseToFloat(String s, Hashtable curScope) {
	
		Double tmpDouble;
		float f = 0.0f;
		ScopeObject so = null;
		
		if (inScope(s,curScope)) {
			so = (ScopeObject) curScope.get(s);
			if (so == null) {
				System.out.println("ERROR: undefined variable " + getToken(s,1));
				return 0.0f;
			}
			if (so.objType != 12) {
				System.out.println("ERROR: variable " + s + " is not of type double");
				return 0.0f;
			}
			tmpDouble = (Double) so.data;	
			f = tmpDouble.floatValue();
		}
		else
			f = new Float(s).floatValue();
	
		return f;
	}


	/******************************/
	// This function parses s into an integer, either it is a constant
	//  or it is an integer variable
	private int parseToInt(String s, Hashtable curScope) {
	
		Integer tmpInt;
		int d = 0;
		ScopeObject so = null;
		
		if (inScope(s,curScope)) {
			so = (ScopeObject) curScope.get(s);
			if (so == null) {
				System.out.println("ERROR: undefined variable " + getToken(s,1));
				return 0;
			}
			if (so.objType != 11) {
				System.out.println("ERROR: variable " + s + " is not of type integer");
				return 0;
			}
			tmpInt = (Integer) so.data;	
			d = tmpInt.intValue();
		}
		else
			d = new Integer(s).intValue();
	
		return d;
	}

	
	/******************************/
	// This function parses s into an boolean
	private boolean parseToBoolean(String s, Hashtable curScope) {
	
		Boolean tmpB = new Boolean(s);
		return tmpB.booleanValue();
	}

	
	/******************************/
	// This function returns the number of tokens in string s
	private int numTokens(String s) {
		
		int curNum = 0;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (st.hasMoreTokens()) {
			curNum++;
		  st.nextToken();
		}
		return(curNum);
	}


	/******************************/
	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				// System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
		
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

	// Computes the factorial of input n
	public BigInteger factorial(int n){
	
		if (n==0)
			return BigInteger.valueOf(1);
	
		return (factorial(n-1).multiply(BigInteger.valueOf(n)));
	}

	public void saveMolecule(Molecule m, String fname, float energy){
	
		try{
			FileOutputStream fileOutputStream = new FileOutputStream(fname);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( 
	fileOutputStream);
			PrintStream printStream = new PrintStream(bufferedOutputStream);
			Hashtable params = new Hashtable(7);
			params.put("printSegID",new Boolean(printSegID));
			params.put("comment","");
			params.put("energy", energy);
			params.put("showConnect",new Boolean(false));					
			new SaveMolecule(m, printStream, params); 
			printStream.close();
		}
		catch (IOException e) {
			System.out.println("ERROR: An io exception occurred while writing file");
			System.exit(0);
		}
		catch ( Exception e ){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing file");
			System.exit(0);
		}
	}

	
	// Uses the VolModule class to calculate the volume of each rotamer
	//  of each amino acid
	public void handleCalculateAAVolumes(String s, Hashtable curScope) {

		Molecule m = new Molecule();

		Amber96PolyPeptideResidue ppr = new Amber96PolyPeptideResidue();
		StrandRotamers LR = null;

		Residue res = ppr.getResidue("ala");
		//res.fullName = "ALA  ";
		m.addResidue(0,res);
		VolModule sm = new VolModule(m);
		sm.setResidueTreatment(0,1);
		
		LR = new StrandRotamers(rl,m.strand[sysStrNum]);		

		PrintStream printStream = setupOutputFile("aavolumes2.dat");

		String aanames[] = rl.getAAtypesAllowed();
		int numAAs = rl.getNumAAallowed();
		
		for(int i=0;i<numAAs;i++){
			LR.changeResidueType(m,0,aanames[i],true);
			printStream.print(aanames[i] + " ");
			System.out.println(aanames[i] + " ");
			if(rl.getNumRotamers(aanames[i])==0){		// this is the ala case
				float vol = sm.getMoleculeVolume(0.25f,0.0f);
				printStream.print(vol + " ");
				System.out.println(vol + " ");
			}			
			for(int j=0;j<rl.getNumRotamers(aanames[i]);j++){
				LR.applyRotamer(m,0,j);
				float vol = sm.getMoleculeVolume(0.25f,0.0f);
				printStream.print(vol + " ");
				System.out.println(vol + " ");
			}
			printStream.println();
		}
		printStream.close();
		
	}

	// This function generates all possible combinations of n choose m
	public void generateCombinations(int residueMutatable[][], int n, int m) {
	
		int curIndex[] = new int[1];
		int curComb[] = new int[n];
		curIndex[0] = 0;
		generateCombHelper(0,n,curIndex,residueMutatable,curComb,0,m);
	}
	private void generateCombHelper(int depth, int maxDepth, int curIndex[], int
		residueMutatable[][], int curComb[], int numUsed, int maxToUse){
		
		if (depth >= maxDepth){
			if (numUsed == maxToUse) {
				for (int i=0; i<maxDepth; i++) {
					residueMutatable[curIndex[0]][i] = curComb[i];
				}
				curIndex[0]++;
			}
			return;
		}
		
		curComb[depth] = 0;
		generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed,maxToUse);

		if (numUsed < maxToUse) {
			curComb[depth] = 1;
			generateCombHelper(depth+1,maxDepth,curIndex,residueMutatable,curComb,numUsed+1,maxToUse);
		}
	}
	// end combination code
	
	
	//Sets up the molecule system and returns the number of ligand rotamers
	private int setupMolSystem(Molecule m, Hashtable sParams, boolean ligPresent, String ligType){
		
		ligStrNum = -1;
		cofStrNum = -1;

		try{
			FileInputStream is = new FileInputStream((String)sParams.get("PDBNAME"));
			new PDBChemModel(m, is);
		}
		catch (Exception e){
			System.out.println("WARNING: An error occurred while reading file");
			System.out.println(e);
			return -1;
		}

		//Get the ligand (if present and if will be used) and cofactor (if present)
		int ligNum = (new Integer((String)sParams.get("LIGNUM"))).intValue();
		
		int numCofactorRes = (new Integer((String)sParams.get("COFACTORNUM"))).intValue();		
		String cofMapString = (String)sParams.get("COFMAP");
		int cofactorRes[] = new int[numCofactorRes];
		for(int i=0;i<numCofactorRes;i++)
			cofactorRes[i] = (new Integer(getToken(cofMapString,i+1))).intValue();
		
		Residue cof[] = null;
		if (numCofactorRes > 0) {  // pull out the cofactor residues
			cof = new Residue[numCofactorRes];
			for (int i=(numCofactorRes-1); i>=0; i--){ //backwards order, since we need to delete each of these residues from the molecule once they are extracted
				cof[i] = m.residue[cofactorRes[i]];
				m.deleteResidue(cofactorRes[i]);
				cof[i].renumberResidue();
			}
		}
		
		int strNum = sysStrNum + 1; //the current strand number; 0 is reserved for the protein strand, the ligand strand is 1 (if present)
		
		if (ligNum>=0){ //with ligand in PDB
			Residue lig = m.residue[ligNum]; // pull out the ligand
			m.deleteResidue(ligNum);
			if (ligPresent) { //ligand will be used in design
				lig.renumberResidue();
				m.addStrand("L");
				ligStrNum = strNum;
				m.addResidue(ligStrNum,lig);
				strNum++;
			}
		}
		if (numCofactorRes > 0){ //with cofactor
			m.addStrand("M");
			cofStrNum = strNum;
			for (int i=0; i<numCofactorRes; i++)
				m.addResidue(cofStrNum,cof[i]);
			strNum++;
		}

		m.strand[sysStrNum].isProtein = true;		// main active-site
		if (ligPresent)
			m.strand[ligStrNum].isProtein = true;	// ligand
		if (numCofactorRes > 0)
			m.strand[cofStrNum].isProtein = false;  // cofactor

		// Setup ligand rotamers and ligand type (if ligand exists)
		StrandRotamers ligLR = null;
		int numLigRotamers = 0;
		if ((ligNum>=0)&&(ligPresent)) {
			ligLR = new StrandRotamers(rl,m.strand[ligStrNum]); 

			// change ligand to the specified residue type
			if (!rl.getThreeLetAACode(ligLR.getCurRotType(0)).equalsIgnoreCase(ligType)) //not the same ligand type
				ligLR.changeResidueType(m,0,ligType,true);
				
			numLigRotamers = rl.getNumRotamers(ligLR.getCurRotType(0));
			if (numLigRotamers == 0)
				numLigRotamers = 1;
		}
		return numLigRotamers;
	}

	// Computes the bound or unbound partition function and
	//  can compute the energyminimized structure for a specified
	//  rotameric conformation
	public void handleKSTest(String s, Hashtable curScope) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Run name (for output files)
		// 3: Pairwise min energy matrix filename 
		//    (excluding the .dat suffix)
		//    (ie. AlaMatrix -> AlaMatrix.dat)
		// 4: Pairwise max energy matrix filename (write anything if not using MinDEE)
		// 5: Perform energy minimization (boolean)
		// 6: nDEE initial Ew (float)
		// 7: MinBounds cutoff energy for pruning (float)
		// 8: Ligand (boolean), is true if present
		// 9: Amino acid type for ligand (if ligand is absent, write none or anything)
		// 10...(n+9): Amino acid types of the n active site residues (1-letter codes)

		int numInAS = -1;

		// The system paramaters
		Hashtable sParams = new Hashtable(hashSize);

		// Read System Parameters
		readConfigFile(sParams,getToken(s,2));

		// Pull search parameters
		numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		String runName = getToken(s,3);
		String eMatrixNameMin = getToken(s,4);
		boolean doMinimize = parseToBoolean(getToken(s,6),curScope);
		boolean minimizeBB = parseToBoolean(getToken(s,7),curScope);
		if (minimizeBB)
			doMinimize = true;
		String eMatrixNameMax = null;
		if (doMinimize)
			eMatrixNameMax = getToken(s,5);
		float initEw = parseToFloat(getToken(s,8),curScope);
		float pruningE = parseToFloat(getToken(s,9),curScope);
		boolean ligPresent = parseToBoolean(getToken(s,10),curScope);
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,11);
		
		boolean scaleInt = parseToBoolean(getToken(s,12),curScope);
		float maxIntScale = parseToFloat(getToken(s,13),curScope);
		
		boolean saveConfs = parseToBoolean(getToken(s,14),curScope);
		String fName = "pdbs/"+getToken(s,15);
		
		double stericE = parseToFloat(getToken(s,16),curScope);
		
		boolean repeatSearch = true;
		
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			System.out.print(" "+residueMap[i]);
		}
		System.out.println();
		
		if (ligPresent)
			getLigPartFn(m,numInAS,eMatrixNameMin+".dat"); //compute the ligand partition function
		
		float epsilon = (new Float((String)sParams.get("EPSILON"))).floatValue();
	
		RotamerSearch rs = null;
		if (ligPresent)
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,true,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl);
		else
			rs = new RotamerSearch(m,sysStrNum,hElect,hVDW,hSteric,true,true,epsilon,stericThresh,distDepDielect,dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl);

		// Define the current amino-acid sequence
		String curSeq[] = new String[numInAS];
		for(int i=0;i<numInAS;i++){
			curSeq[i] = rl.getThreeLetAACode(getToken(s,i+17).charAt(0));
		}

		// Check for HIS
		for(int i=0;i<numInAS;i++) {
			if( (curSeq[i].equalsIgnoreCase("HIS")) || (curSeq[i].equalsIgnoreCase("HID")) || (curSeq[i].equalsIgnoreCase("HIE")) ){
				System.out.println(curSeq[i] + " changed to HIP");
				curSeq[i] = new String("HIP");
			}
		}
		
		System.out.print("Current Sequence:");
		for(int w=0;w<numInAS;w++)
			System.out.print(" "+curSeq[w]);
		System.out.println();

		System.out.println("Beginning setAllowables");
		for(int i=0;i<numInAS;i++){
			rs.setAllowable(residueMap[i],curSeq[i]);
		}

		System.out.print("Loading precomputed min energy matrix...");
		rs.loadPairwiseEnergyMatrices(new String(eMatrixNameMin+".dat"));
		System.out.println("done");		
		
		if (doMinimize){
			System.out.print("Loading precomputed max energy matrix...");
			rs.loadPairwiseEnergyMatricesMax(new String(eMatrixNameMax+".dat"));
			System.out.println("done");
		}
		
		System.out.println("Before start");		
			
		boolean prunedRotAtRes[] = new boolean[rs.arpMatrix.length];
		for (int i=0; i<prunedRotAtRes.length; i++)
			prunedRotAtRes[i] = false;
		
		//Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
		prunedRotAtRes = rs.DoPruneStericTemplate(numInAS, totalNumRotamers, numLigRotamers, 
				residueMap, rotamerIndexOffset, prunedRotAtRes, stericE);
		
		if (doMinimize) //precompute the interval terms in the MinDEE criterion
			rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, 
					rotamerIndexOffset, prunedRotAtRes, scaleInt, maxIntScale);
		
		prunedRotAtRes = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
				rotamerIndexOffset, initEw, prunedRotAtRes, doMinimize, false, minimizeBB);
		
		//Prune with MinBounds
		prunedRotAtRes = rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
				residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, false);
		
		//Compute the Ec value and prunedIsSteric[]
		rs.DoMinBounds(numInAS,totalNumRotamers,numLigRotamers,
				residueMap,rotamerIndexOffset,pruningE,prunedRotAtRes,initEw, false, true);
		
		//Do the rotamer search
		rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
				residueMap,false,(new BigDecimal("0.0")),null,minimizeBB,saveConfs,fName);
		
		if ((repeatSearch)&&(rs.repeatSearch)){ //the desired accuracy was not achieved, so repeat the search: the setup is already done
			
			System.out.println();
			System.out.println("Repeating search..");
			rs.repeatSearch = false; //reset the flag
			rs.slaveDoRotamerSearch(true,doMinimize,numInAS,numAAallowed,totalNumRotamers,rotamerIndexOffset,resAllowed,
					residueMap,false,(new BigDecimal("0.0")),null,minimizeBB,saveConfs,fName);
		}
	}
	
	// Finds the energy for a given input system (a molecule with specified flexible residues)
	public void handleComputeEnergyMol(String s, Hashtable curScope) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Ligand (boolean), is true if present
		// 3: Amino acid type for ligand (if ligand is absent, write none or anything)

		// The system paramaters
		Hashtable sParams = new Hashtable(hashSize);

		// Read System Parameters
		readConfigFile(sParams,getToken(s,2));

		// Pull search parameters
		int numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		int ligNum = (new Integer((String)sParams.get("LIGNUM"))).intValue();
		boolean ligPresent = parseToBoolean(getToken(s,3),curScope);
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,4);
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]);
		}
		System.out.println();
		
		System.out.println("Starting energy computation");
		Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
		
		StrandRotamers sysLR = null;
		StrandRotamers ligLR2 = null;
		
		sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
		if (ligPresent)
			ligLR2 = new StrandRotamers(rl,m.strand[ligStrNum]);

		for(int i=0;i<numInAS;i++){
			sysLR.setAllowable(residueMap[i],resDefault[i]);
			m.residue[residueMap[i]].flexible = true;
		}
		if (ligPresent){
			ligLR2.setAllowable(0,ligType);
			m.strand[ligStrNum].residue[0].flexible = true;
		}
		a96ff.calculateTypesWithTemplates();
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		if (ligPresent)
			a96ff.setLigandNum(ligNum);
		
		double energy[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1);
		System.out.println("System energy: " + energy[0]+" (elect: "+energy[1]+" vdW: "+energy[2]+" solvation: "+energy[3]+")");
	}

	
	// This function waits the specified number of milliseconds before returning
	//  for some reason I can't get the generic wait() method to work
	public void spinwait(long timetowait) {
		
		long startTime = System.currentTimeMillis();
		while ((System.currentTimeMillis() - startTime) < timetowait);

		return;
	}


	// Mutation search Master function
	public void handleKSMaster(String s, Hashtable curScope) {
	
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)
		
		// The system and mutation search paramaters
		Hashtable sParams = new Hashtable(hashSize);

		// Read System parameters
		readConfigFile(sParams,getToken(s,2));
		// Read Mutation search parameters
		readConfigFile(sParams,getToken(s,3));
		
		int numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		int numResAllowed = 19;
		int numTotalRotamers = 152;
		int numMutations = (new Integer((String)sParams.get("NUMMUTATIONS"))).intValue();
		int ligNum = (new Integer((String)sParams.get("LIGNUM"))).intValue();
		String runName = (String)sParams.get("RUNNAME");
		String mutFileName = (String)sParams.get("MUTFILENAME");
		String eMatrixNameMin = (String)sParams.get("MINENERGYMATRIXNAME");
		String eMatrixNameMax = (String)sParams.get("MAXENERGYMATRIXNAME");
		boolean minimizeBB = (new Boolean((String)sParams.get("MINIMIZEBB"))).booleanValue();
		boolean repeatSearch = (new Boolean((String)sParams.get("REPEATSEARCH"))).booleanValue();
		int algOption = (new Integer((String)sParams.get("ALGOPTION"))).intValue();
		int numSplits = (new Integer((String)sParams.get("NUMSPLITS"))).intValue();
		boolean approxMinGMEC = (new Boolean((String)sParams.get("APPROXMINGMEC"))).booleanValue();
		float lambda = (new Float((String)sParams.get("LAMBDA"))).floatValue();
		boolean scaleInt = (new Boolean((String)sParams.get("SCALEINT"))).booleanValue();
		float maxIntScale = (new Float((String)sParams.get("MAXINTSCALE"))).floatValue();
		float initEw = (new Float((String)sParams.get("INITEW"))).floatValue();
		double pruningE = (new Double((String)sParams.get("PRUNINGE"))).doubleValue();
		boolean ligPresent = (new Boolean((String)sParams.get("LIGPRESENT"))).booleanValue();
		String ligType = (String)sParams.get("LIGTYPE");
		float targetVol = (new Float((String)sParams.get("TARGETVOLUME"))).floatValue();
		float volWindow = (new Float((String)sParams.get("VOLUMEWINDOW"))).floatValue();
		boolean resumeSearch = (new Boolean((String)sParams.get("RESUMESEARCH"))).booleanValue();
		String resumeFilename = (String)sParams.get("RESUMEFILENAME");
		double gamma = (new Double((String)sParams.get("GAMMA"))).doubleValue();
		float epsilon = (new Float((String)sParams.get("EPSILON"))).floatValue();
		
		if (!mpiRun){
			System.out.println("ERROR: Distributed computation requires MPI");
			System.exit(1);
		}
		
		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Min Energy Matrix: "+eMatrixNameMin);
		System.out.println("Precomputed Max Energy Matrix: "+eMatrixNameMax);
		System.out.println("Ligand Type: "+ligType);
		System.out.println("Volume Center: "+targetVol);
		System.out.println("Volume Window Size: "+volWindow);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);
		
		if(resumeSearch) {
			System.out.println("** Resuming Search **");
			System.out.println("     resuming from file: "+resumeFilename);
		}
		
		// Create the mutation list with estimated energies
		OneMutation mutArray[] = new OneMutation[200000];
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		// Generate all combinations (include (n choose m), (n choose m-1), ... , (n choose 1), and (n choose 0) )
		int numCombAll = 0;
		for (int i=numMutations; i>=0; i--)
			numCombAll += factorial(numInAS).divide(factorial(numInAS-i).multiply(factorial(i))).intValue();
		int residueMutatableAll[][] = new int[numCombAll][numInAS];
		int curInd = 0;
		for (int i=numMutations; i>=0; i--){
			int numCombCur = factorial(numInAS).divide(factorial(numInAS-i).multiply(factorial(i))).intValue();
			int residueMutatableCur[][] = new int[numCombCur][numInAS];
			generateCombinations(residueMutatableCur,numInAS,i);
			for (int j=0; j<numCombCur; j++){
				residueMutatableAll[curInd] = residueMutatableCur[j];
				curInd++;
			}
		}
		
		// At this point each row of residueMutatble is a 0/1 array, 1 indicates
		//  that that residues can mutate

		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+resDefault[i]+")");
		}
		System.out.println();

		//Load mutation list for distribution
		mutArray = handleHybridKSLoadMutList(mutArray, mutFileName, numInAS, m, numCombAll, residueMutatableAll,
				sParams, residueMap, resDefault, resAllowed, numResAllowed, targetVol, volWindow);

		BigDecimal bestScore = new BigDecimal("0.0"); //for the resume results		
		// If doing a resume, read the initial results into a bunch of OneMutations
		if (resumeSearch) {
		
			bestScore = new BigDecimal("0.0"); //high scores are better
			
			OneMutation resumeResults[] = new OneMutation[mutArray.length];
			for(int q=0;q<mutArray.length;q++)
				resumeResults[q] = new OneMutation();
			resumeResults = readResumeFile(resumeResults,resumeFilename,numInAS,false,false);
			System.out.println("Read "+resumeResults.length+" completed mutations");
			
			// Now filter removed mutations (already computed results
			//  are NOT written to file since you already have them)
			// We do need to maintain the best score
			int newIndex = 0;
			OneMutation newArray2[] = new OneMutation[mutArray.length];
			for(int q=0;q<mutArray.length;q++) {
				int w = findMutationIndex(resumeResults,mutArray[q].resTypes);
				if (w>=0)
					bestScore = bestScore.max(resumeResults[w].score); //higher scores are better for Hybrid MinDEE-K*
				else 
					newArray2[newIndex++] = mutArray[q];
			}
			mutArray = new OneMutation[newIndex];
			System.arraycopy(newArray2,0,mutArray,0,newIndex);
			System.out.println("Length of mutArray after removing already computed mutations: "+mutArray.length);
		}
		
		BigDecimal q_L = getLigPartFn(m,numInAS,eMatrixNameMin+".dat");
	
		MutationManager mutMan = new MutationManager(runName,mutArray,false);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigNum(ligNum);
		mutMan.setLigType(ligType);
		mutMan.setLigPresent(ligPresent);
		mutMan.setNumMutations(numMutations);
		mutMan.setarpFilenameMin(eMatrixNameMin+".dat");
		mutMan.setarpFilenameMax(eMatrixNameMax+".dat");
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setRepeatSearch(repeatSearch);
		mutMan.setAlgOption(algOption);
		mutMan.setNumSplits(numSplits);
		mutMan.setInitEw(initEw);
		mutMan.setGamma(gamma);
		mutMan.setEpsilon(epsilon);
		mutMan.setApproxMinGMEC(approxMinGMEC);
		mutMan.setLambda(lambda);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.numTotalRotamers(numTotalRotamers);
		mutMan.numResAllowed(numResAllowed);
		mutMan.setIsLigAA(true);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(true);
		mutMan.setCalculateVolumes(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setLigPartFn(q_L);
		
		if (resumeSearch)
			mutMan.setBestScore(bestScore);	// Set the current best score from the partial results
		else
			mutMan.setBestScore(new BigDecimal("0.0")); //the initial best score is 0.0
		
		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
	}
	
	//Computes the partition function for the ligand using the rotamers from the rotamer library
	private BigDecimal getLigPartFn(Molecule m, int numInAS, String eMatrixNameMin){
		
		float minMatrix[][] = (float [][])readObject(eMatrixNameMin);
		
		int numRot = rl.getNumRotamers(m.strand[ligStrNum].residue[0].name);
		if (numRot==0) //ALA or GLY
			numRot = 1;
		
		int initInd = 1+numInAS*totalNumRotamers;
		
		ExpFunction ef = new ExpFunction();
		
		BigDecimal q_L = new BigDecimal("0.0");
		
		for (int i=0; i<numRot; i++)
			q_L = q_L.add(ef.exp(-minMatrix[initInd+i][0]/constRT));
		
		System.out.println("Ligand partition function (double): "+q_L.doubleValue());
		
		return q_L;
	}
	
	//Loads the mutation sequence list for Hybrid MinDEE-K*; computes a list if one cannot be loaded
	private OneMutation[] handleHybridKSLoadMutList (OneMutation mutArray[], String mutFileName, int numInAS,
			Molecule m, int numComb, int residueMutatable[][], Hashtable sParams,int residueMap[], String resDefault[],
			String resAllowed[], int numResAllowed, float targetVol,float volWindow){
		
		// Look for previous mutation file
		System.out.println();
		System.out.print("Looking for mutation list file ");
		mutArray = loadMutationList(mutFileName,numInAS,false);

		if (mutArray == null) {
			// Create the mutation list with estimated energies
			mutArray = new OneMutation[200000];
			RotamerSearch rs = new RotamerSearch(m,sysStrNum,hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst,doDihedE,doSolvationE,solvScale,softvdwMultiplier,rl);

			int curConfNum = 0;
			for(int i=0; i<numComb; i++) {
				// Reset each amino acid type
				System.out.print("Starting mutation combination " + i + " ... ");
				rs.refreshSystemStrand(); // clears allowables and does some other stuff

				for(int j=0; j<numInAS; j++) {
					if (residueMutatable[i][j] == 1) {
						String tempResAllow = (String)sParams.get("RESALLOWED"+j);
						for(int q=0;q<numTokens(tempResAllow);q++) {
							rs.setAllowable(residueMap[j],getToken(tempResAllow,q+1));
						}
					}
					else {
						rs.setAllowable(residueMap[j],resDefault[j]);
					}
				}
				
				// Perform simple mutation search for this set of mutatable residues
				curConfNum = rs.simpleMasterMutationSearch(residueMap,numInAS,
					resAllowed,numResAllowed,curConfNum,mutArray,targetVol-volWindow,
					targetVol+volWindow);
				System.out.println("finished");
			}
			
			System.out.println("Conformations remaining after volume filter "+curConfNum);
			
			// We now have all the mutations in mutArray, collapse the mutArray
			//  to the actual number of mutations we have.
			OneMutation newArray[] = new OneMutation[curConfNum];
			System.out.println("Allocated newArray");
			System.out.println("Initial Length of mutArray: "+mutArray.length);
			System.arraycopy(mutArray,0,newArray,0,curConfNum);
			mutArray = newArray;
			System.out.println("Trimmed Length of mutArray: "+mutArray.length);
			
			System.out.print("Removing duplicates...");
			mutArray = removeDuplicates(mutArray);
			System.out.println("done");
			
			System.out.println(mutArray.length+" unique mutations found in volume range "+(targetVol-volWindow)+" to "+(targetVol+volWindow));
			// Save mutation list
			saveMutationList(mutArray,mutFileName,false);
		}

	
		// Sort the mutation list
		// System.out.print("Sorting mutation list ... ");
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		// System.out.println("done");
		
		return mutArray;
	}


	// Reads the results of a partially completed run into an array
	//  of CommucObj. The MutationManager then queries this array
	//  before sending out a task.
	public OneMutation[] readResumeFile(OneMutation resumeResults[],
							   String resumeFilename, int numInAS, boolean distrDACS, boolean PEMcomp) {
	
		BufferedReader bufread = null;
		try {
			File file = new File(resumeFilename);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Resume File Not Found");
			return(null);
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).equalsIgnoreCase("Completed")) {
				if (PEMcomp) {//PEM computation
					resumeResults[resultNum].resMut = new int[numInAS];
					for(int q=0;q<numInAS;q++) 
						resumeResults[resultNum].resMut[q] = new Integer(getToken(str,8+q)).intValue();
					resumeResults[resultNum].flagMutType = getToken(str,8+numInAS);
				}
				else { //mutation search
					if (!distrDACS){ //Hybrid-K* or MinDEE/A* resume
						resumeResults[resultNum].score = new BigDecimal(getToken(str,5));
						resumeResults[resultNum].resTypes = new String[numInAS];
						for(int q=0;q<numInAS;q++) {
							resumeResults[resultNum].resTypes[q] = getToken(str,17+q);	
						}
					}
					else {//distributed DACS resume
						resumeResults[resultNum].mutNum = new Integer(getToken(str,5)).intValue();
						resumeResults[resultNum].score = new BigDecimal(getToken(str,9));
					}
				}
				resultNum++;
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		// Resize completed mutation array
		OneMutation temp[] = new OneMutation[resultNum];
		System.arraycopy(resumeResults,0,temp,0,resultNum);
		resumeResults = temp;
		return(resumeResults);
	}
	
	
	// Finds the index of the mutation in resumeResults with the same
	//  mutation sequence as the targetMutation. If none are found, -1
	//  is returned.
	public int findMutationIndex(OneMutation resumeResults[],
								 String targetMutation[]) {
		
		for(int q=0;q<resumeResults.length;q++) {
			if (resumeResults[q].isSame(targetMutation))
				return(q);
		}
		return(-1);
	}

	
	// Attempts to read a list of mutations from file
	public OneMutation[] loadMutationList(String fName, int numInAS, boolean PEMcomp) {
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... no mutation list file found. Computing one.");
			return(null);
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;
		OneMutation mutList[] = new OneMutation[1];
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else {
				if (PEMcomp) {//PEM computation
					mutList[resultNum] = new OneMutation();
					mutList[resultNum].resMut = new int[numInAS];
					
					for(int q=0;q<numInAS;q++) 
						mutList[resultNum].resMut[q] = new Integer(getToken(str,1+q)).intValue();
					
					mutList[resultNum].flagMutType = getToken(str,1+numInAS);
				}
				else {//mutation search
					mutList[resultNum] = new OneMutation();
					mutList[resultNum].score = new BigDecimal(getToken(str,1));
					mutList[resultNum].vol = new Float(getToken(str,2)).floatValue();
					mutList[resultNum].resTypes = new String[numInAS];
					for(int q=0;q<numInAS;q++) {
						mutList[resultNum].resTypes[q] = getToken(str,3+q);	
					}					
				}
				
				resultNum++;
				if (resultNum >= mutList.length){
					OneMutation newArray[] = new OneMutation[mutList.length+1000];
					System.arraycopy(mutList,0,newArray,0,resultNum);
					mutList = newArray;
				}
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}

		// Resize completed mutation array
		OneMutation temp[] = new OneMutation[resultNum];
		System.arraycopy(mutList,0,temp,0,resultNum);
		mutList = temp;
		System.out.println(" ... read "+mutList.length+" mutations from mutation list "+fName);
		return(mutList);
	}	

	// Saves the list of mutations so that a PEM computation/mutation search
	//  doesn't need to recompute these during a resume. Thus
	//  the resume can go more quickly.
	public void saveMutationList(OneMutation mutList[], String fName, boolean PEMcomp) {

		if (mutList.length == 0)
			return;
		
		int numInAS = 0;
		if (PEMcomp)
			numInAS = mutList[0].resMut.length;
		else
			numInAS = mutList[0].resTypes.length;

		PrintStream printStream = setupOutputFile(fName);
		for(int q=0;q<mutList.length;q++) {
			if (PEMcomp) {//PEM computation
				for(int w=0;w<numInAS;w++) {
					printStream.print(" "+mutList[q].resMut[w]);
				}
				printStream.print(" "+mutList[q].flagMutType);
				printStream.println();
			}
			else { //mutation search
				printStream.print(mutList[q].score + " " + mutList[q].vol);
				for(int w=0;w<numInAS;w++) {
					printStream.print(" "+mutList[q].resTypes[w]);
				}
				printStream.println();
			}
		}
		printStream.close();
	}

	//Removes duplicate mutations (for which the mutation sequence is the same) from a given list
	public OneMutation [] removeDuplicates(OneMutation mutArray[]){
		
		//First, sort the list alphabetically, according to the mutation sequence
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		
		//Copy mutArray into nArray, excluding duplicate entries
		OneMutation nArray[] = new OneMutation[mutArray.length];
		
		//Copy the first element
		nArray[0] = mutArray[0];
		int nAIndex = 1;
		
		//Compare each mutation with the previous one in the list
		for (int i=1; i<mutArray.length; i++){ //for each mutation
			if (!(mutArray[i].isSame(mutArray[i-1].resTypes))){ //different sequence
				nArray[nAIndex] = mutArray[i];
				nAIndex++;
			}
		}
		
		mutArray = new OneMutation[nAIndex];
		System.arraycopy(nArray,0,mutArray,0,mutArray.length);
		
		return mutArray;//return the reduced list
	}

	// Mutation search Slave function
	public CommucObj handleKSSlave(CommucObj cObj) {

		if (cObj.PEMcomp){ //PEM computation
			if (!cObj.entropyComp) //PEM computation
				cObj = handleComputeAllPairwiseRotamerEnergiesSlave(cObj);
			else //entropy E matrix computation
				handleDoResEntropySlave(cObj);
		}		
		else { //distributed mutation search
			if (cObj.distrDACS){ //running distributed DACS
				cObj = doDistrDACSSlave(cObj);
			}
			else if (cObj.distrDEE){ //running distributed DEE
				cObj = doDistrDEESlave(cObj);
			}
			else { //running Hybrid MinDEE-K*
				cObj = hybridKScompute(cObj);
			}
		}
		return cObj;
	}
	
	//Handles the mutation search for Hybrid MinDEE-K*
	private CommucObj hybridKScompute(CommucObj cObj){
		
		for(int runNum = 0; runNum<2; runNum++) {
			long startTime = System.currentTimeMillis();

			// First compute q_E, then q_EL
			boolean ligPresent = (runNum == 1);

			//Setup the molecule system
			Molecule m = new Molecule();
			int numLigRotamers = setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);

			RotamerSearch rs = null;
			if (ligPresent)
				rs = new RotamerSearch(m,sysStrNum,ligStrNum,cObj.isLigAA, hElect, hVDW, hSteric, true,
						true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst,cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl);		
			else
				rs = new RotamerSearch(m,sysStrNum, hElect, hVDW, hSteric, true,
						true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst,cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl);
			
			System.out.print("Loading precomputed min energy matrix...");
			rs.loadPairwiseEnergyMatrices(cObj.arpFilenameMin);
			System.out.println("done");
			
			if (cObj.doMinimization){ //MinDEE, so load the max matrix
				System.out.print("MinDEE: Loading precomputed max energy matrix...");
				rs.loadPairwiseEnergyMatricesMax(cObj.arpFilenameMax);
				System.out.println("done");
			}

			System.out.println("Beginning setAllowables");
			// Ligand allowable set in the RotamerSearch() constructor
			for(int q=0; q<cObj.numInAS; q++)
				rs.setAllowable(cObj.residueMap[q],cObj.currentMutation[q]);
			
			
			//Initially, no rotamers have been pruned
			boolean prunedRotAtRes[] = new boolean[rs.arpMatrix.length-1];
			for (int i=0; i<prunedRotAtRes.length; i++)
				prunedRotAtRes[i] = false;
			
			
			//Prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
			double stericE = Math.pow(10,10);
			prunedRotAtRes = rs.DoPruneStericTemplate(cObj.numInAS, cObj.numTotalRotamers, numLigRotamers, 
					cObj.residueMap, rotamerIndexOffset, prunedRotAtRes, stericE);			
		
			//Perform the DEE pruning
			if (cObj.doMinimization) //compute the MinDEE interval terms
				rs.doCompMinDEEIntervals(cObj.numInAS, cObj.numTotalRotamers, numLigRotamers, cObj.residueMap, 
						rotamerIndexOffset, prunedRotAtRes, cObj.scaleInt, cObj.maxIntScale);
			
			prunedRotAtRes = rs.DoDEEGoldstein(cObj.numInAS, cObj.numTotalRotamers, numLigRotamers, 
					cObj.residueMap,cObj.rotamerIndexOffset, cObj.initEw, prunedRotAtRes, cObj.doMinimization, false, cObj.minimizeBB);
			
			//Prune with MinBounds (last parameter is false)
			prunedRotAtRes = rs.DoMinBounds(cObj.numInAS,cObj.numTotalRotamers,numLigRotamers,
					cObj.residueMap,cObj.rotamerIndexOffset,stericE,prunedRotAtRes,cObj.initEw, false, false);
			
			//Compute the Ec value and prunedIsSteric[] (last parameter is true)
			rs.DoMinBounds(cObj.numInAS,cObj.numTotalRotamers,numLigRotamers,
					cObj.residueMap,cObj.rotamerIndexOffset,stericE,prunedRotAtRes,cObj.initEw, false, true);
		
			
			boolean usingInitialBest = ligPresent;
			BigDecimal initialBest = (new BigDecimal("0.0"));
			
			if (usingInitialBest)
				initialBest =  cObj.q_E.multiply(cObj.bestScore).multiply(new BigDecimal(cObj.gamma * cObj.epsilon));

			rs.slaveDoRotamerSearch(cObj.computeEVEnergy,cObj.doMinimization,cObj.numInAS,
				cObj.numResAllowed,cObj.numTotalRotamers,cObj.rotamerIndexOffset,resAllowed,
				cObj.residueMap,usingInitialBest,initialBest,cObj,cObj.minimizeBB,false,null);
			
			if ((cObj.repeatSearch)&&(rs.repeatSearch)){ //the desired accuracy was not achieved, so repeat the search: the setup is already done
				
				rs.repeatSearch = false; //reset the flag
				if (ligPresent)
					cObj.EL_repeatEw = true; //set the flag
				else
					cObj.E_repeatEw = true;
				
				rs.slaveDoRotamerSearch(cObj.computeEVEnergy,cObj.doMinimization,cObj.numInAS,
						cObj.numResAllowed,cObj.numTotalRotamers,cObj.rotamerIndexOffset,resAllowed,
						cObj.residueMap,usingInitialBest,initialBest,cObj,cObj.minimizeBB,false,null);
			}
				
			long stopTime = System.currentTimeMillis();
			if(runNum == 0)
				cObj.q_E_Time = Math.round((stopTime - startTime) / 1000.0f);
			else
				cObj.q_EL_Time = Math.round((stopTime - startTime) / 1000.0f);
		} // end for(runNum)
		
		return cObj;
	}
	
/////////////////////////////////////////////////////////////////////////
// MIN and MAX pairwise energy matrices computation
/////////////////////////////////////////////////////////////////////////
	
	//This function computes all min and max pairwise rotamer interaction energies
	// For every possible 2 mutation positions
	//  For every possible 19^2 amino-acid mutations for this position (other res removed)
	//   For every rotamer conformation
	//    Compute the pairwaise interaction energy, save all energies
	// The matrix has four types of entries
	//  M[i][0] = intra-residue energy for rotamer i
	//  M[0][i] = rotamer i to backbone (shell) energy
	//  M[i][j] = rotamer i to rotamer j energy
	//  M[0][0] = backbone (shell) to backbone (shell) energy
	public void handleComputeAllPairwiseRotamerEnergiesMaster(String s, Hashtable curScope, Hashtable sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new Hashtable(hashSize);
	
			// Read System parameters
			readConfigFile(sParams,getToken(s,2));
			
			// Read Mutation search parameters
			readConfigFile(sParams,getToken(s,3));
		}
		
		int numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		int ligNum = (new Integer((String)sParams.get("LIGNUM"))).intValue();
		int numMutations = 2; //pairwise energies are computed
		boolean doMinimize = (new Boolean((String)sParams.get("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.get("MINIMIZEBB"))).booleanValue();
		String runName = (String)sParams.get("RUNNAME");
		String mutFileName = (String)sParams.get("MUTFILENAME");
		String minEMatrixName = (String)sParams.get("MINENERGYMATRIXNAME");
		String maxEMatrixName = (String)sParams.get("MAXENERGYMATRIXNAME");
		
		boolean ligPresent = (new Boolean((String)sParams.get("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)sParams.get("LIGTYPE");

		boolean resumeSearch = (new Boolean((String)sParams.get("RESUMESEARCH"))).booleanValue();
		String resumeFilename = (String)sParams.get("RESUMEFILENAME");
		
		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
		System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
		System.out.println("Ligand Type: "+ligType);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);

		
		System.out.println("Computing _All_ Rotamer-Rotamer Energies");
		
		System.out.println("Starting minimum and maximum bound energy computation");
		
		if(resumeSearch) {
			System.out.println("** Resuming Search **");
			System.out.println("     resuming from file: "+resumeFilename);
		}
		
		// Create the mutation list with estimated energies
		OneMutation mutArray[] = new OneMutation[200000];
		
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]);
		}
		System.out.println();
		
		int mutEnerMatrixSize = totalNumRotamers*numInAS+numLigRotamers+1; // +1 for the backbone
		float mutationEnergiesMin[][] = new float[mutEnerMatrixSize][mutEnerMatrixSize];
		float mutationEnergiesMax[][] = new float[mutEnerMatrixSize][mutEnerMatrixSize];
		for(int i=0; i<mutEnerMatrixSize; i++) {
			for(int j=0; j<mutEnerMatrixSize; j++){
				mutationEnergiesMin[i][j] = 0.0f;
				mutationEnergiesMax[i][j] = 0.0f;
			}
		}
		
		//Look for previous mutation file
		System.out.println();
		System.out.print("Looking for mutation list file ");
		mutArray = loadMutationList(mutFileName,numInAS,true);
		if (mutArray == null) //generate mutation list		
			mutArray = getMutArrayPairEcomp(numInAS,ligPresent,mutFileName,minimizeBB);

		//Sort the mutation list
		System.out.print("Sorting mutation list ... ");
		RyanQuickSort rqs = new RyanQuickSort();
		rqs.Sort(mutArray);
		rqs = null;
		System.out.println("done");
		
		// If doing a resume, read the initial results into a bunch of OneMutations
		if (resumeSearch) {
			OneMutation resumeResults[] = new OneMutation[mutArray.length];
			for(int q=0;q<mutArray.length;q++)
				resumeResults[q] = new OneMutation();
			resumeResults = readResumeFile(resumeResults,resumeFilename,numInAS,false,true);
			System.out.println("Read "+resumeResults.length+" completed mutations");
			
			// Now filter removed mutations (already computed results
			//  are NOT written to file since we already have them)
			int newIndex = 0;
			OneMutation newArray2[] = new OneMutation[mutArray.length];
			for(int q=0;q<mutArray.length;q++) {
				int w = sampFindMutationIndex(resumeResults, mutArray[q].flagMutType, mutArray[q].resMut);
				if (w<0) {//this mutation has not been computed yet
					newArray2[newIndex] = mutArray[q];
					newIndex++;
				}
			}
			mutArray = new OneMutation[newIndex];
			System.arraycopy(newArray2,0,mutArray,0,newIndex);
			System.out.println("Length of mutArray after removing already computed mutations: "+mutArray.length);
			
			System.out.print("Loading partially precomputed energy matrices...");
			
			mutationEnergiesMin = (float[][])readObject(minEMatrixName+".dat");
			if (doMinimize) //energy minimization is performed
				mutationEnergiesMax = (float [][])readObject(maxEMatrixName+".dat");
			
			System.out.println("done");
		}
		
		MutationManager mutMan = new MutationManager(runName,mutArray,true);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigNum(ligNum);
		mutMan.setLigType(ligType);
		mutMan.setarpFilenameMin(minEMatrixName);
		mutMan.setPairEMatrixMin(mutationEnergiesMin);
		if (doMinimize){
			mutMan.setarpFilenameMax(maxEMatrixName);
			mutMan.setPairEMatrixMax(mutationEnergiesMax);
		}
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.numResAllowed(numAAallowed);
		mutMan.setIsLigAA(true);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setCalculateVolumes(false);
		mutMan.setLigPresent(ligPresent);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
	}
	
	// Mutation search Slave function
	public CommucObj handleComputeAllPairwiseRotamerEnergiesSlave(CommucObj cObj) {
				
		long startTime = System.currentTimeMillis();
						
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);
		
		RotamerSearch rs = null;
		if (ligPresent)
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,true, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl);
		else
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl);

		if (((cObj.flagMutType.compareTo("AS-AS")==0)||(cObj.flagMutType.compareTo("SHL-AS")==0)||(cObj.flagMutType.compareTo("TEMPL")==0))&&(ligPresent)){
			m.deleteStrand(ligStrNum);//we do not need the ligand for these runs
			rs = null;
			rs = new RotamerSearch(m,0, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult,rl);					
		}
		
		System.out.println("Beginning setAllowables");
		// Ligand allowable set in the RotamerSearch() constructor		
		for(int j=0; j<cObj.numInAS; j++) {
			if (cObj.resMut[j] == 1) {
				for(int q=0;q<resAllowed.length;q++) {
					rs.setAllowable(cObj.residueMap[j],resAllowed[q]);
				}
			}
			else {
				rs.setAllowable(cObj.residueMap[j],cObj.resDefault[j]);
			}
		}	
		
		// The goal is that the total energy of a system can be bounded by the sum of 
		//  all pairwise active site residue entries plus the entry for each active site
		//  residue's shell run plus each active site residue's self intra-residue energy.
		//  If a ligand is present then one should add the ligand to shell energy, the
		//  ligand to each active site residue pairwise energy, and the ligand self intra-
		//  residue energy.
		
		//Create a matrix that will contain both the min and max pairwise energy matrices;
		//	the min matrix is in the first half of the columns, the max is in the second half:
		//		minMaxEMatrix[i][j] = minMatrix[i][j]
		//		minMaxEMatrix[i][mutEMatrixSize + j] = maxMatrix[i][j]
		int mutEnerMatrixSize = cObj.numTotalRotamers*cObj.numInAS+cObj.numLigRotamers+1; // +1 for the backbone				
		float minMaxEMatrix[][] = new float[mutEnerMatrixSize][mutEnerMatrixSize * 2];
		for(int i=0; i<mutEnerMatrixSize; i++) {
			for(int j=0; j<mutEnerMatrixSize*2; j++){
				minMaxEMatrix[i][j] = 0.0f;
			}
		}
		
		boolean shellRun = false;
		boolean intraRun = false;
		boolean templateOnly = false;
		
		
		if (cObj.flagMutType.compareTo("TEMPL")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = true;
		}
		else if (cObj.flagMutType.compareTo("AS-AS")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = false;ligPresent = false;intraRun = false;templateOnly = false;
		}	
		else if (cObj.flagMutType.compareTo("SHL-AS")==0){
			
			// Then shell runs for the active site residues
			// Computes the active site residue rotamers to shell energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = false;				
		}
		else if (cObj.flagMutType.compareTo("INTRA")==0){
			
			// Compute all intra-residue energies					
			shellRun = false;intraRun = true;templateOnly = false;
		}				
		else if (cObj.flagMutType.compareTo("LIG-AS")==0){
			
			// **** Ligand present runs ****
			// This section computes the inter-residue energies between
			//  active site residues and the ligand
			shellRun = false;intraRun = false;templateOnly = false;			
		}
		else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)
			
			// Computes ligand rotamer to shell energies
			shellRun = true; intraRun = false;templateOnly = false;
		}
		
		//Compute the corresponding matrix entries
		minMaxEMatrix = rs.simplePairwiseMutationAllRotamerSearch(cObj.residueMap,
				cObj.numInAS,resAllowed,resAllowed.length,cObj.rotamerIndexOffset,
				cObj.numTotalRotamers,cObj.doMinimization,ligPresent,shellRun,intraRun,
				cObj.resMut,minMaxEMatrix,mutEnerMatrixSize,cObj.minimizeBB,templateOnly);
		
		
		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
		
		
		//Store the information in less space to allow the master node to buffer several cObj at once
		int numCompEntries = 0;				
		for(int i1=0; i1<mutEnerMatrixSize; i1++) {
			for(int i2=0; i2<mutEnerMatrixSize; i2++){
				if ((minMaxEMatrix[i1][i2]!=0.0f)||(minMaxEMatrix[i1][mutEnerMatrixSize + i2]!=0.0f)){
					//if the min is computed, then so should be the max
					numCompEntries++;
				}
			}
		}
		cObj.compEE = new SamplingEEntries[numCompEntries];//computed new entries
		for (int i1=0; i1<numCompEntries; i1++){
			cObj.compEE[i1] = new SamplingEEntries();
		}
		
		int curE = 0;
		for(int i1=0; i1<mutEnerMatrixSize; i1++) {
			for(int i2=0; i2<mutEnerMatrixSize; i2++){
				if ((minMaxEMatrix[i1][i2]!=0.0f)||(minMaxEMatrix[i1][mutEnerMatrixSize + i2]!=0.0f)){
					cObj.compEE[curE].index1 = i1;
					cObj.compEE[curE].index2 = i2;
					cObj.compEE[curE].minE = minMaxEMatrix[i1][i2];
					cObj.compEE[curE].maxE = minMaxEMatrix[i1][mutEnerMatrixSize + i2];
					curE++;
				}
			}
		}
		
		minMaxEMatrix = null;
		
		return cObj;
	}	
	
	//Generates and saves to file the mutation list for the pairwise energy matrix computation
	private OneMutation[] getMutArrayPairEcomp(int numInAS, boolean ligPresent, String mutFileName, boolean minimizeBB){
		
		final int numMutations = 2; //pairwise energy computation
		
		// Generate all combinations
		int numComb = factorial(numInAS).divide(factorial(numInAS-numMutations).multiply(factorial(numMutations))).intValue();
		int residueMutatable[][] = new int[numComb][numInAS];
		generateCombinations(residueMutatable,numInAS,numMutations);
		// At this point each row of residueMutatble is a 0/1 array which specifies a mutation 
		//  pair combination, 1 indicates that that residue can mutate in the specified combination
		
		System.out.println("Number of possible mutation combinations: "+numComb);
		
		// Create the mutation list with estimated energies
		OneMutation mutArray[] = new OneMutation[numComb];			
		
		//Set the AS-AS mutations
		int curMutNum = 0;
		for(int i=0; i<numComb; i++) {
			
			mutArray[i] = new OneMutation();
			mutArray[i].flagMutType = "AS-AS";
			mutArray[i].resMut = new int[numInAS];
			for(int j=0; j<numInAS; j++) {
				mutArray[i].resMut[j] = residueMutatable[i][j];
			}
			
			// Perform simple mutation search for this set of mutatable residues
			curMutNum++;
		}
		
		//Add the runs for template only, AS-shell, AS-ligand, intra-residue energies, and ligand-shell
		int t=0;
		if (minimizeBB)
			t = 1;
		int numOtherMut;
		if (ligPresent)
			numOtherMut = 2+t+2*numInAS;
		else
			numOtherMut = 1+t+numInAS;
		OneMutation otherMutArray[] = new OneMutation[numOtherMut];
		
		for (int i=0; i<numOtherMut; i++){
			otherMutArray[i] = new OneMutation();
			otherMutArray[i].resMut = new int[numInAS];
		}
		
		//Set the AS-shell mutations
		for (int i=0; i<numInAS; i++){
			otherMutArray[i].flagMutType = "SHL-AS";
			for (int j=0; j<numInAS; j++){
				if (i==j)
					otherMutArray[i].resMut[j] = 1;
				else
					otherMutArray[i].resMut[j] = 0;
			}
		}
		
		//Set the intra-residue energies run
		otherMutArray[numInAS].flagMutType = "INTRA";
		for (int j=0; j<numInAS; j++)
			otherMutArray[numInAS].resMut[j] = 1;
		
		//Set the template energy run
		if (minimizeBB){
			otherMutArray[numInAS+1].flagMutType = "TEMPL";
			for (int j=0; j<numInAS; j++)
				otherMutArray[numInAS+1].resMut[j] = 0;
		}
		
		if (ligPresent){//if the ligand is present, set the corresponding runs
			
			//Set the AS-ligand mutations
			for (int i=1+t+numInAS; i<=2*numInAS+t; i++){
				otherMutArray[i].flagMutType = "LIG-AS";
				for (int j=0; j<numInAS; j++){
					if ((i-1-t-numInAS)==j)
						otherMutArray[i].resMut[j] = 1;
					else
						otherMutArray[i].resMut[j] = 0;
				}
			}
			
			//Set the ligand-shell run
			otherMutArray[1+t+2*numInAS].flagMutType = "LIG-SHL";
			for (int j=0; j<numInAS; j++)
				otherMutArray[1+t+2*numInAS].resMut[j] = 0;
		}
		
		// We now have all the mutations in mutArray, collapse the mutArray
		//  to the actual number of mutations we have.
		OneMutation newArray[] = new OneMutation[curMutNum+numOtherMut];
		System.arraycopy(otherMutArray,0,newArray,0,numOtherMut);//add the other mutations first
		System.arraycopy(mutArray,0,newArray,numOtherMut,curMutNum);//then add the AS-AS mutations
		
		mutArray = newArray;
		System.out.println("Length of mutArray: "+mutArray.length);
		// Save mutation list
		saveMutationList(mutArray,mutFileName,true);
		
		return mutArray;
	}
	
	// Finds the index of the mutation in resumeResults with the same
	//  mutation sequence as resMut. If none are found, -1 is returned.
	public int sampFindMutationIndex(OneMutation resumeResults[], String flMutType, int mutResidues[]) {
		
		for(int q=0;q<resumeResults.length;q++) 
			if ((resumeResults[q].flagMutType.compareTo(flMutType)==0) && (sameSeq(resumeResults[q].resMut,mutResidues)))
				return(q);
			
		return(-1);
	}

	//Determines if the residues that can mutate are the same for two mutation sequences
	private boolean sameSeq (int computedResMut[], int allResMut[]){
		
		boolean found = true;
		for (int i=0; i<computedResMut.length; i++){
			if (computedResMut[i]!=allResMut[i])
				found = false;
		}
		return found;
	}
///////////////////////////////////////////////////////////////////////////
//	End of MIN and MAX Pairwise Energy Precomputation
///////////////////////////////////////////////////////////////////////////
	
////////////////////////////////////////////////////////////////
//	 Compute minimized-GMEC section
////////////////////////////////////////////////////////////////	
	// Computes the energy-minimized structure for the rotameric conformations
	//		specified by the input file
	public void handleMinDEEApplyRot(String s, Hashtable curScope) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Input results filename (string)
		// 3: Run name (for output files)
		// 4: Ligand (boolean), is true if present
		// 5: Amino acid type for ligand (if ligand is absent, write none or anything)
		// 6: Backbone minimization (boolean)
		// 7: Output minimized structures (PDB files) (boolean)

		int numInAS = -1;

		// The system paramaters
		Hashtable sParams = new Hashtable(hashSize);

		// Read System Parameters
		readConfigFile(sParams,getToken(s,2));

		// Pull search parameters
		numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		int ligNum = (new Integer((String)sParams.get("LIGNUM"))).intValue();
		String confResFile = getToken(s,3);
		String runName = getToken(s,4);
		boolean ligPresent = parseToBoolean(getToken(s,5),curScope);
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,6);
		int numResults = parseToInt(getToken(s,7),curScope);
		boolean minimizeBB = parseToBoolean(getToken(s,8),curScope);
		boolean outputPDB = parseToBoolean(getToken(s,9),curScope);
		
		int numResidues;
		if (ligPresent)
			numResidues = numInAS+1;
		else
			numResidues = numInAS;
		
		//Read the results file into the AA and rot matrices
		String AAtypes[] = new String[numResults*numResidues];
		int rotNums[] = new int[numResults*numResidues];
		readRotResFile(confResFile,AAtypes,rotNums,numResults,numResidues);
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			System.out.print(" "+residueMap[i]);
		}
		System.out.println();
		
		String curSeq[] = new String[numInAS];
		
		PrintStream logPS = setupOutputFile(runName);
		
		int numSaved = 0;
		for (int curResult=0; curResult<numResults; curResult++){
			System.out.print("Starting minimization of result "+(curResult+1)+"..");
			
			Amber96ext a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, softvdwMultiplier);
			
			StrandRotamers sysLR = null;
			StrandRotamers ligLR2 = null;
			
			sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
			if(ligPresent)
				ligLR2 = new StrandRotamers(rl,m.strand[ligStrNum]);
			
			//the starting index for the current result within AAtypes[] and rotNums[]
			int startInd = numResidues*curResult;
		
			//Set the allowables for the AS residues
			for(int j=0;j<numInAS;j++){
				curSeq[j] = AAtypes[startInd+j];
				sysLR.setAllowable(residueMap[j],curSeq[j]);
				sysLR.changeResidueType(m,residueMap[j],curSeq[j],true,true);
				m.residue[residueMap[j]].flexible = true;
			}
			if(ligPresent){
				ligLR2.setAllowable(0,ligType); //the only allowable AA type for the ligand
				ligLR2.changeResidueType(m,0,ligType,true,true);
				m.strand[ligStrNum].residue[0].flexible = true;
			}
			
			a96ff.calculateTypesWithTemplates();
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			a96ff.setLigandNum(ligNum);
			
			int curAA[] = new int[m.numberOfResidues];
			for(int j=0;j<m.strand[sysStrNum].numberOfResidues;j++){
				curAA[j] = sysLR.getIndexOfNthAllowable(j,0);
			}
			
			SimpleMinimizer simpMin = null;
			BBMinimizer bbMin = null;
			if (!minimizeBB){
				simpMin = new SimpleMinimizer();
				if(ligPresent)
					simpMin.initialize(m,0,1,a96ff,sysLR,ligLR2,curAA,ligLR2.getIndexOfNthAllowable(0,0),doDihedE,rl);
				else
					simpMin.initialize(m,0,a96ff,sysLR,curAA,doDihedE,rl);
			}
			else {
				bbMin = new BBMinimizer();
				if (ligPresent)
					bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
				else
					bbMin.initialize(m, a96ff, residueMap, sysStrNum);
			}
			
			m.backupAtomCoord();
			//Apply the corresponding rotamers
			for (int j=0; j<numInAS; j++){										
				if (rl.getNumRotForAAtype(curAA[residueMap[j]])!=0){//not GLY or ALA
					int curRot = rotNums[startInd+j];
					sysLR.applyRotamer(m, residueMap[j], curRot);
				}
			}				
			if (ligPresent){ //apply the ligand rotamer
				if (rl.getNumRotForAAtype(ligLR2.getIndexOfNthAllowable(0,0))!=0){//not GLY or ALA
					int curRot = rotNums[startInd+numInAS]; //the ligand rotamer
					ligLR2.applyRotamer(m, 0, curRot);//the ligand level
				}
			}	
			
			double unMinE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy before minimization
			
			if (!minimizeBB)
				simpMin.minimize(35,false);
			else
				bbMin.minimizeFull(false);
			
			double minE[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //the energy after minimization
			if ((!minimizeBB)&&(doDihedE)) //add dihedral energies
				minE[0] += simpMin.computeDihedEnergy();
			
			if (outputPDB){ //save molecule
				saveMolecule(m,"pdbs/savedMol"+(numSaved+1)+".pdb",(float)minE[0]);
				numSaved++;
			}
			
			m.restoreAtomCoord();
			m.updateCoordinates();	
			
			if (minE[0]>unMinE[0])
				minE = unMinE;
			
			logPS.print((curResult+1)+" ");
			for (int j=0; j<numInAS; j++){
				logPS.print(AAtypes[startInd+j]+" ");
			}
			logPS.println(unMinE[0]+" "+minE[0]+" ( "+minE[1]+" "+minE[2]+" "+minE[3]+" )");logPS.flush();
			
			System.out.println("done");
		}
		logPS.close();
		System.out.println("done");
	}
	
	private void readRotResFile (String resFile, String AAtypes[], int rotNums[], int numResults, int numResidues){
		
		BufferedReader bufread = null;
		try {
			File file = new File(resFile);		
			FileReader fr = new FileReader(file);  
			bufread = new BufferedReader(fr); 
		}
		catch (FileNotFoundException e) {
			System.out.println("ERROR: Results File Not Found");
		}

		boolean done = false;
		String str = null;
		int resultNum = 0;

		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else {
				for (int i=0; i<numResidues; i++){
					AAtypes[resultNum*numResidues+i] = getToken(str,2+i);
					rotNums[resultNum*numResidues+i] = new Integer((String)getToken(str,2+numResidues+i)).intValue();
				}
				resultNum++;
			}
		}
		
		if (numResults!=resultNum){
			System.out.println("Error: Not all results available for reading");
			System.exit(0);
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
	}
////////////////////////////////////////////////////////////////
//	 End of Compute minimized-GMEC section
////////////////////////////////////////////////////////////////

	
	
	
///////////////////////////////////////////////////////////////////////////
//	DEE section
///////////////////////////////////////////////////////////////////////////
	//Executes DEE (Traditional or MinDEE)
	public void handleDoDEE(String s, Hashtable curScope) {
		
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: DEE config filename (string)
		
		long startTime = System.currentTimeMillis();
		
		// The system paramaters
		Hashtable sParams = new Hashtable(hashSize);

		System.out.println("Performing DEE");

		// Read System Parameters
		readConfigFile(sParams,getToken(s,2));
		readConfigFile(sParams,getToken(s,3));
		
		// Pull search parameters
		int numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		int ligNum = (new Integer((String)sParams.get("LIGNUM"))).intValue();
		String runName = ((String)sParams.get("RUNNAME"));
		int numMaxMut = (new Integer((String)sParams.get("NUMMAXMUT"))).intValue();
		int algOption = (new Integer((String)sParams.get("ALGOPTION"))).intValue();
		int numSplits = (new Integer((String)sParams.get("NUMSPLITS"))).intValue();
		boolean useFlags = (new Boolean((String)sParams.get("SPLITFLAGS"))).booleanValue();
		boolean distrDACS = (new Boolean((String)sParams.get("DISTRDACS"))).booleanValue();
		boolean distrDEE = (new Boolean((String)sParams.get("DISTRDEE"))).booleanValue();
		boolean doMinimize = (new Boolean((String)sParams.get("DOMINIMIZE"))).booleanValue();
		boolean minimizeBB = (new Boolean((String)sParams.get("MINIMIZEBB"))).booleanValue();
		boolean approxMinGMEC = (new Boolean((String)sParams.get("APPROXMINGMEC"))).booleanValue();
		boolean preprocPairs = (new Boolean((String)sParams.get("PREPROCPAIRS"))).booleanValue();
		boolean scaleInt = (new Boolean((String)sParams.get("SCALEINT"))).booleanValue();
		String runNameEMatrixMin = (String)(sParams.get("MINENERGYMATRIXNAME"));
		String runNameEMatrixMax = (String)(sParams.get("MAXENERGYMATRIXNAME"));
		float initEw = (new Float((String)sParams.get("INITEW"))).floatValue();
		float lambda = (new Float((String)sParams.get("LAMBDA"))).floatValue();
		float pruningE = (new Float((String)sParams.get("PRUNINGE"))).floatValue();
		double stericE = (new Double((String)sParams.get("STERICE"))).doubleValue();
		float pairSt = (new Float((String)sParams.get("PAIRST"))).floatValue();
		float maxIntScale = (new Float((String)sParams.get("MAXINTSCALE"))).floatValue();
		double minRatioDiff = (new Double((String)sParams.get("MINRATIODIFF"))).doubleValue();
		int diffFact = (new Integer((String)sParams.get("DIFFFACT"))).intValue();
		int maxDepth = (new Integer((String)sParams.get("MAXDEPTH"))).intValue();
		String outputConfInfo = (String)(sParams.get("OUTPUTCONFINFO"));
		String outputPruneInfo = (String)(sParams.get("OUTPUTPRUNEINFO"));
		boolean ligPresent = (new Boolean((String)sParams.get("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.get("LIGTYPE"));
		boolean useEref = (new Boolean((String)sParams.get("USEEREF"))).booleanValue();
		boolean resumeSearch = (new Boolean((String)sParams.get("RESUMESEARCH"))).booleanValue();
		String resumeFilename = ((String)sParams.get("RESUMEFILENAME"));
		
		
		if ((!mpiRun)&&((distrDACS)||distrDEE)){
			System.out.println("ERROR: Distributed computation requires MPI");
			System.exit(1);
		}
		
		if (minimizeBB) //backbone minimization
			doMinimize = true;
			
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+resDefault[i]+")");
		}
		System.out.println();
		
		RotamerSearch rs = null;
		if (ligPresent)
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,true, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);
		else
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);

		System.out.print("Loading precomputed energy matrix...");

		rs.loadPairwiseEnergyMatrices(new String(runNameEMatrixMin+".dat"));
		if (doMinimize)
			rs.loadPairwiseEnergyMatricesMax(new String(runNameEMatrixMax+".dat"));
		System.out.println("done");
		
		float eRef[] = null;
		if (useEref) { //add the reference energies to the min (and max) intra-energies
			eRef = getResEntropyEmatricesEref(sParams, useEref);
			rs.addEref(eRef, doMinimize, ligPresent, numInAS);
		}
		

		/////////////////////////////////////////////////////////////
		// DEE section
		
		//Set the allowable AAs for each AS residue
		for (int j=0; j<numInAS; j++){
			String tempResAllow = (String)sParams.get("RESALLOWED"+j);
			for(int q=0;q<numTokens(tempResAllow);q++) {
				rs.setAllowable(residueMap[j],getToken(tempResAllow,q+1));
			}
			
			rs.setAllowable(residueMap[j],resDefault[j]); //the default type is set last
		}
		
		boolean prunedRotAtRes[] = new boolean [numInAS*totalNumRotamers + numLigRotamers];
		
		int msp[] = null; //major split positions (for DACS)
		
		final String rotFile = ("rot_out"+System.currentTimeMillis()); //output the pruned rotamers (for distributed DEE)
		final String sfFile = ("sf_matrix"+System.currentTimeMillis()); //output the split flags (for distributed DEE)
		
		
		//first prune all rotamers that are incompatible with the template (intra E + re-to-template E >= stericE)
		System.out.println();
		System.out.println("Pruning all rotamers incompatible with the template..");			
		prunedRotAtRes = rs.DoPruneStericTemplate(numInAS, totalNumRotamers, numLigRotamers, residueMap,
			rotamerIndexOffset, prunedRotAtRes, stericE);
		System.out.println();
		
		//preprocess pairs of rotamers (mark pairwise energies greater than the cutoff as steric clashes)
		if (preprocPairs){
			System.out.println("Preprocessing pairs of rotamers, cutoff of "+pairSt);
			rs.preprocessPairs(pairSt);
			System.out.println();
		}
		
		
		//Setup and do the DEE pruning
		if ((useFlags)||(algOption>=3))
			rs.setSplitFlags(prunedRotAtRes.length);//initialize the split flags
		
		int numPrunedRot = 0;
		int numPrunedPairs = 0;
		int numPrunedRotThisRun = 0;
		int numPrunedPairsThisRun = 0;
		boolean done = false;
		int numRuns = 1;
		
		while (!done){ //repeat the pruning cycle until no more rotamers are pruned	
			
			numPrunedRotThisRun = 0; 
			numPrunedPairsThisRun = 0;
			
			System.out.println("Starting DEE cycle run: "+numRuns);
			
			if (doMinimize) //precompute the interval terms in the MinDEE criterion
				rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, rotamerIndexOffset, 
						prunedRotAtRes, scaleInt, maxIntScale);
		
			//Depending on the chosen algorithm option, apply the corresponding pruning criteria;			
			System.out.println("Starting pruning with DEE (simple Goldstein)");		
			prunedRotAtRes = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
					rotamerIndexOffset, initEw, prunedRotAtRes, doMinimize, useFlags, minimizeBB);			
			System.out.println();
			
			if ((algOption>=3)){ //simple Goldstein pairs
				System.out.println("Starting pruning with DEE (mb pairs)");			
				rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
						rotamerIndexOffset, initEw, prunedRotAtRes, null, doMinimize, useFlags, true, false, minimizeBB, scaleInt, maxIntScale);
				System.out.println();
			}
			
			if ((useFlags)||(algOption>=3)){
				System.out.println("Starting pruning with Bounding Flags");			
				rs.DoBoundFlags(numInAS, totalNumRotamers, numLigRotamers, residueMap,
					rotamerIndexOffset, pruningE, prunedRotAtRes, initEw, useFlags);
				System.out.println();
			}
			
			System.out.println("Starting pruning with DEE (1-sp split-DEE)");			
			prunedRotAtRes = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
					rotamerIndexOffset, initEw, prunedRotAtRes, null, msp, -1, doMinimize, useFlags, 1, false, minimizeBB);
			System.out.println();
			
			System.out.println("Starting pruning with Bounds");			
			prunedRotAtRes = rs.DoMinBounds(numInAS, totalNumRotamers, numLigRotamers, residueMap,
				rotamerIndexOffset, pruningE, prunedRotAtRes, initEw, useFlags, false);
			System.out.println();
			
			//check how many rotamers/pairs are pruned this run
			int numTotalPrunedRot = countPrunedRot(prunedRotAtRes);
			int numTotalPrunedPairs = 0;
			if ((useFlags)||(algOption>=3)) //pairs pruning is performed
				numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
			
			if ((numTotalPrunedRot!=numPrunedRot)||(numTotalPrunedPairs!=numPrunedPairs)) { //new rotamers/pairs pruned this run
				
				numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
				numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
				
				numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
				numPrunedPairs = numTotalPrunedPairs;
				numRuns++;
			}
			else { //no more rotamers pruned, so perform the computationally-expensive 2-sp split-DEE and pairs
				
				if (numSplits<=1){ //DACS will not be performed
					
					if ((algOption>=3)){ //simple Goldstein pairs
						System.out.println("Starting pruning with DEE (full pairs)");	
						
						if (distrDEE){ //distributed full Goldstein pairs DEE
							doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
									ligNum, ligType, sParams, ligPresent, resDefault,
									runNameEMatrixMin, runNameEMatrixMax, initEw, prunedRotAtRes, doMinimize, 
									rs, sfFile, rotFile, useFlags, -1, msp, optPairs, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
							
							rs.setSplitFlags((boolean [][])readObject(sfFile)); //get the DE pairs from the distributed run
						}
						else { //perform on a single processor
							rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
									rotamerIndexOffset, initEw, prunedRotAtRes, null, doMinimize, useFlags, 
									false, false, minimizeBB, scaleInt, maxIntScale);
						}
						System.out.println();
					}
					
					if ((algOption>=2)){ //2-sp conf splitting
						System.out.println("Starting pruning with DEE (2-sp split-DEE)");
						
						if (distrDEE){ //distributed 2-sp split-DEE
							doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
									ligNum, ligType, sParams, ligPresent, resDefault,
									runNameEMatrixMin, runNameEMatrixMax, initEw, prunedRotAtRes, doMinimize, 
									rs, sfFile, rotFile, useFlags, 2, msp, optSplit, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
							
							prunedRotAtRes = ((boolean [])readObject(rotFile));
						}
						else { //perform on a single processor			
							prunedRotAtRes = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
									rotamerIndexOffset, initEw, prunedRotAtRes, null, msp, -1, doMinimize, useFlags, 2, false, minimizeBB);
						}
						System.out.println();
					}
					
					//check if 2-sp split-DEE and pairs pruned new rotamers
					numTotalPrunedRot = countPrunedRot(prunedRotAtRes);
					numTotalPrunedPairs = 0;
					if ((useFlags)||(algOption>=3)) //pairs pruning is performed
						numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
					
					if ((numTotalPrunedRot==numPrunedRot)&&(numTotalPrunedPairs==numPrunedPairs)) //no more rotamers/pairs pruned
						done = true;
					else { //new rotamers pruned this run
						numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
						numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
						
						numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
						numPrunedPairs = numTotalPrunedPairs;
						numRuns++;
					}
				}
				else //DACS will be performed
					done = true;
			}
			System.out.println("Num pruned rot this run: "+numPrunedRotThisRun);
			System.out.println("Num pruned pairs this run: "+numPrunedPairsThisRun);
			System.out.println();
		}
			
		if (numSplits<=1){ //DACS will not be performed, so do A* search here
			double bestScore = Math.pow(10,38); //for DEE/A*, the initial best score is the highest possible
			rs.doAStarGMEC(outputConfInfo,true,doMinimize,
					numInAS,totalNumRotamers,rotamerIndexOffset,residueMap,resDefault,numMaxMut,initEw,
					bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,eRef);
		}
		else { //DACS			
			int numRotForRes[] = compNumRotForRes(numInAS, rs, numLigRotamers, residueMap);
			BigInteger numInitUnprunedConfs = compNumUnprunedConfs(numInAS, totalNumRotamers, prunedRotAtRes, numLigRotamers,
					numRotForRes, ligPresent);			
			
			msp = new int[maxDepth];
			for (int i=0; i<maxDepth; i++)
				msp[i] = -1;
			
			if (distrDACS) { //distributed DACS (only for level 0)
				
				distrDEE = false; //do not perform both distributed DACS and distributed DEE
				
				//choose the major split position at level 0
				msp[0] = chooseSplitPos(numInAS,prunedRotAtRes,totalNumRotamers,numRotForRes,
						msp, 0, minRatioDiff);//the splitting position
				
				OneMutation resumeResults[] = null;			
				if (resumeSearch){ //read resume results
					System.out.println("Reading resume results..");
					resumeResults = new OneMutation[numRotForRes[msp[0]]];
					for(int q=0;q<resumeResults.length;q++)
						resumeResults[q] = new OneMutation();
					resumeResults = readResumeFile(resumeResults,resumeFilename,numInAS,true,false);
					System.out.println("Read "+resumeResults.length+" completed partitions.");
					System.out.println();
				}
				
				doDistrDACSMaster(runName, numInAS, rs, numLigRotamers, residueMap, ligPresent, resDefault,
						rotFile, prunedRotAtRes, algOption, sfFile, useFlags, initEw, pruningE, maxDepth, msp,
						numInitUnprunedConfs, diffFact, outputPruneInfo, outputConfInfo, minRatioDiff, doMinimize,
						runNameEMatrixMin, runNameEMatrixMax,ligNum, ligType, sParams, approxMinGMEC, 
						lambda, numRotForRes, resumeResults, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef);
			}
			else { //single-processor DACS
				
				PrintStream logPS = setupOutputFile(outputPruneInfo);
			
				doDACS(numInAS, rs, numLigRotamers, residueMap, ligPresent, resDefault,
						rotFile, prunedRotAtRes, algOption, sfFile, useFlags, initEw, pruningE, maxDepth, 0, logPS, msp,
						numInitUnprunedConfs, diffFact, outputConfInfo, minRatioDiff, doMinimize,runNameEMatrixMin,
						runNameEMatrixMax, distrDEE, ligNum, ligType, sParams, approxMinGMEC, lambda, false, -1, 
						null, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef);
			}
		}
		
		long stopTime = System.currentTimeMillis();
		
		System.out.println("Execution time: "+((stopTime-startTime)/(60.0*1000.0)));
		System.out.println("DEE done");
		//end of DEE section
		/////////////////////////////////////////////////////////////
	}
	
	//Do DACS; rs and sysLR in rs are assumed to be valid
	private void doDACS(int numInAS, RotamerSearch rs, int numLigRotamers, int residueMap[],
			boolean ligPresent, String resDefault[], String rotFile,
			boolean prunedRotAtRes[], int algOption, String sfFile, boolean useFlags, float initEw, float pruningE,
			int maxDepth, int curDepth, PrintStream logPS, int majorSplitPos[], BigInteger numInitUnprunedConfs,
			int diffFact, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM, 
			String maxPEM, boolean distrDEE, int ligNum, String ligType, Hashtable sParams, 
			boolean approxMinGMEC, float lambda, boolean distrDACS, int distrPartIndex, CommucObj cObj, 
			boolean minimizeBB, int numMaxMut, boolean scaleInt, float maxIntScale, boolean useEref, float eRef[]){
		
		if (curDepth>=maxDepth)
			return;
		
		System.out.println("Starting pruning with DEE (DACS).");
		
		//prunedRotAtRes[] should not be modified here, in order to be able to distinguish
		//	newly pruned rotamers and rotamers pruned by Bounds or DEE
		
		//the num rotamers for each AS residue and the ligand (if present);
		//	sysLR in rs must be valid (with all the possible AA's for each residue position)
		int numRotForRes[] = compNumRotForRes(numInAS, rs, numLigRotamers, residueMap);
		
		if (!distrDACS) {//not distributed DACS, so majorSplitPos[curDepth] is unknown; compute it
			majorSplitPos[curDepth] = chooseSplitPos(numInAS,prunedRotAtRes,totalNumRotamers,numRotForRes,
					majorSplitPos, curDepth, minRatioDiff);//the splitting position
		}
		
		int numPartitions = numRotForRes[majorSplitPos[curDepth]]; //the number of partitions
		int numPrunedPartitions = 0; //the number of partitions with no unpruned confs
		
		System.out.println("Current depth: "+curDepth);
		System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
		System.out.println();
		
		//map rotamer index to rot num for the splitting residue
		int indexMap[] = getIndexMap(numPartitions,rs,majorSplitPos[curDepth],residueMap,totalNumRotamers,rotamerIndexOffset);
		
		//count the total num confs and the num unpruned confs
		BigInteger numConfsTotalForPartition[] = new BigInteger[numPartitions];
		BigInteger numConfsUnprunedForPartition[] = new BigInteger[numPartitions];
		//BigInteger numInitConfsUnprunedForPartition[] = new BigInteger[numPartitions];
		BigInteger numTotalConfs = new BigInteger("0");
		BigInteger numUnprunedConfs = new BigInteger("0");
		BigInteger numEvaluatedConfs = new BigInteger("0");
		
		//update the best energy found
		pruningE = (float)Math.min(pruningE, rs.getBestE().doubleValue());
		double bestScore = pruningE;
		
		boolean savedSpFlags[][] = null;
		if (useFlags)
			savedSpFlags = rs.getSplitFlags(savedSpFlags);
		
		//determine the prunings for each of the sub-solutions (the partitions)
		boolean prunedForPartition[][] = new boolean[numPartitions][prunedRotAtRes.length];
		for (int i=0; i<prunedForPartition.length; i++){
			
			if ((!distrDACS)||(indexMap[i]==distrPartIndex)) { //not distributed DACS or current partition is the partition distributed for computation
			
				//copy the prunings from before conf splitting (from Bounds and simple traditional DEE)
				System.arraycopy(prunedRotAtRes,0,prunedForPartition[i],0,prunedRotAtRes.length);
				
				int partIndex = indexMap[i]; //the rotamer index of the partitioning rotamer
				
				//artificially set all rotamers at the splitting position, other than the rotamer 
				//	for the current partition, to pruned, so that there will be only one rotamer at
				//	that residue position when the conf splitting criterion is applied;
				//	when done, subtract the artifcially pruned rotamers from the total number of pruned
				boolean indToUnprune[] = new boolean[numPartitions];
				for (int j=0; j<indToUnprune.length; j++)
					indToUnprune[j] = false;
				
				//check the partition only if the current partitioning rotamer is not already pruned
				if (!prunedRotAtRes[partIndex]){
					
					for (int j=0; j<numPartitions; j++){
						if (j!=i) {//not the rotamer for the current partition
							int curInd = indexMap[j]; //the index of the current rotamer
							if (!prunedForPartition[i][curInd]){ //not pruned by the other DEE methods
								prunedForPartition[i][curInd] = true;
								indToUnprune[j] = true;
							}
						}
					}
					
					if (useFlags)
						rs.setSplitFlags(savedSpFlags); //for each partition, reset the flags to the globally valid ones
					
					int numPrunedRot = countPrunedRot(prunedForPartition[i]);
					int numPrunedPairs = countPrunedPairs(rs.getSplitFlags());
					int numPrunedRotThisRun = 0;
					int numPrunedPairsThisRun = 0;
					boolean done = false;
					int numRuns = 1;
					
					while (!done){ //repeat the pruning cycle until no more rotamers are pruned	
						
						numPrunedRotThisRun = 0; 
						numPrunedPairsThisRun = 0;
						
						System.out.println("Starting DEE cycle run: "+numRuns);
						
						if (doMinimize) //precompute the interval terms in the MinDEE criterion
							rs.doCompMinDEEIntervals(numInAS, totalNumRotamers, numLigRotamers, residueMap, rotamerIndexOffset, 
									prunedForPartition[i], scaleInt, maxIntScale);
						
						System.out.println("Starting pruning with DEE (simple Goldstein)");		
						prunedForPartition[i] = rs.DoDEEGoldstein(numInAS, totalNumRotamers, numLigRotamers, residueMap,
								rotamerIndexOffset, initEw, prunedForPartition[i], doMinimize, useFlags, minimizeBB);			
						System.out.println();
						
						if ((algOption>=3)){ //simple Goldstein pairs
							System.out.println("Starting pruning with DEE (mb pairs)");			
							rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
									rotamerIndexOffset, initEw, prunedForPartition[i], null, doMinimize, 
									useFlags, true, false, minimizeBB, scaleInt, maxIntScale);
							System.out.println();
						}
						
						if ((useFlags)||(algOption>=3)){
							System.out.println("Starting pruning with Bounding Flags");			
							rs.DoBoundFlags(numInAS, totalNumRotamers, numLigRotamers, residueMap,
								rotamerIndexOffset, pruningE, prunedForPartition[i], initEw, useFlags);
							System.out.println();
						}
							
						System.out.println("Starting pruning with DEE (1-sp split-DEE)");
						prunedForPartition[i] = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
								rotamerIndexOffset, initEw, prunedForPartition[i], null, majorSplitPos, curDepth, doMinimize, 
								useFlags, 1, false, minimizeBB);
						System.out.println();
						
						System.out.println("Starting pruning with Bounds");			
						prunedForPartition[i]= rs.DoMinBounds(numInAS, totalNumRotamers, numLigRotamers, residueMap,
							rotamerIndexOffset, pruningE, prunedForPartition[i], initEw, useFlags, false);
						System.out.println();
						
						//check how many rotamers/pairs are pruned this run
						int numTotalPrunedRot = countPrunedRot(prunedForPartition[i]);
						int numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
						
						if ((numTotalPrunedRot!=numPrunedRot)||(numTotalPrunedPairs!=numPrunedPairs)) { //new rotamers/pairs pruned this run
							
							numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
							numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
							
							numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
							numPrunedPairs = numTotalPrunedPairs;
							numRuns++;
						}
						else { //no more rotamers pruned, so perform the computationally-expensive 2-sp split-DEE and pairs
							
							if ((algOption>=3)&&(curDepth>=(maxDepth-1))){ //simple Goldstein pairs					
								System.out.println("Starting pruning with DEE (full pairs)");	
								
								if (distrDEE){ //distributed full Goldstein pairs DEE
									doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
											ligNum, ligType, sParams, ligPresent, resDefault, 
											minPEM, maxPEM, initEw, prunedForPartition[i], doMinimize, rs, sfFile, rotFile, 
											useFlags, -1, majorSplitPos, optPairs, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
									
									rs.setSplitFlags((boolean [][])readObject(sfFile)); //get the DE pairs from the distributed run
								}
								else { //perform on a single processor
									rs.DoDEEPairs(numInAS, totalNumRotamers, numLigRotamers, residueMap,
											rotamerIndexOffset, initEw, prunedForPartition[i], null, doMinimize, 
											useFlags, false, false, minimizeBB, scaleInt, maxIntScale);
								}
								System.out.println();
							}
							
							if ((algOption>=2)){ //2-sp conf splitting					
								System.out.println("Starting pruning with DEE (2-sp split-DEE)");
								
								if (distrDEE) { //distributed 2-sp split-DEE
									doDistrDEEMaster(numInAS, numLigRotamers, residueMap, 
											ligNum, ligType, sParams, ligPresent, resDefault,
											minPEM, maxPEM, initEw, prunedForPartition[i], doMinimize, rs, sfFile, rotFile, 
											useFlags, 2, majorSplitPos, optSplit, minimizeBB, scaleInt, maxIntScale, useEref, eRef);
									
									prunedForPartition[i] = ((boolean [])readObject(rotFile));
								}
								else {
									prunedForPartition[i] = rs.DoDEEConfSplitting(numInAS, totalNumRotamers, numLigRotamers, residueMap,
											rotamerIndexOffset, initEw, prunedForPartition[i], null, majorSplitPos, curDepth, doMinimize, 
											useFlags, 2, false, minimizeBB);
								}
								System.out.println();
							}
							
							//check if 2-sp split-DEE and pairs pruned new rotamers/pairs
							numTotalPrunedRot = countPrunedRot(prunedForPartition[i]);
							numTotalPrunedPairs = countPrunedPairs(rs.getSplitFlags());
							
							if ((numTotalPrunedRot==numPrunedRot)&&(numTotalPrunedPairs==numPrunedPairs)) //no more rotamers/pairs pruned
								done = true;
							else { //new rotamers pruned this run
								numPrunedRotThisRun += numTotalPrunedRot - numPrunedRot;
								numPrunedPairsThisRun += numTotalPrunedPairs - numPrunedPairs;
								
								numPrunedRot = numTotalPrunedRot; //update the number of pruned rotamers
								numPrunedPairs = numTotalPrunedPairs;
								numRuns++;
							}
						}
						System.out.println("Num pruned rot this run: "+numPrunedRotThisRun);
						System.out.println("Num pruned pairs this run: "+numPrunedPairsThisRun);
						System.out.println();
					}
				}
					
				//count the number of pruned rotamers for each residue (except for the splitting residue,
				//	which always has only 1 available rotamer for the current partition)
				int numPrunedRotForRes[] = new int[numInAS]; //after pruning for this partition
				//int numInitPrunedRotForRes[] = new int[numInAS]; //initial prunings, before pruning for this partition
				for (int j=0; j<numInAS; j++){
					numPrunedRotForRes[j] = 0;
					if (j!=majorSplitPos[curDepth]){
						for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
							if (prunedForPartition[i][j*totalNumRotamers + k])
								numPrunedRotForRes[j]++;
							//if (prunedRotAtRes[j*totalNumRotamers + k])
								//numInitPrunedRotForRes[j]++;
						}
					}
					else {// j==majorSplitPos
						numPrunedRotForRes[j] = 0;
						//numInitPrunedRotForRes[j] = 0;
					}
				}
				int numPrunedLigRot = 0;
				//int numInitPrunedLigRot = 0;
				for (int k=0; k<numLigRotamers; k++){
					if (prunedForPartition[i][numInAS*totalNumRotamers + k])
						numPrunedLigRot++;
					//if (prunedRotAtRes[numInAS*totalNumRotamers + k])
					//	numInitPrunedLigRot++;
				}
				
				//count the total num confs and the num unpruned confs
				numConfsTotalForPartition[i] = new BigInteger("1");
				numConfsUnprunedForPartition[i] = new BigInteger("1");
				//numInitConfsUnprunedForPartition[i] = new BigInteger("1");
				if (prunedRotAtRes[partIndex]){ //current partitioning rotamer already pruned, so no unpruned confs for this partition
					numConfsUnprunedForPartition[i] = new BigInteger("0");
					numPrunedPartitions++;
				}
				
				for (int j=0; j<numInAS; j++){
					if (!(isSplitRes(j,majorSplitPos,curDepth))){ //the split residues contribute only 1 rotamer
						numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]));
						numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
						//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numRotForRes[j]-numInitPrunedRotForRes[j]));
					}
				}
				if(ligPresent){
					numConfsTotalForPartition[i] = numConfsTotalForPartition[i].multiply(BigInteger.valueOf(numLigRotamers));
					numConfsUnprunedForPartition[i] = numConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
					//numInitConfsUnprunedForPartition[i] = numInitConfsUnprunedForPartition[i].multiply(BigInteger.valueOf(numLigRotamers-numInitPrunedLigRot));
				}			
	
				numTotalConfs = numTotalConfs.add(numConfsTotalForPartition[i]);
				numUnprunedConfs = numUnprunedConfs.add(numConfsUnprunedForPartition[i]);
				
				BigInteger pruneDiffFact = BigInteger.valueOf(10).pow(diffFact);
				
				System.out.println("Num unpruned confs: "+numConfsUnprunedForPartition[i]+" diffFact: "+pruneDiffFact);
				System.out.println();
				System.out.println();
				
				//if ((curDepth+1<maxDepth)&&(numConfsUnprunedForPartition[i].compareTo(numInitUnprunedConfs.divide(pruneDiffFact))==1)){ //not enough pruned, so partition at new depth
				if ((curDepth+1<maxDepth)&&(numConfsUnprunedForPartition[i].compareTo(pruneDiffFact)==1)){ //not enough pruned, so partition at new depth
					doDACS(numInAS, rs, numLigRotamers, residueMap, ligPresent, resDefault,
							rotFile, prunedForPartition[i], algOption, sfFile, useFlags, initEw, pruningE, maxDepth, curDepth+1, 
							logPS, majorSplitPos, numConfsUnprunedForPartition[i], diffFact, outputConfInfo, minRatioDiff,
							doMinimize, minPEM, maxPEM, distrDEE, ligNum, ligType, sParams, approxMinGMEC, 
							lambda, false, -1, null, minimizeBB, numMaxMut, scaleInt, maxIntScale, useEref, eRef);
				}
				else if (!prunedRotAtRes[partIndex]){ //if enough pruned or maxDepth partitioning reached, do the rotamer search
					
					bestScore = Math.min(bestScore,rs.getBestE().doubleValue());//best E for the partitions so far
					
					//Do the rotamer search
					rs.doAStarGMEC(outputConfInfo,true,doMinimize,
							numInAS,totalNumRotamers,rotamerIndexOffset,residueMap,
							resDefault,numMaxMut,initEw,bestScore,null,approxMinGMEC,lambda,minimizeBB,useEref,eRef);
					
					numEvaluatedConfs = numEvaluatedConfs.add(rs.numConfsEvaluated); //add the evaluated confs for this partition
					pruningE = (float)Math.min(pruningE,rs.getBestE().doubleValue());//update cutoff energy for MinBounds
				}
				
				//unprune the artificially pruned indices
				for (int j=0; j<indToUnprune.length; j++){
					if (indToUnprune[j]){
						prunedForPartition[i][indexMap[j]] = false;
					}
				}
				
				//output pruning info to file
				logPS.print("curDepth: "+curDepth+" curPartition: "+i+" majorSplitPos: "+majorSplitPos[curDepth]+" ");
				logPS.print(numConfsTotalForPartition[i]+" "+numInitUnprunedConfs+" "+numConfsUnprunedForPartition[i]);
				logPS.println();
				logPS.println();logPS.flush();
			}
		}
		
		if (!distrDACS){ //not distributed DACS
		
			System.out.println("numTotalConfs: "+numTotalConfs+"; numUnprunedConfs: "+numUnprunedConfs+"; numEvaluatedConfs: "+numEvaluatedConfs);
			for (int i=0; i<numPartitions; i++)System.out.print(numConfsTotalForPartition[i]+" ");System.out.println();
			for (int i=0; i<numPartitions; i++)System.out.print(numConfsUnprunedForPartition[i]+" ");
			System.out.println();
			System.out.println("Major splitting residue number: "+majorSplitPos[curDepth]);
			System.out.println("Number of partitions: "+numPartitions);
			System.out.println("Number of non-zero partitions: "+(numPartitions-numPrunedPartitions));
			System.out.println("Additional pruned rotamers: ");
			
			//count the number of partitions for which a rotamer is pruned (counting only the rotamers
			//	not pruned by Bounds or simple traditional DEE)
			int countNumPartitionsPrunedRot[] = new int[prunedRotAtRes.length];
			for (int i=0; i<prunedRotAtRes.length; i++){//for each rotamer
				countNumPartitionsPrunedRot[i] = 0;
				if (!prunedRotAtRes[i]){ //only if not pruned by the other two methods
					for (int j=0; j<numPartitions; j++){ //check for each partition
						if (prunedForPartition[j][i])
							countNumPartitionsPrunedRot[i]++;
					}
				}
				
				//output information
				if (countNumPartitionsPrunedRot[i]>0)
					System.out.println("index: "+i+"; num partitions in which pruned: "+countNumPartitionsPrunedRot[i]);
			}
		}
		else { //distributed DACS
			cObj.bestScore = cObj.bestScore.min(rs.getBestE());
		}
	}
	
	//Compute the number of unpruned conformations 
	private BigInteger compNumUnprunedConfs(int numInAS, int totalNumRotamers, boolean prunedRotAtRes[], int numLigRotamers,
			int numRotForRes[], boolean ligPresent) {
		
		int numPrunedRotForRes[] = new int[numInAS]; //after pruning for this partition
		for (int j=0; j<numInAS; j++){
			numPrunedRotForRes[j] = 0;
			for (int k=0; k<totalNumRotamers; k++){ //pruned rot are true and must be in the current set of allowed AA
				if (prunedRotAtRes[j*totalNumRotamers + k])
					numPrunedRotForRes[j]++;
			}
		}
		int numPrunedLigRot = 0;
		for (int k=0; k<numLigRotamers; k++){
			if (prunedRotAtRes[numInAS*totalNumRotamers + k])
				numPrunedLigRot++;
		}
		
		//count the total num confs and the num unpruned confs
		BigInteger numConfsTotal = new BigInteger("1");
		BigInteger numConfsUnpruned = new BigInteger("1");
		
		for (int j=0; j<numInAS; j++){
			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numRotForRes[j]));
			numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(numRotForRes[j]-numPrunedRotForRes[j]));
		}
		if(ligPresent){
			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numLigRotamers));
			numConfsUnpruned = numConfsUnpruned.multiply(BigInteger.valueOf(numLigRotamers-numPrunedLigRot));
		}	
		
		return numConfsUnpruned;
	}
	
	//Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position)
	private int chooseSplitPosRandom(int numInAS){		
		Random randNum = new Random();
		return randNum.nextInt(numInAS);
	}
	
	//Choose the splitting position (use only AS positions; the ligand is chosen not to be a splitting position);
	//	Choose the AS residue position with the smallest fraction of pruned rotamers (from MinBounds and simple MinDEE)
	private int chooseSplitPos(int numInAS, boolean prunedRotAtRes[], int numTotalRotamers, int numRotForRes[],
			int majorSplitPos[], int curDepth, double minRatioDiff){		
		
		final int minPartitions = 5; //the min number of rotamers that a splitting residue can have
		
		double pruneRatio[] = new double[numInAS];
		int minPos = -1;
		double minRatio = (double)Math.pow(10,38);
		for (int curRes=0; curRes<numInAS; curRes++){
			
			if (!(isSplitRes(curRes,majorSplitPos,curDepth))){
				if (numRotForRes[curRes]>=minPartitions){ //do not split at residues with very small number of rotamers
					int curPruned = 0;
					for (int curRot=0; curRot<numTotalRotamers; curRot++){
						if (prunedRotAtRes[curRes*numTotalRotamers+curRot]){ //cur rot is pruned (pruned rotamers are necessarily in the cur set of allowed AAs)
							curPruned++;
						}
					}
					pruneRatio[curRes] = (double)curPruned/numRotForRes[curRes];
					if (minRatio>=pruneRatio[curRes]){
						if ((minPos==-1)||(curRes<minPos)||(minRatio>=pruneRatio[curRes]+minRatioDiff)) {//preference to split at lower-numbered residues
							minRatio = pruneRatio[curRes];
							minPos = curRes;
						}
					}
				}
			}
		}
		
		if (minPos!=-1){
			//System.out.println("minPos: "+minPos);
			//for (int i=0;i<numInAS;i++)System.out.print(pruneRatio[i]+" ");System.out.println();
			return minPos;
		}
		else //if split position not hosen, choose randomly
			return chooseSplitPosRandom(numInAS);
	}
	
	//Check if the residue curRes is one of the splitRes
	private boolean isSplitRes(int curRes, int majorSplitPos[], int curDepth){
		for (int i=0; i<=curDepth; i++){
			if (curRes==majorSplitPos[i])
				return true;
		}
		return false;
	}
	
	//Compute the number of rotamers for each residue position (assign to numRotForRes[])
	private int [] compNumRotForRes(int numInAS, RotamerSearch rs, int numLigRot, int residueMap[]){
		
		int numRotForRes[] = new int[numInAS];
		boolean ligPresent = (numLigRot==0); //ligand present
		int treeLevels = numInAS;
		if (ligPresent)
			treeLevels++;
		
		numRotForRes = new int[treeLevels];
		
		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				curNumRot = numLigRot;
			}
			else { //AS residue				
				curNumRot = 0;
				for (int i=0; i<rs.sysLR.getNumAllowable(residueMap[curLevel]); i++){ //add the rot for all allowable AA at this residue
					int newRot = rl.getNumRotForAAtype(rs.sysLR.getIndexOfNthAllowable(residueMap[curLevel],i));
					if (newRot==0) //GLY or ALA
						newRot = 1;
					curNumRot += newRot; 
				}
			}
			numRotForRes[curLevel] = curNumRot;
		}
		return numRotForRes;
	}
	
	//Get the mapping between rotamer indices (into the pruning matrix) and the number of the
	//	current rotamer for the giveen residue; assumes sysLR in rs is valid (all allowables for the AS residues)
	private int [] getIndexMap(int numPartitions, RotamerSearch rs, int curRes, int residueMap[], int numTotalRot,
			int rotamerIndexOffset[]){
		
		int indexMap[] = new int[numPartitions];
		int indNum = 0;
		for (int AA=0; AA<rs.sysLR.getNumAllowable(residueMap[curRes]); AA++){ //for each AA for the given AS residue
			int curAA = rs.sysLR.getIndexOfNthAllowable(residueMap[curRes],AA);
			int numRotForAA = rl.getNumRotForAAtype(curAA);
			if (numRotForAA==0) //GLY or ALA
				numRotForAA = 1;
			
			for (int curRot=0; curRot<numRotForAA; curRot++){ //for each rot for the given AA
				indexMap[indNum] = curRes*numTotalRot + rotamerIndexOffset[curAA] + curRot;
				indNum++;
			}
		}
		return indexMap;
	}
	
	//Setup the file with name filename for output
	private PrintStream setupOutputFile(String fileName){
		PrintStream logPS = null; //the output file for conf info
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream(fileName);
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		return logPS;
	}	
///////////////////////////////////////////////////////////////////////////
//	End of DEE section
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//	Distributed DACS section
///////////////////////////////////////////////////////////////////////////
	//Handles the distribution of the DEE computation to the slave nodes
	private void doDistrDACSMaster(String runName, int numInAS, RotamerSearch rs, int numLigRotamers, int residueMap[],
			boolean ligPresent, String resDefault[], String rotFile,
			boolean prunedRotAtRes[], int algOption, String sfFile, boolean useFlags, float initEw, float pruningE,
			int maxDepth, int majorSplitPos[], BigInteger numInitUnprunedConfs,
			int diffFact, String outputPruneInfo, String outputConfInfo, double minRatioDiff, boolean doMinimize, String minPEM, 
			String maxPEM, int ligNum, String ligType, Hashtable sParams, 
			boolean approxMinGMEC, float lambda, int numRotForRes[], OneMutation resumeResults[], boolean minimizeBB, int numMaxMut,
			boolean scaleInt, float maxIntScale, boolean useEref, float eRef[]){
		
		System.out.println("Starting DACS (distributed)");
		
		int numPartitions = numRotForRes[majorSplitPos[0]]; //the number of partitions
		
		//map rotamer index to rot num for the splitting residue
		int indexMap[] = getIndexMap(numPartitions,rs,majorSplitPos[0],residueMap,totalNumRotamers,rotamerIndexOffset);
	
		//get only the non-zero partitions
		OneMutation mutArray[] = new OneMutation[numPartitions];
		int curMut = 0;
		for (int i=0; i<numPartitions; i++){
			if (!prunedRotAtRes[indexMap[i]]){
				mutArray[curMut] = new OneMutation();
				mutArray[curMut].mutNum = indexMap[i];
				curMut++;
			}
		}		
		
		OneMutation tmpArray[] = new OneMutation[curMut];
		System.arraycopy(mutArray, 0, tmpArray, 0, curMut);
		mutArray = tmpArray;
		
		System.out.println();
		System.out.println("Partitioning residue: "+majorSplitPos[0]);
		System.out.println("Number of partitions: "+numPartitions);
		System.out.println("Number of non-zero partitions: "+curMut);
		
		if (resumeResults!=null){ //remove completed partitions
			curMut = 0;
			OneMutation tmpArray2[] = new OneMutation[mutArray.length];
			boolean partFound = false;
			for (int i=0; i<mutArray.length; i++){
				partFound = false;
				for (int j=0; j<resumeResults.length; j++){
					if (mutArray[i].mutNum==resumeResults[j].mutNum){
						partFound = true;
						break;
					}
				}
				if (!partFound){
					tmpArray2[curMut] = new OneMutation();
					tmpArray2[curMut].mutNum = mutArray[i].mutNum;
					curMut++;
				}
			}
			tmpArray = new OneMutation[curMut];
			System.arraycopy(tmpArray2, 0, tmpArray, 0, curMut);
			mutArray = tmpArray;
			
			System.out.println("Number non-zero partitions after removing completed results: "+mutArray.length);
			
			//update the best energy so far
			for (int j=0; j<resumeResults.length; j++){
				pruningE = (float)Math.min(pruningE, resumeResults[j].score.doubleValue());
			}
		}
		
		System.out.println();
		
		
		//output the rs object
		outputObject(rs,rotFile);
		
		MutationManager mutMan = new MutationManager(runName,mutArray,false);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigNum(ligNum);
		mutMan.setLigType(ligType);
		mutMan.setarpFilenameMin(minPEM);
		mutMan.setarpFilenameMax(maxPEM);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.setIsLigAA(true);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setCalculateVolumes(false);
		mutMan.setNumMutations(numMaxMut);
		mutMan.setInitEw(initEw);
		mutMan.setPruningE(pruningE);
		mutMan.setLigPresent(ligPresent);
		mutMan.setUseSF(useFlags);
		mutMan.setPrunedRot(prunedRotAtRes);
		mutMan.setSfFile(sfFile);
		mutMan.setRotFile(rotFile);
		mutMan.setDistrDACS(true);
		mutMan.setDistrDEE(false);
		mutMan.setBestScore(new BigDecimal(pruningE)); //the best E initially is the pruningE read from the parameter file
		mutMan.setAlgOption(algOption);
		mutMan.setMaxDepth(maxDepth);
		mutMan.setDiffFact(diffFact);
		mutMan.setMinRatioDiff(minRatioDiff);
		mutMan.setNumInitUnprunedConf(numInitUnprunedConfs);
		mutMan.setOutputPruneInfo(outputPruneInfo);
		mutMan.setOutputConfInfo(outputConfInfo);
		mutMan.setMSP(majorSplitPos);
		mutMan.setApproxMinGMEC(approxMinGMEC);
		mutMan.setLambda(lambda);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setUseEref(useEref);
		mutMan.setEref(eRef);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
	}
	
	// Distributed DEE Slave function
	private CommucObj doDistrDACSSlave(CommucObj cObj) {
				
		long startTime = System.currentTimeMillis();
						
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);
		
		RotamerSearch rs = (RotamerSearch)readObject(cObj.rotFileIn); //load the saved rs from master
		
		String outputConfInfo = ("./conf_info/"+cObj.outputConfInfo+cObj.mutationNumber+"_"+cObj.partIndex);
		PrintStream logPS = setupOutputFile("./conf_info/"+cObj.outputPruneInfo+cObj.mutationNumber+"_"+cObj.partIndex);
		
		
		//Perform DACS
		doDACS(cObj.numInAS, rs, cObj.numLigRotamers, cObj.residueMap,
				cObj.ligPresent, cObj.resDefault, null,
				cObj.prunedRot, cObj.algOption, cObj.sfFileIn, cObj.useSF, cObj.initEw, cObj.pruningE,
				cObj.maxDepth, 0, logPS, cObj.msp, cObj.numInitUnprunedConf,
				cObj.diffFact, outputConfInfo, cObj.minRatioDiff, cObj.doMinimization, null, 
				null, false, cObj.ligNum, cObj.ligType, cObj.params,  
				cObj.approxMinGMEC, cObj.lambda, true, cObj.partIndex, cObj, cObj.minimizeBB, cObj.numMutations,
				cObj.scaleInt, cObj.maxIntScale, cObj.useEref, cObj.eRef);	
		
		
		logPS.flush();
		logPS.close();		
		
		rs = null;
		cObj.prunedRot = null;//smaller object, for sending back		
		
		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
		
		return cObj;
	}		
///////////////////////////////////////////////////////////////////////////
//	End of Distributed DACS section
///////////////////////////////////////////////////////////////////////////
	
	
///////////////////////////////////////////////////////////////////////////
//	Distributed DEE section
///////////////////////////////////////////////////////////////////////////
	//Handles the distribution of the DEE computation to the slave nodes
	private void doDistrDEEMaster(int numInAS, int numLigRotamers, int residueMap[],
			int ligNum, String ligType, Hashtable sParams, boolean ligPresent, String resDefault[],
			String minPEM, String maxPEM, float initEw, boolean prunedRotAtRes[],
			boolean doMinimize, RotamerSearch rs, String sfFile, String rotFile, boolean useSF, 
			int numSpPos, int msp[], int typeDEE, boolean minimizeBB, boolean scaleInt, float maxIntScale, boolean useEref, float eRef[]){
		
		//the total number of residues (active site + ligand, if present)
		int totalNumRes = numInAS;
		if (ligPresent) //ligand is present
			totalNumRes++;
		
		// Generate all combinations
		int numMutRes = 1; //singles
		if (typeDEE==optPairs)
			numMutRes = 2; //pairs criterion
		int numComb = factorial(totalNumRes).divide(factorial(totalNumRes-numMutRes).multiply(factorial(numMutRes))).intValue();
		int residueMutatable[][] = new int[numComb][totalNumRes];
		generateCombinations(residueMutatable,totalNumRes,numMutRes);
		
		OneMutation mutArray[] = new OneMutation[numComb];
		for (int curMut=0; curMut<mutArray.length; curMut++){
			mutArray[curMut] = new OneMutation();
			mutArray[curMut].mutNum = curMut;
			mutArray[curMut].resMut = new int[totalNumRes];
			for (int curRes=0; curRes<totalNumRes; curRes++)
				mutArray[curMut].resMut[curRes] = residueMutatable[curMut][curRes];
		}
		
		boolean splitFlags[][] = null;
		splitFlags = rs.getSplitFlags(splitFlags);//get the current DE pairs
		outputObject(splitFlags,sfFile);
		
		String AAallowed[] = new String[numInAS]; //the AA's to which each residue can mutate
		for (int i=0; i<numInAS; i++)
			AAallowed[i] = (String)sParams.get("RESALLOWED"+i);
		
		MutationManager mutMan = new MutationManager(null,mutArray,false);
		mutMan.setResidueMap(residueMap);
		mutMan.setResDefault(resDefault);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigNum(ligNum);
		mutMan.setLigType(ligType);
		mutMan.setarpFilenameMin(minPEM);
		mutMan.setarpFilenameMax(maxPEM);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.setnumLigRotamers(numLigRotamers);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.setIsLigAA(true);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setCalculateVolumes(false);
		mutMan.setInitEw(initEw);
		mutMan.setLigPresent(ligPresent);
		mutMan.setUseSF(useSF);
		mutMan.setSpFlags(splitFlags);
		mutMan.setPrunedRot(prunedRotAtRes);
		mutMan.setSfFile(sfFile);
		mutMan.setRotFile(rotFile);
		mutMan.setDistrDACS(false);
		mutMan.setDistrDEE(true);
		mutMan.setAAallowed(AAallowed);
		mutMan.setNumSpPos(numSpPos);
		mutMan.setMSP(msp);
		mutMan.setTypeDEE(typeDEE);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setScaleInt(scaleInt);
		mutMan.setMaxIntScale(maxIntScale);
		mutMan.setUseEref(useEref);
		mutMan.setEref(eRef);

		try{
			handleDoMPIMaster(mutMan,mutArray.length);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
	}
	
	// Distributed DEE Slave function
	private CommucObj doDistrDEESlave(CommucObj cObj) {
				
		long startTime = System.currentTimeMillis();
						
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,ligPresent,cObj.ligType);
		
		RotamerSearch rs = null;
		if (ligPresent)
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,true, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE, cObj.doSolvationE, cObj.solvScale, cObj.vdwMult, rl);
		else
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE, cObj.doSolvationE, cObj.solvScale, cObj.vdwMult, rl);
		
		
		System.out.print("Loading precomputed energy matrix...");

		rs.loadPairwiseEnergyMatrices(new String(cObj.arpFilenameMin+".dat"));
		if (cObj.doMinimization)
			rs.loadPairwiseEnergyMatricesMax(new String(cObj.arpFilenameMax+".dat"));
		System.out.println("done");
		
		if (cObj.useEref) { //add the reference energies to the min (and max) intra-energies
			rs.addEref(cObj.eRef, cObj.doMinimization, ligPresent, cObj.numInAS);
		}
		
		
		int totalNumRes = cObj.numInAS;
		if (ligPresent)
			totalNumRes++;				
		
		//determine the two residues in the pair (for pairs DE) or the one residue (split-DEE)
		boolean resInMut[] = new boolean[totalNumRes];
		for (int i=0; i<totalNumRes; i++){
			if (cObj.resMut[i]==1)
				resInMut[i] = true;
			else
				resInMut[i] = false;
		}
		
		
		//Set the allowable AA for each residue;
		// 		the ligand allowable set in the RotamerSearch() constructor				
		for (int j=0; j<cObj.numInAS; j++){
			String tempResAllow = (String)cObj.params.get("RESALLOWED"+j);
			for(int q=0;q<numTokens(tempResAllow);q++) {
				rs.setAllowable(cObj.residueMap[j],getToken(tempResAllow,q+1));
			}
			
			rs.setAllowable(cObj.residueMap[j],cObj.resDefault[j]); //the default type is set last
		}
	
		
		//Perform DEE pairs
		boolean splitFlags[][] = (boolean [][])readObject(cObj.sfFileIn);
		rs.setSplitFlags(splitFlags.length);
		rs.setSplitFlags(splitFlags);//initialize the split flags
		
		if (cObj.typeDEE==optPairs){ //simple Goldstein pairs	
			rs.DoDEEPairs(cObj.numInAS, cObj.numTotalRotamers, cObj.numLigRotamers, cObj.residueMap,
					rotamerIndexOffset, cObj.initEw, cObj.prunedRot, resInMut, cObj.doMinimization,
					cObj.useSF, false, true, cObj.minimizeBB, cObj.scaleInt, cObj.maxIntScale);
		}
		else if (cObj.typeDEE==optSplit){ //1- or 2-sp split-DEE
			//Precompute the MinDEE intervals
			rs.doCompMinDEEIntervals(cObj.numInAS, cObj.numTotalRotamers, cObj.numLigRotamers, cObj.residueMap, 
					rotamerIndexOffset, cObj.prunedRot, cObj.scaleInt, cObj.maxIntScale);
			
			cObj.prunedRot = rs.DoDEEConfSplitting(cObj.numInAS, cObj.numTotalRotamers, cObj.numLigRotamers, 
					cObj.residueMap, rotamerIndexOffset, cObj.initEw, cObj.prunedRot, resInMut, cObj.msp, -1, 
					cObj.doMinimization, true, cObj.numSpPos, true, cObj.minimizeBB);
		}
		
		long stopTime = System.currentTimeMillis();
		cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
		
		
		//We cannot send back the whole splitFlags matrix, since it is too big; even just sending the 
		//		newly-identified DE pairs can be infeasible for large systems, so we output the 
		//		new splitFlags matrix to a temp file, which is then read by the master node
		if (cObj.typeDEE==optPairs){
			cObj.sfFileOut = ("tmp_"+cObj.mutationNumber+"_"+stopTime);
			splitFlags = rs.getSplitFlags(splitFlags);
			outputObject(splitFlags,cObj.sfFileOut);
		}
		
		return cObj;
	}	
///////////////////////////////////////////////////////////////////////////
//	End of Distributed DEE section
///////////////////////////////////////////////////////////////////////////
	
	private Object readObject(String inFile){
		return readObject(inFile,true);
	}
	
	private Object readObject(String inFile, boolean repeat){
		Object inObj = null;
		boolean done = false;
		while (!done){
			try{
				ObjectInputStream in = new ObjectInputStream(new FileInputStream(inFile));
				inObj = in.readObject();
				in.close();
				done = true;
			}
			catch (Exception e){
				//System.out.println(e.toString());
				//System.out.println("ERROR: An exception occurred while reading from object file");
				if (repeat)
					done = false;
				else
					done = true;
			}
		}
		return inObj;
	}
	
	private void outputObject(Object outObj, String outFile){
		try{
			FileOutputStream fout = new FileOutputStream(outFile);
			ObjectOutputStream out = new ObjectOutputStream(fout);
			out.writeObject(outObj);
			out.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while writing object file");
			System.exit(0);
		}
	}
	
	private int countPrunedRot(boolean prunedRot[]){
		int countPruned = 0;
		for (int i=0; i<prunedRot.length; i++){
			if (prunedRot[i])
				countPruned++;
		}
		return countPruned;
	}
	
	private int countPrunedPairs(boolean prunedPairs[][]){
		int countPruned = 0;
		for (int i=0; i<prunedPairs.length; i++){
			for (int j=i+1; j<prunedPairs.length; j++){
				if (prunedPairs[i][j])
					countPruned++;
			}
		}
		return countPruned;
	}
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// MPI section
//		The tutorial "One-step Tutorial: MPI: It's easy to get started" 
//			(http://www.lam-mpi.org/tutorials/one-step/ezstart.php ; accessed Oct 23, 2006) was used as MPI code reference
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Setup and do MPI
	//	The parameter is only used by the master node
	public void handleDoMPI(String args[]) throws MPIException {
		
		MPI.Init(args);
		
		int procRank = MPI.COMM_WORLD.Rank();
		numProc = MPI.COMM_WORLD.Size();	
		
		System.out.println("Node rank: "+procRank+" of "+numProc);
		MPI.COMM_WORLD.Barrier();
		
		setConfigPars();
		MPI.COMM_WORLD.Barrier();
		
		if (procRank==0){ //master node
			
			outputProgInfo();
			parse(args);
		}
		else {//slave node
			handleDoMPISlave();
		}
		
		MPI.Finalize();
	}
	
	//Do MPI for the master node
	public void handleDoMPIMaster(MutationManager mutMan, int size) throws MPIException {
		
		CommucObj cObjArray[] = new CommucObj[size];
		int numFinished = 0;
		
		int curMut = 0;
		for (int curProc=1; curProc<numProc; curProc++){ //distribute a single mutation per processor, for all processors
			
			if (curMut<cObjArray.length){ //more mutations to distribute
				
				System.out.println("Retrieving "+curMut+" of "+(cObjArray.length));
				cObjArray[curMut] = mutMan.getNextComObj(curMut);
				
				MPI.COMM_WORLD.Send(cObjArray, curMut, 1, MPI.OBJECT, curProc, regTag);
				curMut++;
				
				System.out.println("Sent to proc "+curProc);
				System.out.println();
			}
			else
				break;
		}
		
		while (numFinished<cObjArray.length){ //distribute and receive all remaining mutations
			
			CommucObj cObj[] = new CommucObj[1];
			cObj[0] = new CommucObj();
			
			Status s = MPI.COMM_WORLD.Recv(cObj, 0, 1, MPI.OBJECT, MPI.ANY_SOURCE, regTag);
			mutMan.processFinishedMutation(cObj[0]);
			numFinished++;
			
			System.out.println("Finished: "+cObj[0].mutationNumber+", Time: "+(cObj[0].elapsedTime/60.0));
			
			if (curMut<cObjArray.length){
				
				System.out.print("Retrieving "+curMut+" of "+(cObjArray.length));
				cObjArray[curMut] = mutMan.getNextComObj(curMut);
				
				MPI.COMM_WORLD.Send(cObjArray, curMut, 1, MPI.OBJECT, s.source, regTag);
				curMut++;
				
				System.out.println(", Sent to proc "+s.source);
				System.out.println();
			}
			/*else { //no more jobs to distribute, so send a done signal to the current node
				cObj[0] = new CommucObj();
				MPI.COMM_WORLD.Send(cObj, 0, 1, MPI.OBJECT, s.source, killTag);
			}*/
		}
	}
	
	//Do MPI for a slave node
	public void handleDoMPISlave() throws MPIException {
		
		while (true){
			
			CommucObj cObj[] = new CommucObj[1];
			cObj[0] = new CommucObj();
			
			Status s = MPI.COMM_WORLD.Recv(cObj, 0, 1, MPI.OBJECT, 0, MPI.ANY_TAG);
			
			if (s.tag==killTag) //done
				return;
			
			else { //new computation
				cObj[0] = handleKSSlave(cObj[0]); //perform computation
				MPI.COMM_WORLD.Send(cObj, 0, 1, MPI.OBJECT, 0, regTag); //send back result
			}
		}
	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	 End of MPI section
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////////
//	Backbone Flexibility Section
///////////////////////////////////////////////////////////////////////////
	public void generateBackbones(String s, Hashtable curScope){
		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Backbone config filename (string)
		
		// The system paramaters
		Hashtable sParams = new Hashtable(hashSize);

		System.out.println("Performing Backbone Generation");

		// Read System Parameters
		readConfigFile(sParams,getToken(s,2));
		readConfigFile(sParams,getToken(s,3));
		
		// Pull search parameters
		int numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();	
		String runName = ((String)sParams.get("RUNNAME"));
		boolean sysSampling = (new Boolean((String)sParams.get("SYSSAMPLING"))).booleanValue();
		double theta = (new Float((String)sParams.get("THETA"))).floatValue();
		double alpha = (new Float((String)sParams.get("ALPHA"))).floatValue();
		int numSamples = (new Integer((String)sParams.get("NUMSAMPLES"))).intValue();
		boolean ligPresent = (new Boolean((String)sParams.get("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.get("LIGTYPE"));
		
		if (theta%alpha!=0){
			System.out.println("ERROR: Choose theta = k*alpha, for k - an integer.");
			System.exit(1);
		}
			
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]+"("+resDefault[i]+")");
		}
		System.out.println();
		
		RotamerSearch rs = null;
		if (ligPresent)
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,true, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);
		else
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);
		
		rs.doGenBackbones(runName, numInAS, residueMap, theta, alpha, numSamples, sysSampling);
	}
///////////////////////////////////////////////////////////////////////////
//	End of Backbone Flexibility Section
///////////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////////
//	Generate Random Conformations Section
///////////////////////////////////////////////////////////////////////////
	//Generates a random set of mutations/conformations for a given system
	public void generateRandConfs(String s, Hashtable curScope) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Run name (for output files)
		// 3: Ligand is present (boolean)
		// 4: Ligand type (if present)
		// 5: Number of conformations to be generated

		int numInAS = -1;

		// The system paramaters
		Hashtable sParams = new Hashtable(hashSize);

		// Read System Parameters
		readConfigFile(sParams,getToken(s,2));

		// Pull search parameters
		numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		String runName = getToken(s,3);
		boolean ligPresent = parseToBoolean(getToken(s,4),curScope);
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,5);
		int num = parseToInt(getToken(s,6),curScope);
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);

		int residueMap[] = new int[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			System.out.print(" "+residueMap[i]);
		}
		System.out.println();
		
		StrandRotamers sysLR = null;
		StrandRotamers ligLR2 = null;
		
		sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
		if (ligPresent)
			ligLR2 = new StrandRotamers(rl,m.strand[ligStrNum]);
		
		
		PrintStream logPS = setupOutputFile(runName);
		
		Random r = new Random();
		for (int i=0; i<num; i++){
			String AAs[] = new String[numInAS];
			int rot[] = new int[numInAS+1];
			for (int a=0; a<numInAS; a++){
				AAs[a] = rl.getAAName(r.nextInt(numAAallowed));
				int n = rl.getNumRotamers(AAs[a]);
				if (n<=0)
					n = 1;
				rot[a] = r.nextInt(n);
			}
			if (ligPresent){
				int n = rl.getNumRotamers(ligType);
				if (n<=0)
					n = 1;
				rot[numInAS] = r.nextInt(n);
			}
			
			logPS.print(i+" ");
			for (int a=0; a<numInAS; a++)
				logPS.print(AAs[a]+" ");
			if (ligPresent)
				logPS.print(ligType+" ");
			
			for (int a=0; a<numInAS; a++)
				logPS.print(rot[a]+" ");
			if (ligPresent)
				logPS.print(rot[numInAS]+" ");
			
			logPS.println();
			logPS.flush();
		}
		logPS.close();
	}
///////////////////////////////////////////////////////////////////////////
//	End of Generate Random Conformations Section
///////////////////////////////////////////////////////////////////////////
	
	public void doSinglePairE(String s, Hashtable curScope, Hashtable sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new Hashtable(hashSize);
	
			// Read System parameters
			readConfigFile(sParams,getToken(s,2));
			
			// Read Mutation search parameters
			readConfigFile(sParams,getToken(s,3));
		}
		
		int numInAS = (new Integer((String)sParams.get("NUMINAS"))).intValue();
		int numMutations = 2; //pairwise energies are computed
		boolean minimizeBB = (new Boolean((String)sParams.get("MINIMIZEBB"))).booleanValue();
		String runName = (String)sParams.get("RUNNAME");
		String minEMatrixName = (String)sParams.get("MINENERGYMATRIXNAME");
		String maxEMatrixName = (String)sParams.get("MAXENERGYMATRIXNAME");
		
		boolean ligPresent = (new Boolean((String)sParams.get("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)sParams.get("LIGTYPE");

		boolean resumeSearch = (new Boolean((String)sParams.get("RESUMESEARCH"))).booleanValue();
		String resumeFilename = (String)sParams.get("RESUMEFILENAME");
		
		System.out.println("Run Name: "+runName);
		System.out.println("Precomputed Minimum Energy Matrix: "+minEMatrixName);
		System.out.println("Precomputed Maximum Energy Matrix: "+maxEMatrixName);
		System.out.println("Ligand Type: "+ligType);
		System.out.println("Num Residues Allowed to Mutate: "+numMutations);

		
		System.out.println("Computing _All_ Rotamer-Rotamer Energies");
		
		System.out.println("Starting minimum and maximum bound energy computation");
		
		if(resumeSearch) {
			System.out.println("** Resuming Search **");
			System.out.println("     resuming from file: "+resumeFilename);
		}
		
		//Setup the molecule system
		Molecule m = new Molecule();
		int numLigRotamers = setupMolSystem(m,sParams,ligPresent,ligType);
					
		int residueMap[] = new int[numInAS];
		String resDefault[] = new String[numInAS];
		String resMapString = (String)sParams.get("RESIDUEMAP");
		System.out.print("ResidueMap:");
		for(int i=0;i<numInAS;i++){
			residueMap[i] = (new Integer(getToken(resMapString,i+1))).intValue();
			resDefault[i] = m.strand[sysStrNum].residue[residueMap[i]].name;
			System.out.print(" "+residueMap[i]);
		}
		System.out.println();
		
		RotamerSearch rs = null;
		if (ligPresent)
			rs = new RotamerSearch(m,sysStrNum,ligStrNum,true, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);
		else
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);
		
		
		
		
		int resMut[] = new int[numInAS];
		for (int i=0; i<numInAS; i++)
			resMut[i] = 0;
		
		String flagMutType = "AS-AS";
		resMut[0] = 1; resMut[1] = 1;
		
		

		if (((flagMutType.compareTo("AS-AS")==0)||(flagMutType.compareTo("SHL-AS")==0)||(flagMutType.compareTo("TEMPL")==0))&&(ligPresent)){
			m.deleteStrand(1);//we do not need the ligand for these runs
			rs = null;
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE, doSolvationE, solvScale, softvdwMultiplier, rl);					
		}
		
		System.out.println("Beginning setAllowables");
		// Ligand allowable set in the RotamerSearch() constructor		
		for(int j=0; j<numInAS; j++) {
			/*if (resMut[j] == 1) {
				for(int q=0;q<resAllowed.length;q++) {
					rs.setAllowable(residueMap[j],resAllowed[q]);
				}
			}*/
			//else {
				rs.setAllowable(residueMap[j],resDefault[j]);
			//}
		}	
		
		// The goal is that the total energy of a system can be bounded by the sum of 
		//  all pairwise active site residue entries plus the entry for each active site
		//  residue's shell run plus each active site residue's self intra-residue energy.
		//  If a ligand is present then one should add the ligand to shell energy, the
		//  ligand to each active site residue pairwise energy, and the ligand self intra-
		//  residue energy.
		
		//Create a matrix that will contain both the min and max pairwise energy matrices;
		//	the min matrix is in the first half of the columns, the max is in the second half:
		//		minMaxEMatrix[i][j] = minMatrix[i][j]
		//		minMaxEMatrix[i][mutEMatrixSize + j] = maxMatrix[i][j]
		int mutEnerMatrixSize = totalNumRotamers*numInAS+numLigRotamers+1; // +1 for the backbone				
		float minMaxEMatrix[][] = new float[mutEnerMatrixSize][mutEnerMatrixSize * 2];
		for(int i=0; i<mutEnerMatrixSize; i++) {
			for(int j=0; j<mutEnerMatrixSize*2; j++){
				minMaxEMatrix[i][j] = 0.0f;
			}
		}
		
		boolean shellRun = false;
		boolean intraRun = false;
		boolean templateOnly = false;
		
		if (flagMutType.compareTo("TEMPL")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = true;
		}
		else if (flagMutType.compareTo("AS-AS")==0){
			
			// **** Normal Active Site residue runs ****
			// Computes active site residue to active site residue pair energies					
			shellRun = false;ligPresent = false;intraRun = false;templateOnly = false;
		}	
		else if (flagMutType.compareTo("SHL-AS")==0){
			
			// Then shell runs for the active site residues
			// Computes the active site residue rotamers to shell energies					
			shellRun = true;ligPresent = false;intraRun = false;templateOnly = false;					
		}
		else if (flagMutType.compareTo("INTRA")==0){
			
			// Compute all intra-residue energies					
			shellRun = false;intraRun = true;templateOnly = false;
		}				
		else if (flagMutType.compareTo("LIG-AS")==0){
			
			// **** Ligand present runs ****
			// This section computes the inter-residue energies between
			//  active site residues and the ligand
			shellRun = false;intraRun = false;templateOnly = false;					
		}
		else { //(cObj.flagMutType.compareTo("LIG-SHL")==0)
			
			// Computes ligand rotamer to shell energies
			shellRun = true; intraRun = false;templateOnly = false;
		}
		
		//Compute the corresponding matrix entries
		minMaxEMatrix = rs.simplePairwiseMutationAllRotamerSearch(residueMap,
				numInAS,resAllowed,resAllowed.length,rotamerIndexOffset,
				totalNumRotamers,true,ligPresent,shellRun,intraRun,
				resMut,minMaxEMatrix,mutEnerMatrixSize,minimizeBB,templateOnly);
		
	}
	
	//Computes conformation energies for different combinations of the energy function parameters
	private void fitEparams(String s, Hashtable curScope){
		
		String firstParam = getToken(s,1);
		
		String sysFile = getToken(s,2);

		// Pull search parameters
		String confResFile = getToken(s,3);
		String runName = getToken(s,4);
		boolean ligPresent = parseToBoolean(getToken(s,5),curScope);
		String ligType = null;
		if (ligPresent)
			ligType = getToken(s,6);
		int numResults = parseToInt(getToken(s,7),curScope);
		boolean minimizeBB = parseToBoolean(getToken(s,8),curScope);
		
		int numSteps = 10;
		double maxVdwMult = 1.05;
		double maxSolvScale = 0.3;
		
		solvScale = 0.0;
		double initVdwMult = 0.63;
				
		double vdwDelta = (maxVdwMult-initVdwMult)/numSteps;
		double solvDelta = (maxSolvScale-solvScale)/numSteps;
		
		for (int i1=0; i1<numSteps; i1++){
			solvScale += solvDelta;
			for (int i2=0; i2<2; i2++){
				if (i2==0)
					distDepDielect = true;
				else
					distDepDielect = false;
				
				for (int i3=0; i3<=numSteps; i3++){
					if (i3==0)
						dielectConst = 1.0;
					else
						dielectConst = 4*i3;
					
					softvdwMultiplier = initVdwMult-vdwDelta;
					for (int i4=0; i4<=numSteps; i4++){
						softvdwMultiplier += vdwDelta;
						
						String runNameParams = (runName+"_"+solvScale+"_"+distDepDielect+"_"+dielectConst+"_"+softvdwMultiplier);
						
						String s1 = (firstParam+" "+sysFile+" "+confResFile+" "+runNameParams+" false none "+numResults+" "+minimizeBB);
						handleMinDEEApplyRot(s1,curScope);
						
						if (ligPresent){
							runNameParams = (runNameParams+"_lig");
							s1 = (firstParam+" "+sysFile+" "+(confResFile+"_lig")+" "+runNameParams+" true "+ligType+" "+numResults+" "+minimizeBB);
							handleMinDEEApplyRot(s1,curScope);
						}
					}
				}
			}
		}
	}
	
	
////////////////////////////////////////////////////////////////
// Compute Residue Entropy Section
////////////////////////////////////////////////////////////////
	public void handleDoResEntropy(String s, Hashtable curScope, Hashtable sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new Hashtable(hashSize);
	
			// Read System parameters
			readConfigFile(sParams,getToken(s,2));
			
			// Read Mutation search parameters
			readConfigFile(sParams,getToken(s,3));
		}
		
		String runName = (String)sParams.get("RUNNAME");
		
		String matrixName = (String)sParams.get("MATRIXNAME");
		boolean usePDBstats = (new Boolean((String)sParams.get("USEPDBSTATS"))).booleanValue();
		String pdbStatsFile = (String)sParams.get("PDBSTATSFILE");
		boolean useEref = (new Boolean((String)sParams.get("USEEREF"))).booleanValue();
		float dist = (new Float((String)sParams.get("DIST"))).floatValue();
		String rotProbFile = (String)sParams.get("ROTPROBFILE");
		float stericE = (new Float((String)sParams.get("STERICE"))).floatValue();
		float maxPairE = (new Float((String)sParams.get("MAXPAIRE"))).floatValue();
		
		
		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		
		int numRes = m.strand[sysStrNum].numberOfResidues;
		
		String resDefault[] = new String[numRes];
		int defResNum[] = new int[numRes];
		for(int i=0;i<numRes;i++){
			resDefault[i] = m.strand[sysStrNum].residue[i].name;
			defResNum[i] = m.strand[sysStrNum].residue[i].getResNumber();
		}
		
		
		rotProbFile = (rotProbFile+".dat");
		float rotProb[][] = (float [][])readObject(rotProbFile,false);
		if (rotProb==null){ //Perform SCMF to compute the rotamer probabilities		
			
			//read in or compute all of the energy matrices
			float mutationEnergiesMin[][] = getResEntropyEmatricesTempl(matrixName, sParams, numRes, resDefault, m, runName, dist);
			float asasE[][][][] = getResEntropyEmatricesPair(matrixName, sParams, numRes, resDefault, m, runName, dist);
			float intraEnergies[][] = getResEntropyEmatricesIntra(matrixName, sParams, numRes, m, runName);
			float eRef[] = getResEntropyEmatricesEref(sParams, useEref);
			
			rotProb = compRotProbSCMF(numRes, mutationEnergiesMin, intraEnergies, asasE, eRef, rotProbFile, resDefault, stericE, maxPairE);
		}
		
		int numProx[] = new int[numRes]; //get the number of proximate residues for each residue position
		for (int i=0; i<numProx.length; i++)
			numProx[i] = 0;
		String asDistFile = matrixName+"_dist.dat";
		boolean as[][] = (boolean [][])readObject(asDistFile,false);
		for (int i=0; i<numRes; i++){
			for (int j=i+1; j<numRes; j++){
				if (as[i][j]){
					numProx[i]++;
					numProx[j]++;
				}
			}
		}
		
		double aaProbPDB[][] = null;
		if (usePDBstats)
			aaProbPDB = getPDBstats(pdbStatsFile,m,numRes);	
		
		m = null;		
		
		PrintStream logPS = setupOutputFile(runName);
		
		logPS.print("resNum pdbResNum resDefault entropy"+" ");
		for (int j=0; j<resAllowed.length; j++)
			logPS.print(resAllowed[j]+" ");
		logPS.println("numProx");
		logPS.flush();
		
		
		//Compute the AA probabilities for each residue position (as a function of the rotamer probabilities);
		//Compute the entropy at each position as a function of the amino acid probabilities for that position
		final double kB = 1.0;
		for (int i=0; i<numRes; i++){
			
			if ( !resDefault[i].equalsIgnoreCase("PRO") ){
			
				float aaProbBBE[] = new float[resAllowed.length];
				int curInd = 0;
				for (int j=0; j<aaProbBBE.length; j++){
					
					int numCurRot = rl.getNumRotamers(resAllowed[j]);
					if (numCurRot==0)
						numCurRot = 1;
					
					aaProbBBE[j] = 0.0f;
					for (int r=0; r<numCurRot; r++){
						aaProbBBE[j] += rotProb[i][curInd];
						curInd++;
					}
				}
				
				//Compute the unnormalized AA probabilities as a weighted function of energies and PDB stats (if included)
				double aaProbUnNorm[] = new double[resAllowed.length];
				double aaNorm = 0.0;
				for (int j=0; j<resAllowed.length; j++){
					
					aaProbUnNorm[j] = 1.0;
					aaProbUnNorm[j] *= aaProbBBE[j];
					if (usePDBstats)
						aaProbUnNorm[j] *= aaProbPDB[i][j];
					
					aaNorm += aaProbUnNorm[j];
				}
				
				//Normalize the probabilities
				double aaProb[] = new double[resAllowed.length];
				for (int j=0; j<resAllowed.length; j++){
					if (aaNorm!=0.0)
						aaProb[j] = aaProbUnNorm[j]/aaNorm;
					else
						aaProb[j] = 0.0;
				}
				
				//Compute the entropy for the current residue position
				double sumAA = 0.0;
				for (int j=0; j<aaProb.length; j++){
					if (aaProb[j]>0.0)
						sumAA += aaProb[j]*Math.log(aaProb[j]);
				}
				
				double entropy = -kB * sumAA;
				
				logPS.print(i+" "+defResNum[i]+" "+resDefault[i]+" "+entropy+" ");
				for (int j=0; j<aaProb.length; j++)
					logPS.print(aaProb[j]+" ");
				
				logPS.println(numProx[i]);
				logPS.flush();
			}
			else {
				logPS.println(i+" "+defResNum[i]+" "+resDefault[i]+" "+0.0); //only for residue positions with wildtype Pro
			}
		}
		logPS.close();
	}
	
	//Computes the rotamer probabilities for all rotamers at all residue positions using SCMF
	private float[][] compRotProbSCMF(int numRes, float mutationEnergiesMin[][], float intraEnergies[][],
			float asasE[][][][], float eRef[], String rotProbFile, String resDefault[], float stericE, float maxPairE){
		
		final float constR = (float)(1.9891/1000.0);//the gas constant
		float T = 50000; //initial temperature
		final float endT = 298.15f; //the minimum temperature for annealing
		float tStepSize = 100.0f; //the temperature step size for annealing
		final float eps = 0.0001f; //the convergence threshold
		final float lambda = 0.5f; //scaling factor for updating the rotamer probabilities
		
		
		for (int i=0; i<asasE.length; i++){//Set the max energy for any element in asasE[][][][] to maxPairE
			if (asasE[i]!=null){
				for (int j=0; j<asasE[i].length; j++){
					if (asasE[i][j]!=null){
						for (int k=0; k<asasE[i][j].length; k++){
							if (asasE[i][j][k]!=null){
								for (int l=0; l<asasE[i][j][k].length; l++){
									if (asasE[i][j][k][l]>maxPairE)
										asasE[i][j][k][l] = maxPairE;
								}
							}
						}
					}
				}
			}
		}
		
		int numPrunedRot = 0;
		boolean prunedRot[][] = new boolean[numRes][totalNumRotamers];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<totalNumRotamers; j++){
				if ( (mutationEnergiesMin[0][1+i*totalNumRotamers+j] + intraEnergies[1+i*totalNumRotamers+j][0]) > stericE){
					prunedRot[i][j] = true;
					numPrunedRot++;
				}
				else
					prunedRot[i][j] = false;
			}
		}
		System.out.println("Num rotamers pruned due to incompatibility with the template: "+numPrunedRot);
		
		
		//For each residue, compute the probability of each rotamer for that residue
		float Emf[][] = new float[numRes][totalNumRotamers];
		float rotProb[][] = new float[numRes][totalNumRotamers];
		float oldProb[][] = new float[numRes][totalNumRotamers];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<totalNumRotamers; j++) {
				if (!prunedRot[i][j])
					rotProb[i][j] = 1.0f/totalNumRotamers;
				else
					rotProb[i][j] = 0.0f;
				
				oldProb[i][j] = rotProb[i][j];
			}
		}
		
		while (T>=endT){ //perform annealing
			
			System.out.println("Starting run at T = "+T);
			
			boolean done = false;
			while (!done){
				
				//Compute the new mean-field energy for each rotamer
				for (int i=0; i<numRes; i++){
					
					if ( !resDefault[i].equalsIgnoreCase("PRO") ){
					
						for (int j=0; j<totalNumRotamers; j++){
							
							if (!prunedRot[i][j]){
								
								Emf[i][j] = mutationEnergiesMin[0][1+i*totalNumRotamers+j] + intraEnergies[1+i*totalNumRotamers+j][0] - eRef[getAAindFromRotNum(j)];
								
								if (asasE[i]!=null){
									
									for (int k=0; k<numRes; k++){
										if (k!=i){ //for all residues with which i_j has contact
											if ( (i<k) && (asasE[i][k]!=null) ){
												for (int l=0; l<totalNumRotamers; l++){
													if (!prunedRot[k][l])
														Emf[i][j] += asasE[i][k][j][l]*rotProb[k][l];
												}
											}
											else if ( (i>k) && (asasE[k][i]!=null) ){
												for (int l=0; l<totalNumRotamers; l++){
													if (!prunedRot[k][l])
														Emf[i][j] += asasE[k][i][l][j]*rotProb[k][l];
												}
											}
										}
									}
								}
							}
						}
					}
				}
				
				//Update the rotamer probabilities			
				for (int i=0; i<numRes; i++){
					
					if ( !resDefault[i].equalsIgnoreCase("PRO") ){
					
						float normFactor = 0.0f;
						for (int j=0; j<totalNumRotamers; j++){
							
							if (!prunedRot[i][j])
								normFactor += (float)Math.exp( -Emf[i][j] / (constR*T));	
						}
						
						for (int j=0; j<totalNumRotamers; j++){
							
							if (!prunedRot[i][j]){
								
								oldProb[i][j] = rotProb[i][j]; //the probability before the update
								
								if (normFactor!=0.0f)
									rotProb[i][j] = lambda*((float)Math.exp( -Emf[i][j] / (constR*T)) / normFactor) + (1-lambda)*oldProb[i][j];
								else
									rotProb[i][j] = 0.0f;
							}
						}
					}
				}
				
				float rms = checkRotProbConvEntropy(rotProb,oldProb,resDefault,prunedRot);
				
				if (rms>eps)
					done = false;
				else
					done = true;
			}
			
			T -= tStepSize;
		}
		
		outputObject(rotProb,rotProbFile);
		
		return rotProb;
	}
	
	//Checks if the rotamer probabilities for the entropy computation have converged
	private float checkRotProbConvEntropy(float rotProb[][], float oldProb[][], String resDefault[], boolean prunedRot[][]){
		
		float sum = 0.0f;
		for (int i=0; i<rotProb.length; i++){
			if ( !resDefault[i].equalsIgnoreCase("PRO") ){
				for (int j=0; j<rotProb[i].length; j++){
					if (!prunedRot[i][j])
						sum += (float)Math.pow( (rotProb[i][j]-oldProb[i][j]) , 2.0);
				}
			}
		}
		float rms = (float)Math.sqrt(sum);
		
		System.out.println("RMS: "+rms);
		
		return rms;
	}
	
	//Reads in (if computed) or computes the energy matrices for rot-to-template energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private float [][] getResEntropyEmatricesTempl(String matrixName, Hashtable sParams, int numRes, 
			String resDefault[], Molecule m, String runName, float dist){
		
		String origPDB = (String)sParams.get("PDBNAME");
		int origLigNum = new Integer((String)sParams.get("LIGNUM")).intValue();
		int origNumCof = (new Integer((String)sParams.get("COFACTORNUM"))).intValue();		
		String origCofMapString = (String)sParams.get("COFMAP");
		int cofactorRes[] = null;
		if (origNumCof>0){
			cofactorRes = new int[origNumCof];
			for(int i=0;i<origNumCof;i++)
				cofactorRes[i] = (new Integer(getToken(origCofMapString,i+1))).intValue();
		}
		
		//Check for the rot-to-backbone file
		String bbName = matrixName+"_shll.dat";
		float mutationEnergiesMin[][] = (float [][])readObject(bbName,false); //check if already computed
		
		if (mutationEnergiesMin==null){ //compute the rot-to-template energy matrix
			
			String allGlyFile = (origPDB+".allgly");
			sParams.put("PDBNAME", allGlyFile);
			sParams.put("LIGNUM","-1");
			if (origLigNum>=0){
				String cm = "";
				for(int i=0;i<origNumCof;i++)
					cm += (cofactorRes[i]-1+" ");
				sParams.put("COFMAP", cm);
			}
			
			RotamerSearch rs = null;
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
						true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE,doSolvationE,solvScale,softvdwMultiplier, rl);
			
			Molecule m1 = new Molecule();
			if (setupMolSystem(m1,sParams,false,null)<0){				
				
				for(int i=0;i<numRes;i++){ //change all residues to GLY here for faster evaluation by the slave nodes
					if ( !resDefault[i].equalsIgnoreCase("GLY") && !resDefault[i].equalsIgnoreCase("PRO") )
						rs.sysLR.changeResidueType(m,i,"gly",true,true);
				}
				saveMolecule(m,allGlyFile,0.0f);
			}
			
			else {
				m = m1;
				m1 = null;
			}
			
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
					true, 0.0f, stericThresh, distDepDielect, dielectConst, doDihedE,doSolvationE,solvScale,softvdwMultiplier, rl);
			
			computeEntropyEmatrixMaster(numRes,runName+".log",bbName,sParams,false,false,
					true,mutationEnergiesMin,null,false,dist,null,false,false); 
			
			mutationEnergiesMin = (float [][])readObject(bbName,false);
		}
		
		sParams.put("PDBNAME", origPDB);
		sParams.put("LIGNUM", ""+origLigNum);
		sParams.put("COFMAP",origCofMapString);
		m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		
		if (ligStrNum>=0) //the ligand is not used here
			m.deleteStrand(ligStrNum);
		
		return mutationEnergiesMin;
	}
	
	//Reads in (if computed) or computes the energy matrices for rot-to-rot pairwise energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private float [][][][] getResEntropyEmatricesPair(String matrixName, Hashtable sParams, int numRes, String resDefault[], Molecule m,
			String runName, float dist){
		
		String origPDB = (String)sParams.get("PDBNAME");
		
		//Check for the pairwise rot-to-rot file;
		String pairName = matrixName+"_pair.dat";
		
		float asasE[][][][] = readPairMatrixEntropy(pairName,numRes);
		if (asasE==null){ //compute the rot-to-template energy matrix
			
			//For each residue position i, get all residue positions that are within dist
			String asDistFile = matrixName+"_dist.dat";
			boolean as[][] = (boolean [][])readObject(asDistFile,false);
			if (as==null){
				as = new boolean[numRes][numRes];
				computeEntropyEmatrixMaster(numRes,runName+".log",asDistFile,sParams,false,false,
						false,null,null,true,dist,as,false,false); //the distances are returned in as[][]
			}
			
			int numPairs = 0;
			asasE = new float[numRes][numRes][][];
			for (int i=0; i<numRes; i++){
				for (int j=i+1; j<numRes; j++){
					if (as[i][j]){
						asasE[i][j] = new float[1][];
						numPairs++;
					}
				}
			}
						
			computeEntropyEmatrixMaster(numRes,runName+".log",pairName,sParams,false,false,
					false,null,asasE,false,0.0f,null,false,false);
			
			sParams.put("PDBNAME", origPDB);
			m = new Molecule();
			setupMolSystem(m,sParams,false,null);
			
			if (ligStrNum>=0) //the ligand is not used here
				m.deleteStrand(ligStrNum);
			
			asasE = readPairMatrixEntropy(pairName,numRes);
		}		
		
		return asasE;
	}
	
	//Reads in (if computed) or computes the energy matrices for intra-rot energies;
	//This is not setup to handle energy minimization properly (mainly due to the approach used for computing only a small
	//		subset of the pairwise rot-to-rot energies), so doMinimize and minimizeBB should be false
	private float [][] getResEntropyEmatricesIntra(String matrixName, Hashtable sParams, int numRes, Molecule m, String runName){
		
		String origPDB = (String)sParams.get("PDBNAME");
		
		//Check for the intra energies file
		String intraName = matrixName+"_intra.dat";
		
		float intraEnergies[][] = (float [][])readObject(intraName,false); //check if already computed			
		
		if (intraEnergies==null) //compute the intra-rotamer energy matrix
			computeEntropyEmatrixMaster(numRes,runName+".log",intraName,sParams,false,false,
					false,intraEnergies,null,false,0.0f,null,false,true); 
		
		sParams.put("PDBNAME", origPDB);
		m = new Molecule();
		setupMolSystem(m,sParams,false,null);
		
		if (ligStrNum>=0) //the ligand is not used here
			m.deleteStrand(ligStrNum);
		
		intraEnergies = (float [][])readObject(intraName,false);
		
		return intraEnergies;
	}
	
	//Reads in the amino acid reference energies (if used);
	private float [] getResEntropyEmatricesEref(Hashtable sParams, boolean useEref){
		
		float eRef[] = new float[resAllowed.length];
		
		if (useEref){ //use AA reference energies		
			for (int i=0; i<numTokens((String)sParams.get("AAREF")); i++){
				String aaRef = getToken((String)sParams.get("AAREF"),i+1);
				float e = new Double(getToken((String)sParams.get("EREF"),i+1)).floatValue();
				int ind = getAllowedAAind(aaRef);
				if (ind>=0)
					eRef[ind] = e;
			}
		}
		else {
			for (int i=0; i<eRef.length; i++)
				eRef[i] = 0.0f;
		}
		
		return eRef;
	}
	
	//Reads in the pairwise energy matrix for the entropy computation
	private float [][][][] readPairMatrixEntropy(String fName, int numRes){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			return null;
		}
		
		float asasE[][][][] = new float[numRes][][][];

		boolean done = false;
		String str = null;
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {				
				int i = new Integer(getToken(str,1)).intValue();
				String name = getToken(str,2);
				asasE[i] = (float [][][])readObject(name,false);
				if (asasE[i]==null){
					System.out.println("ERROR: Could not read data from file "+name);
					System.exit(1);
				}
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
		
		return asasE;
	}
	
	//Determines how many residue positions in the system strand numbered from (pos+1...) (pos is strand-relative numbering)
	//		are within dist from residue position pos
	public boolean [] getProxAS(Molecule m, int pos, float dist, boolean as[]){
		
		int residueMap[] = new int[2];
		residueMap[0] = pos;
		
		for (int i=pos+1; i<m.strand[sysStrNum].numberOfResidues; i++){
			
			if (i!=pos){
				
				boolean done = false;
				
				residueMap[1] = i;
				
				Molecule m1 = getASASMolEntropy(m, residueMap);
				
				StrandRotamers sysLR = new StrandRotamers(rl,m1.strand[0]);
				
				for (int j=0; j<m1.numberOfResidues; j++){
					for (int r=0; r<resAllowed.length; r++){
						sysLR.setAllowable(j,resAllowed[r]);
						m1.residue[j].flexible = true;
					}
				}
				
				for(int q1=0;q1<sysLR.getNumAllowable(0);q1++) {
					
					int AAindex1 = sysLR.getIndexOfNthAllowable(0,q1);
					sysLR.changeResidueType(m1,0,rl.getAAName(AAindex1),true,true);
					
					for(int q2=0;q2<sysLR.getNumAllowable(1);q2++) {
						
						int AAindex2 = sysLR.getIndexOfNthAllowable(1,q2);
						sysLR.changeResidueType(m1,1,rl.getAAName(AAindex2),true,true);
						
						int numRot1 = rl.getNumRotForAAtype(AAindex1);
						
						int w1 = 0;
						if (numRot1<=0)
							w1 = -1;
						
						while ((w1<numRot1)&&(!done)){
							
							if (w1!=-1)
								sysLR.applyRotamer(m1, 0, w1);
							
							int numRot2 = rl.getNumRotForAAtype(AAindex2);
							
							int w2 = 0;
							if (numRot2<=0)
								w2= -1;
							
							while ((w2<numRot2)&&(!done)){
								
								if (w2!=-1)
									sysLR.applyRotamer(m1, 1, w2);
								
								Residue r1 = m1.strand[0].residue[0];
								Residue r2 = m1.strand[0].residue[1];
								
								if (r1.getDist(r2,false)<=dist){
									as[i] = true;
									done = true;
								}
								
								w2++;
							}
							
							w1++;
						}
						if (done)
							break;
					}
					if (done)
						break;
				}
				if (!done)
					as[i] = false;
			}
		}
		
		return as;
	}
	
	//Returns the backbone-dependent AA probabilities for each residue position in the protein
	private double [][] getPDBstats(String fName, Molecule m, int numRes){		
		
		double aaProbAll[][][] = readAAprobFromPDB(fName); //read the backbone-dependent AA probabilities from file fName
			
		Backbone bb = new Backbone();
		
		int bbBins[][] = new int[numRes][2];
		
		//Get the phi/psi bins to which each of the protein residue positions is mapped
		for (int i=0; i<numRes; i++){
			double phi = bb.getFiPsi(m, sysStrNum, i, 0);
			double psi = bb.getFiPsi(m, sysStrNum, i, 1);
			
			if (phi<0.0)
				phi = 360.0 + phi;
			if (psi< 0.0)
				psi = 360.0 + psi;
			
			int phiBin = (int)(phi/10.0);
			int psiBin = (int)(psi/10.0);
			
			bbBins[i][0] = phiBin;
			bbBins[i][1] = psiBin;
		}
		
		//Get the AA probabilities from PDB for each residue position in the protein
		double aaProbForSystem[][] = new double[numRes][resAllowed.length];
		for (int i=0; i<numRes; i++){
			for (int j=0; j<resAllowed.length; j++){
				aaProbForSystem[i][j] = aaProbAll[bbBins[i][0]][bbBins[i][1]][j];
			}
		}
		
		return aaProbForSystem;
	}
	
	//Reads in the backbone-dependent AA probabilities from pdbStatsFile
	private double [][][] readAAprobFromPDB(String fName){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... backbone-dependent amino acid probability file not found");
			System.exit(1);
		}
		
		double aaProbAll[][][] = null;

		boolean done = false;
		String str = null;
		boolean firstLine = true;
		int numBins = 0;
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {
				if (firstLine){
					numBins = new Integer(getToken(str,1)).intValue();
					aaProbAll = new double[numBins][numBins][resAllowed.length];
					firstLine = false;
				}
				else {					
					int phiBin = new Integer(getToken(str,1)).intValue();
					int psiBin = new Integer(getToken(str,2)).intValue();
					String name = getToken(str,3);
					double prob = new Double(getToken(str,4)).doubleValue();
					
					int aaInd = getAllowedAAind(name);
					if (aaInd>=0)					
						aaProbAll[phiBin][psiBin][aaInd] = prob;
				}
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
		
		return aaProbAll;
	}
	
	//Returns the index into resAllowed[] of the amino acid name; returns -1 if name is not found in resAllowed[]
	private int getAllowedAAind(String name){
		
		if ( name.equalsIgnoreCase("HIS") || name.equalsIgnoreCase("HID") || name.equalsIgnoreCase("HIE") )
			name = "HIP";
		
		for (int i=0; i<resAllowed.length; i++){
			if (resAllowed[i].equalsIgnoreCase(name))
				return i;
		}
		return -1;
	}
	
	//Returns the AA index into rotamerIndexOffset to which rotNum belongs
	private int getAAindFromRotNum(int rotNum){
		for (int i=0; i<rotamerIndexOffset.length-1; i++){
			if ( (rotNum>=rotamerIndexOffset[i]) && (rotNum<rotamerIndexOffset[i+1]) )
				return i;
		}
		if (!(rotNum>=totalNumRotamers))
			return (rotamerIndexOffset.length-1);
		else
			return -1;
	}
	
	//Distributes the different types of energy computation for the entropy calculation
	private void computeEntropyEmatrixMaster(int numRes, String runName, String matrixName, Hashtable sParams, boolean doMinimize, boolean minimizeBB, boolean shellRun, 
			float bbEnergies[][], float asasE[][][][], boolean compASASdist, float dist, boolean asDist[][], boolean compEref, boolean intraRun){
		
		
		int mutEnerMatrixSize = 0;
		int residueMap[] = null;
		OneMutation mutArray[] = null;
		int numInAS = numRes;
		
		int numMut = 0;
		
		if (compASASdist){ //compute the min distance between any pair of rotamers for each pair of residue positions
			mutArray = new OneMutation[numInAS];
			for (int i=0; i<asDist.length; i++)			
				mutArray[i] = new OneMutation();
			
			numMut = numRes;
		}
		else {
			if (shellRun) {//computing rotamer-to-backbone energies	
				
				System.out.println("Starting rot-to-backbone energy computation..");
				
				mutEnerMatrixSize = 1 + totalNumRotamers*numInAS;
				
				bbEnergies = new float[1][mutEnerMatrixSize];
				for(int i=0; i<1; i++) {
					for(int j=0; j<mutEnerMatrixSize; j++){
						bbEnergies[i][j] = 0.0f;
					}
				}
				
				if (!compEref)
					residueMap = new int[numRes];
				else {
					numInAS = 1;
					residueMap = new int[numInAS];
				}
				
				mutArray = new OneMutation[numRes];
				for (int i=0; i<mutArray.length; i++){
					mutArray[i] = new OneMutation();
					mutArray[i].flagMutType = "SHL-AS";
				}
				
				numMut = numRes;
			}
			
			else if (intraRun) { //computing intra energies
				
				System.out.println("Starting intra-rot energy computation..");
				
				mutEnerMatrixSize = 1 + totalNumRotamers*numInAS;
				
				bbEnergies = new float[mutEnerMatrixSize][1];
				for(int i=0; i<mutEnerMatrixSize; i++) {
					for(int j=0; j<1; j++){
						bbEnergies[i][j] = 0.0f;
					}
				}
				
				numInAS = 1;
				residueMap = new int[numInAS];
				
				mutArray = new OneMutation[numRes];
				for (int i=0; i<mutArray.length; i++){
					mutArray[i] = new OneMutation();
					mutArray[i].flagMutType = "INTRA";
				}
				
				numMut = numRes;
			}
			
			else { //AS-AS energies
				
				System.out.println("Starting rot-to-rot energy computation..");
				
				numInAS = 2;
				
				int numPairs = 0;
				for (int i=0; i<asasE.length; i++){
					for (int j=i+1; j<asasE[0].length; j++){
						if (asasE[i][j]!=null)
							numPairs++;
					}
				}
				mutArray = new OneMutation[numPairs];
				
				int curPair = 0;
				for (int i=0; i<asasE.length; i++){
					for (int j=i+1; j<asasE[0].length; j++){
						if (asasE[i][j]!=null){				
							mutArray[curPair] = new OneMutation();
							mutArray[curPair].flagMutType = "AS-AS";
							mutArray[curPair].resMut = new int[2];
							mutArray[curPair].resMut[0] = i; //strand-relative numbering (system strand)
							mutArray[curPair].resMut[1] = j;
							curPair++;
							
							asasE[i][j] = new float[totalNumRotamers][totalNumRotamers];
						}
					}
				}
				
				numMut = numPairs;
			}
		}
		
		
		MutationManager mutMan = new MutationManager(runName,mutArray,true);
		mutMan.setResidueMap(residueMap);
		mutMan.setRotamerIndexOffset(rotamerIndexOffset);
		mutMan.setLigType(null);
		mutMan.setarpFilenameMin(matrixName);
		mutMan.setPairEMatrixMin(bbEnergies);
		mutMan.setParams(sParams);
		mutMan.setStericThresh(stericThresh);
		mutMan.setNumInAS(numInAS);
		mutMan.numTotalRotamers(totalNumRotamers);
		mutMan.numResAllowed(numAAallowed);
		mutMan.setIsLigAA(true);
		mutMan.setComputeEVEnergy(true);
		mutMan.setDoMinimization(doMinimize);
		mutMan.setMinimizeBB(minimizeBB);
		mutMan.setCalculateVolumes(false);
		mutMan.setLigPresent(false);
		mutMan.setDistDepDielect(distDepDielect);
		mutMan.setDielectConst(dielectConst);
		mutMan.setDoDihedE(doDihedE);
		mutMan.setDoSolvationE(doSolvationE);
		mutMan.setSolvScale(solvScale);
		mutMan.setVdwMult(softvdwMultiplier);
		mutMan.setEntropyComp(true);
		mutMan.setPairEntropyMatrix(asasE);
		mutMan.setASdistMatrix(asDist);
		mutMan.setASdist(dist);
		mutMan.setCompASdist(compASASdist);
		mutMan.setCompEref(compEref);

		
		try{
			handleDoMPIMaster(mutMan,numMut);
		}
		catch (Exception e){
			System.out.println("ERROR: "+e);
			System.exit(1);
		}
		
		if (compASASdist){
			asDist = mutMan.getASdistMatrix();
			outputObject(asDist,matrixName);
		}
		else {
			if (shellRun) {			
				bbEnergies = mutMan.getMinEmatrix();
				
				int numCompEntries = 0;				
				for(int i2=0; i2<mutEnerMatrixSize; i2++){
					if ((bbEnergies[0][i2]!=0.0f))
						numCompEntries++;
				}
				System.out.println("Num computed entries: "+numCompEntries);
				
				outputObject(bbEnergies,matrixName);
			}
			else if (intraRun){
				bbEnergies = mutMan.getMinEmatrix();
				
				int numCompEntries = 0;				
				for(int i1=0; i1<mutEnerMatrixSize; i1++){
					if ((bbEnergies[i1][0]!=0.0f))
						numCompEntries++;
				}
				System.out.println("Num computed entries: "+numCompEntries);
				
				outputObject(bbEnergies,matrixName);
			}
			else {
				asasE = mutMan.getPairEntropyEmatrix();		
				
				PrintStream logPS = setupOutputFile(matrixName);
				for (int i=0; i<asasE.length; i++){
					if (asasE[i]!=null){
						String fn = ("peme/pem_entr_"+i);
						logPS.println(i+" "+fn);
						outputObject(asasE[i],fn);
						logPS.flush();
					}
				}
				logPS.close();
			}
		}
	}
	
	private CommucObj handleDoResEntropySlave(CommucObj cObj){
		
		long startTime = System.currentTimeMillis();
		
		boolean ligPresent = cObj.ligPresent;

		//Setup the molecule system
		Molecule m = new Molecule();
		setupMolSystem(m,cObj.params,false,null);
		
		
		if (cObj.compASdist){ //AS-AS distance computation
			cObj.asDist = getProxAS(m,cObj.mutationNumber,cObj.dist,cObj.asDist);
			cObj.compEE = new SamplingEEntries[0];
		}
		else { //AS-AS, INTRA, or SHL-AS energy computation
			
			int numInAS = cObj.numInAS;
		
			float mutationEnergiesMin[][] = null;
			int resMut[] = new int[numInAS];
			int residueMap[] = new int[cObj.residueMap.length];
			
			boolean shellRun = false; ligPresent = false; boolean intraRun = false; boolean templateOnly = false;
			
			int mutEnerMatrixSize = -1;
			
			if (cObj.flagMutType.compareTo("AS-AS")==0){ //AS-AS run
				
				mutEnerMatrixSize = 1 + cObj.numTotalRotamers*numInAS;				
				mutationEnergiesMin = new float[mutEnerMatrixSize][mutEnerMatrixSize * 2];
				for(int i=0; i<mutEnerMatrixSize; i++) {
					for(int j=0; j<2*mutEnerMatrixSize; j++){
						mutationEnergiesMin[i][j] = 0.0f;
					}
				}
				
				for (int i=0; i<numInAS; i++)
					resMut[i] = 1;
				
				m = getASASMolEntropy(m,cObj.residueMap);
				
				residueMap[0] = 0;
				residueMap[1] = 1;
			}
			
			else if (cObj.flagMutType.compareTo("INTRA")==0){
				
				intraRun = true;
				
				mutEnerMatrixSize = 1 + cObj.numTotalRotamers*numInAS;
				mutationEnergiesMin = new float[mutEnerMatrixSize][1];
				for(int i=0; i<mutEnerMatrixSize; i++) {
					for(int j=0; j<1; j++){
						mutationEnergiesMin[i][j] = 0.0f;
					}
				}
				
				m = getASASMolEntropy(m,cObj.residueMap);
				
				residueMap = new int[1];
				residueMap[0] = 0;
				resMut = new int[1];
				resMut[0] = 1;
			}
			
			else { //SHL-AS run
				
				shellRun = true; //we are computing the protein residue-to-template energies
							
				int numProx = 0;
				boolean asProx[] = new boolean[m.strand[sysStrNum].numberOfResidues];
				Residue r1 = m.strand[sysStrNum].residue[cObj.mutationNumber];
				
				for (int i=0; i<asProx.length; i++){
					
					if (i!=cObj.mutationNumber){
						
						Residue r2 = m.strand[sysStrNum].residue[i];
						
						if (r1.getDist(r2,true)<=cObj.dist){
							asProx[i] = true;
							numProx++;
						}
						else
							asProx[i] = false;
					}
					else {
						asProx[i] = true;
						numProx++;
					}
				}
				
				residueMap = new int[numProx];
				int curInd = 0;
				int curRes = -1;
				for (int i=0; i<asProx.length; i++){						
					if (asProx[i]){							
						residueMap[curInd] = i;
						if (i==cObj.mutationNumber)
							curRes = curInd;
						curInd++;
					}
				}
				
				m = getASASMolEntropy(m,residueMap);
				
				if (!cObj.compEref){
					residueMap = new int[m.strand[0].numberOfResidues];
					resMut = new int[residueMap.length];
					for (int i=0; i<residueMap.length; i++){
						residueMap[i] = i;
						resMut[i] = 0;
					}
					resMut[curRes] = 1;
					numInAS = residueMap.length;
				}
				else { //SHL-AS for the reference energy computation
					residueMap = new int[1];
					residueMap[0] = curRes;
					resMut = new int[1];
					resMut[0] = 1;
				}
				
				mutEnerMatrixSize = 1 + cObj.numTotalRotamers*numInAS;
				mutationEnergiesMin = new float[1][mutEnerMatrixSize * 2];
				for(int i=0; i<1; i++) {
					for(int j=0; j<2*mutEnerMatrixSize; j++){
						mutationEnergiesMin[i][j] = 0.0f;
					}
				}
			}
			
			RotamerSearch rs = null;
			rs = new RotamerSearch(m, sysStrNum, hElect, hVDW, hSteric, true,
						true, cObj.epsilon, cObj.stericThresh, cObj.distDepDielect, cObj.dielectConst, cObj.doDihedE,cObj.doSolvationE,cObj.solvScale,cObj.vdwMult, rl);
			
			for(int j=0; j<numInAS; j++) {
				for(int q=0;q<resAllowed.length;q++)
					rs.setAllowable(residueMap[j],resAllowed[q]);
			}
			
			mutationEnergiesMin = rs.simplePairwiseMutationAllRotamerSearch(residueMap,
					numInAS,resAllowed,resAllowed.length,cObj.rotamerIndexOffset,
					cObj.numTotalRotamers,cObj.doMinimization,ligPresent,shellRun,intraRun,
					resMut,mutationEnergiesMin,mutEnerMatrixSize,cObj.minimizeBB,templateOnly);
			
			
			if ((cObj.flagMutType.compareTo("SHL-AS")==0)&&(!cObj.compEref)){
				mutEnerMatrixSize = 1 + cObj.numTotalRotamers;
				float tmpM[][] = new float[1][mutEnerMatrixSize];
				int res = -1;
				for (int i=0; i<resMut.length; i++){
					if (resMut[i]!=0){
						res = i;
						break;
					}
				}
				for (int i=0; i<tmpM[0].length-1; i++){
					tmpM[0][i+1] = mutationEnergiesMin[0][1+res*cObj.numTotalRotamers+i];
				}
				mutationEnergiesMin = tmpM;
			}
			
			
			long stopTime = System.currentTimeMillis();
			cObj.elapsedTime = Math.round((stopTime - startTime) / 1000.0f);
			
			int matrixSize2 = 1;
			if (mutationEnergiesMin[0].length>1)
				matrixSize2 = mutEnerMatrixSize;
			
			//Store the information in less space to allow the master node to buffer several cObj at once
			int numCompEntries = 0;
			for(int i=0; i<mutationEnergiesMin.length; i++) {
				for(int j=0; j<matrixSize2; j++) {
					if (mutationEnergiesMin[i][j]!=0.0f)
						numCompEntries++;
				}
			}
			
			cObj.compEE = new SamplingEEntries[numCompEntries];//computed new entries
			for (int i=0; i<numCompEntries; i++)
				cObj.compEE[i] = new SamplingEEntries();
			
			int curE = 0;
			for(int i=0; i<mutationEnergiesMin.length; i++) {
				for(int j=0; j<matrixSize2; j++) {
					if (mutationEnergiesMin[i][j]!=0.0f){
						cObj.compEE[curE].index1 = i;
						cObj.compEE[curE].index2 = j;
						cObj.compEE[curE].minE = mutationEnergiesMin[i][j];
						curE++;
					}
				}
			}
			
			mutationEnergiesMin = null;
		}
		
		return cObj;
	}
	
	//Returns a molecule that contains only the residues in the system strand (sysStrNum) of molecule m that are specified by residueMap[];
	//	This function is used for the pairwise energy matrix computation in the residue entropy calculations
	private Molecule getASASMolEntropy (Molecule m, int residueMap[]){
		
		Molecule m1 = new Molecule();
		
		for (int i=0; i<residueMap.length; i++){
			
			Residue oldResidue = m.strand[sysStrNum].residue[residueMap[i]];
			
			Residue newResidue = new Residue();
			newResidue.name = oldResidue.name;
			newResidue.fullName = oldResidue.fullName;
			
			for (int j=0; j<oldResidue.numberOfAtoms; j++){
				
				Atom oldAtom = oldResidue.atom[j];
				
				Atom newAtom = new Atom(oldAtom.name,oldAtom.coord[0],oldAtom.coord[1],oldAtom.coord[2]);
				newAtom.modelAtomNumber = oldAtom.modelAtomNumber;
				newAtom.strandNumber = oldAtom.strandNumber;
				newAtom.elementType = oldAtom.elementType;
				newResidue.addAtom(newAtom);
			}
			
			m1.addResidue(0,newResidue);
		}
		
		//Determine the bonds between the atoms in the molecule
		m1.determineBonds();
		
		// Assign the molecule relative atom numbers
		m1.updateMoleculeAtomNumbers();
		
		m1.strand[0].isProtein = true;
		
		return m1;
	}
	
	//Compute the reference energies for all amino acids
	public void handleCompEref (String s, Hashtable curScope, Hashtable sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Mutation search parameter filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new Hashtable(hashSize);
	
			// Read System parameters
			readConfigFile(sParams,getToken(s,2));
			
			// Read Mutation search parameters
			readConfigFile(sParams,getToken(s,3));
		}
		
		String runName = (String)sParams.get("RUNNAME");
		
		String protPDBname = (String)sParams.get("PROTPDBNAME");
		int numPDBfiles = (new Integer((String)sParams.get("NUMPDBFILES"))).intValue();
		boolean usePDBstats = (new Boolean((String)sParams.get("USEPDBSTATS"))).booleanValue();
		String pdbStatsFile = (String)sParams.get("PDBSTATSFILE");
		float h = (new Float((String)sParams.get("H"))).floatValue();
		float a = (new Float((String)sParams.get("A"))).floatValue();
		float dist = (new Float((String)sParams.get("DIST"))).floatValue();
		int maxIt = (new Integer((String)sParams.get("MAXIT"))).intValue();
		float eps = (new Float((String)sParams.get("EPS"))).floatValue();
		String eRefName = (String)sParams.get("INITEREFNAME");
		
		
		String pdbFiles[] = getErefPDBfiles(protPDBname,numPDBfiles); //get the PDB filenames
		numPDBfiles = pdbFiles.length;
		
		float aaResE[][][] = new float[numPDBfiles][][];
		float intraEnergies[][][] = new float[numPDBfiles][][];
		int defAA[][] = new int[numPDBfiles][];
		int numRes[] = new int[numPDBfiles];
		double aaProbPDB[][][] = new double[numPDBfiles][][];
		
		for (int i=0; i<numPDBfiles; i++){ //read in or compute the energy matrices for each of the PDB files
			
			sParams.put("PDBNAME",pdbFiles[i]);
		
			//Setup the molecule system
			Molecule m = new Molecule();
			setupMolSystem(m,sParams,false,null);
			
			numRes[i] = m.strand[sysStrNum].numberOfResidues;
			
			defAA[i] = new int[numRes[i]];
			for (int j=0; j<defAA[i].length; j++)
				defAA[i][j] = getAllowedAAind(m.strand[sysStrNum].residue[j].name);			
			
			if (usePDBstats)
				aaProbPDB[i] = getPDBstats(pdbStatsFile,m,numRes[i]);
			
			intraEnergies[i] = getResEntropyEmatricesIntra(pdbFiles[i], sParams, numRes[i], m, runName);
				
			String aaResEname = pdbFiles[i]+"_aaRefE.dat";
			aaResE[i] = (float [][])readObject(aaResEname,false); //check if already computed
			if (aaResE[i]==null){
				computeEntropyEmatrixMaster(numRes[i],runName+".log",aaResEname,sParams,false,false,true,aaResE[i],null,false,dist,null,true,false);
				aaResE[i] = (float [][])readObject(aaResEname,false);
			}
		}
		
		
		float eRef[] = getInitEref(aaResE, intraEnergies, usePDBstats, aaProbPDB, numRes, defAA, eRefName, maxIt, h, a, eps, pdbFiles);
		
		
		//Output the computed amino acid reference energies
		PrintStream logPS = setupOutputFile(runName);
		
		for (int i=0; i<eRef.length; i++)
			logPS.print(resAllowed[i]+" ");
		logPS.println();
		
		for (int i=0; i<eRef.length; i++)
			logPS.print(eRef[i]+" ");
		logPS.println();
		
		logPS.flush();
		logPS.close();
	}
	
	//Get the initial amino acid reference energies
	private float [] getInitEref(float mutationEnergiesMin[][][], float intraEnergies[][][], boolean usePDBstats, 
			double aaProbPDB[][][], int numRes[], int defAA[][], String eRefName, int maxIt,
			float h, float a, float eps, String pdbFiles[]){
		
		final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)		
		
		float eRef[] = readErefFile(eRefName);		
		
		if (eRef!=null){ //already computed
			return eRef;
		}
		else { //not computed, so compute now
			
			eRef = new float[resAllowed.length];
			float bestEref[] = new float[resAllowed.length];
			
			for (int i=0; i<eRef.length; i++){
				eRef[i] = 0.0f;
				bestEref[i] = 0.0f;
			}		
			
			bestEref = getInitErefHelper(mutationEnergiesMin, intraEnergies, usePDBstats, aaProbPDB,
					numRes, eRef, bestEref, defAA, maxIt, h, a, eps, constRT);
			
			
			PrintStream logPS = setupOutputFile(eRefName);
			logPS.println("% AA Eref");
			
			for (int i=0; i<bestEref.length; i++)
				logPS.println(resAllowed[i]+" "+bestEref[i]);
			
			for (int i=0; i<pdbFiles.length; i++){
				
				logPS.println(pdbFiles[i]);
			
				float defAAprob[] = new float[numRes[i]];
				for (int j=0; j<numRes[i]; j++){				
					double aaProb[] = getAAprob(mutationEnergiesMin[i], intraEnergies[i], j, bestEref, constRT, usePDBstats, aaProbPDB[i]);
					if (defAA[i][j]>=0)
						defAAprob[j] = (float)aaProb[defAA[i][j]];
				}			
				
				for (int j=0; j<numRes[i]; j++)
					logPS.print(defAAprob[j]+" ");
				logPS.println();		
				
				logPS.flush();
			}
			logPS.close();
			
			return bestEref;
		}
	}
	
	//Performs steepest-descent-based minimization for the amino acid reference energies
	private float [] getInitErefHelper(float mutationEnergiesMin[][][], float intraEnergies[][][], boolean usePDBstats, 
			double aaProbPDB[][][], int numRes[], float eRef[], float bestEref[],
			int defAA[][], int maxIt, double h, float a, float eps, double constRT){
		
		double delta_h = h/maxIt;
		
		for (int k=0; k<maxIt; k++){
			
			System.out.print(".");
			
			for (int prot=0; prot<numRes.length; prot++){
				
				System.out.print("*");
		
				for (int res=0; res<numRes[prot]; res++){
					
					if ( defAA[prot][res]>=0 ){ //not Pro (Gly should be computed here)
						
						for (int curAA=0; curAA<eRef.length; curAA++){
							
							double aaProbCur[] = getAAprob(mutationEnergiesMin[prot], intraEnergies[prot], res, eRef, constRT, usePDBstats, aaProbPDB[prot]);
							
							eRef[curAA] += h;
							
							double aaProbNew[] = getAAprob(mutationEnergiesMin[prot], intraEnergies[prot], res, eRef, constRT, usePDBstats, aaProbPDB[prot]);
							
							eRef[curAA] -= h;
							
							float curEref = eRef[curAA];
							
							eRef[curAA] = curEref + a * ( (float)((aaProbNew[defAA[prot][res]]-aaProbCur[defAA[prot][res]]) / h) );
							
							/*double df = (aaProbNew[defAA[prot][res]]-aaProbCur[defAA[prot][res]]);
							if (aaProbCur[defAA[prot][res]]!=0.0)
								df /= aaProbCur[defAA[prot][res]];
							eRef[curAA] = curEref + a * ( (float)(df / h) );*/
						}
					}
				}
			}
			
			//h -= delta_h;
		}
		
		return eRef;
	}
	
	//Get the probability for all amino acid types at residue poistion i
	private double [] getAAprob(float mutationEnergiesMin[][], float intraEnergies[][], int i, float eRef[], double constRT, 
			boolean usePDBstats, double aaProbPDB[][]){
		
		double aaProbBBE[] = null;
		
		double normFactor = 0.0;
		for (int j=0; j<totalNumRotamers; j++){
			float E = mutationEnergiesMin[0][1+i*totalNumRotamers+j] + intraEnergies[1+i*totalNumRotamers+j][0] - eRef[getAAindFromRotNum(j)];
			normFactor += Math.exp( -E / constRT);
		}
		
		double rotProb[] = new double[totalNumRotamers];
		for (int j=0; j<totalNumRotamers; j++){
			float E = mutationEnergiesMin[0][1+i*totalNumRotamers+j] + intraEnergies[1+i*totalNumRotamers+j][0] - eRef[getAAindFromRotNum(j)];
			rotProb[j] = Math.exp( -E / constRT) / normFactor;
		}
		
		aaProbBBE = new double[resAllowed.length];
		int curInd = 0;
		for (int j=0; j<aaProbBBE.length; j++){
			
			int numCurRot = rl.getNumRotamers(resAllowed[j]);
			if (numCurRot==0)
				numCurRot = 1;
			
			aaProbBBE[j] = 0.0;
			for (int r=0; r<numCurRot; r++){
				aaProbBBE[j] += rotProb[curInd];
				curInd++;
			}
		}
		
		
		//Compute the unnormalized AA probabilities as a weighted function of backbone energies (if included) and PDB stats (if included)
		double aaProbUnNorm[] = new double[resAllowed.length];
		double aaNorm = 0.0;
		for (int j=0; j<resAllowed.length; j++){
			
			aaProbUnNorm[j] = 1.0;
			
			aaProbUnNorm[j] *= aaProbBBE[j];
			if (usePDBstats)
				aaProbUnNorm[j] *= aaProbPDB[i][j];
			
			aaNorm += aaProbUnNorm[j];
		}
		
		//Normalize the probabilities
		double aaProb[] = new double[resAllowed.length];
		for (int j=0; j<resAllowed.length; j++){
			if (aaNorm!=0.0)
				aaProb[j] = aaProbUnNorm[j]/aaNorm;
			else
				aaProb[j] = 0.0;
		}
		
		return aaProb;
	}
	
	//Reads the pdb filenames for the amino acid reference energy computation
	private String [] getErefPDBfiles(String fName, int numFiles){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			System.out.println(" ... pdbs config file not found");
			System.exit(1);
		}
		
		String pdbFiles[] = new String[numFiles];

		boolean done = false;
		String str = null;
		int curFile = 0;
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {
				if (curFile>=numFiles)
					break;
				else {
					pdbFiles[curFile] = getToken(str,1);
					curFile++;
				}
			}
		}
		
		if (curFile<numFiles){
			String tmp[] = new String[curFile];
			System.arraycopy(pdbFiles, 0, tmp, 0, tmp.length);
			pdbFiles = tmp;
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
		
		return pdbFiles;
	}
	
	//Reads the amino acid rederence energies from file fName
	private float [] readErefFile(String fName){
		
		BufferedReader bufread = null;
		try {
			File file = new File(fName);
			FileReader fr = new FileReader(file);
			bufread = new BufferedReader(fr);
		}
		catch (FileNotFoundException e) {
			return null;
		}
		
		float eRef[] = new float[resAllowed.length];

		boolean done = false;
		String str = null;
		
		while (!done) {
			try {
				str = bufread.readLine();
			}
			catch ( Exception e ){
			 System.out.println("ERROR: An error occurred while reading input");
			 System.exit(0);
			}

			if (str == null) // stop if we've reached EOF
				done = true;
			else if (getToken(str,1).charAt(0)=='%') //skip comment lines
				continue;
			else {								
				String name = getToken(str,1);
				float E = new Float(getToken(str,2)).floatValue();
				
				int aaInd = getAllowedAAind(name);
				if (aaInd>=0)					
					eRef[aaInd] = E;
			}
		}
					
		// We're done reading them in
		try {
			bufread.close();
		}
		catch(Exception e)
		{}
		
		return eRef;
	}
	
	private void selectResidues(String s, Hashtable curScope, Hashtable sParams) {

		// Takes the following parameters
		// 1: System parameter filename (string)
		// 2: Residue search filename (string)		
		
		//Only read system and mutation files if sParams is null
		
		if (sParams==null){ //parameter files not read yet, so read them in
			sParams = new Hashtable(hashSize);
	
			// Read System parameters
			readConfigFile(sParams,getToken(s,2));
			
			// Read Mutation search parameters
			readConfigFile(sParams,getToken(s,3));
		}
		
		String runName = (String)sParams.get("RUNNAME");
		int numRes = (new Integer((String)sParams.get("NUMRES"))).intValue();
		int residues[] = new int[numRes];
		float dist[] = new float[numRes];
		boolean ligPresent = (new Boolean((String)sParams.get("LIGPRESENT"))).booleanValue();
		String ligType = null;
		if (ligPresent)
			ligType = (String)(sParams.get("LIGTYPE"));
		
		String resString = (String)sParams.get("RESIDUES");
		String distString = ((String)sParams.get("DIST"));
		for (int i=0; i<numRes; i++){
			residues[i] = new Integer((String)getToken(resString,i+1)).intValue();
			dist[i] = new Float((String)getToken(distString,i+1)).floatValue();
		}
		
		
		Molecule m = new Molecule();
		setupMolSystem(m,sParams,ligPresent,ligType);
		
		
		boolean asProx[] = new boolean[m.numberOfResidues];
		for (int i=0; i<asProx.length; i++)
			asProx[i] = false;
		
		for (int res=0; res<numRes; res++){
			Residue r1 = m.residue[residues[res]];
			
			for (int i=0; i<asProx.length; i++){
				
				if (i!=residues[res]){
					
					Residue r2 = m.residue[i];
					
					if (r1.getDist(r2,true)<=dist[res])
						asProx[i] = true;
				}
				else
					asProx[i] = true;
			}
		}
		
		Molecule m1 = new Molecule(); //add all proximate residues; the connectivity/bonds will not be valid
		for (int i=0; i<asProx.length; i++){
			if (asProx[i]){
				m1.addResidue(0, m.residue[i]);
			}
		}
		saveMolecule(m1,runName+".pdb",0.0f);
	}
	
} // end of KSParser class

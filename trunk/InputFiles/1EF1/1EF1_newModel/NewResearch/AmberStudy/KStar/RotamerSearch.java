///////////////////////////////////////////////////////////////////////////////////////////////
// RotamerSearch.java
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
 * This class provides a variety of tools and search algorithms for
 *  doing rotamer based searching over molecular conformations
 *
 * The system consists of one molecule containing two strands, one
 *  strand that the protein and one strand is the ligand
 * The system strand does not have to contain sequential residues
 *  but it must be made of standard amino acids
 * The ligand strand can only be one 'thing'
 *  -if this 'thing' is an AA then the Lovell Rotamer library is used
 *  -if this 'thing' is not an AA then an nonAARotamers object is used
 * 
 * In this class more variables than one might expect are class
 *  members. This is done to reduce runtime at the cost of increased
 *  memory requirements. For example we could declare curEEnergy
 *  whenever we're going to compute electrostatic energies but
 *  we choose to declare and allocate it once.
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


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.math.*;

public class RotamerSearch implements Serializable
{

	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = true;
	final double constRT = 1.9891/1000.0 * 298.15;   // in kCal / kelvin-mole (the T value here should be consistent with T in EEF1)
	
	Molecule m;		// the molecule
	Amber96ext a96ff;	// the forcefield and energy function to use for energy evaluation etc...
	
	SimpleMinimizer simpMin;	// the simple energy minimizer (side-chains)
	BBMinimizer bbMin = null;	//the backbone minimizer
	boolean eliminatedRotAtRes[] = null;	// rotamers pruned by MinDEE
	boolean splitFlags[][] = null;
	boolean prunedIsSteric[] = null;	// pruned rotamers due to unallowed sterics
	double indIntMinDEE[] = null;	//the single-residue interval term in the MinDEE criterion
	double pairIntMinDEE[] = null;	//the pairwise interval term in the MinDEE criterion
	boolean repeatSearch = false;	// determines if the search must be repeated to achieve the desired accuracy
	int curConf[] = null;		// the current conformation returned by A*
	boolean allPruned = false;	// determines if all the rotamers for a given residue have been pruned by MinDEE;
									// sends this information to the master node in the mutation search
	StericCheck stericF = null;		//the steric filter object
	
	final float stericE = (float)Math.pow(10,38);	// the energy stored for an unallowed steric
	double Ec_const = stericE;	// the minimum lower energy bound for a pruned conformation
	
	final int samplesBB = 1000;//number samples for intra-rotamer energy matrix computation with backbone flexibility
	
	boolean distDepDielect = true; //distance-dependent dielectric
	double dielectConst = 1.0; //the dielectric constant
	
	boolean doDihedE = false; 		// if true dihedral energies are computed and used
									//  during energy minimization. Note that if
									//  energy minimization is NOT used then dihedral
									//  energies are not explicitly computed. In
									//  reality the total energy values are the same
									//  because although we're using AMBER dihedral
									//  energy terms we assume that each dihedral of
									//  each rotamer is at the bottom of an energy
									//  well. Thus without minimization the total
									//  dihedral energy is zero.
	boolean doSolvationE = false; //determines if solvation energies should be computed
	
	double solvScale = 1.0; //the solvation energies scaling factor
	
	RotamerLibrary rl = null; //the rotamer library
	StrandRotamers sysLR = null;// the rotamers object for the system strand
	StrandRotamers ligROT = null;	// the rotamers object for the ligand strand
	int sysStrNum = -1;	// the strand number of the system
	int ligStrNum = -1;	// the strand number of the ligand
	boolean haveLigand = false;	// does this system have a ligand
	boolean isTheLigandAA = false;	// is the ligand an AA
	boolean ligandReady = false;	// is the ligand ready, do we
			// know what it is and is the rotamer handler ready
	double overlapThresh = 1.5f;	// the overlap threshold used for checking sterics
	int curAANum[] = null;	// for each residue in the system strand, the
			// index of the current amino acid type; if the residue is not
			// rotamerizable (ie. it's not flexible) then the curAANum entry
			// should be -1
	int curLigNum = -1;	// the index of the current ligand type
	boolean computeEVEnergy = false;
			// do we compute EV energies during a conformation search
	boolean doMinimization = false;
			// do we some EV minimization steps during a conformation search
	boolean hElect = true;
		// should hydrogens be used in electrostatic energy calculations
	boolean hVDW = true;
		// should hydrogens be used in vdw energy calculations
	boolean hSteric = false; //should hydrogens be used in steric checks
	float curVolume = 0.0f;
		//  this is the amount of active site space filled with residues
	double vdwMultiplier = 1.0f;
		// vdw multiplier used in energy evaluation
	boolean addHydrogens = true;
		// during a mutation, should hydrogens be included when
		//  changing residue type
	boolean connectResidues = true;
		// during a mutation, should a new residue be bonded to
		//  the prior and subsequent resiudes if the numbering
		//  is sequential
	int curConfNum = 0;
	int numMinSteps = 35; //140
		// number of minimization steps to perform by simpmin

	float arpMatrix[][] = null;
		// all rotamer pairs lower min energy bound matrix, created with
		//  simplePairwiseMutationAllRotamerSearch and loaded with
		//  loadPairwiseEnergyMatrices()
	float arpMatrixMax[][] = null;
		// all rotamer pairs lower max energy bound matrix, created with
		//  simplePairwiseMutationAllRotamerSearchMax and loaded with
		//  loadPairwiseEnergyMatricesMax()
	float intraEnergyArray[] = null;
		// stores the intra-residue energies. Used when computing the
		//  arpMatrix.
	int ASAANums[] = null;
		// integer array containing the index for each AS residues
		//  that can be used with rotamerIndexOffset and for the arpMatrix
	int curASRotNum[] = null;
		// integer array containing the currently assumed rotamer for
		//  each amino acid in the active site
		// is allocated during a rotamer search.
		// note that it is _not_ the same size as curAANum
	int curLigRotNum = 0;
		// the current rotamer number of the ligand
	float bestEnergyFound = 9999999.0f;
		// the best energy found thus far
	float bestEMin = 9999999.0f;
	float bestEUnMin = 9999999.0f;
		// the best minimized and unminimized energy found thus far
	BigInteger numConfsTotal = new BigInteger("0");
		// the number of total conformations for the current configuration
		// this is created and computed in computeTotalNumConfs()
	BigInteger numConfsLeft = new BigInteger("0");
		// the number of remaining conformations for the current configuration
		//  updated as the search progresses
	BigInteger numConfsBelowLevel[] = null;
		// the number of conformations below the specified level, level 0
		//  refers to the ligand level, if there's no ligand then level 0
		//  and level 1 have the same value
		// this is created and computed in computeTotalNumConfs()
	BigInteger numConfsAboveLevel[] = null; //at level i, num confs from level i+1 to the last level
	BigInteger numConfsPrunedByE = new BigInteger("0");
		// the number of conformations not minimized because their 'best'
		//  energy (as computed from the arpMatrix) was too unfavorable
		//  based on the accuracy threshold below
	BigInteger numConfsPrunedByS = new BigInteger("0");
		// number of conformations pruned due to a steric clash
	BigInteger numConfsPrunedByMinDEE = new BigInteger("0");	//the number of confs pruned by MinDEE
	BigInteger numConfsEvaluated = new BigInteger("0");
		// number of conformations that got all the way down to the energy
		//  evaluation
		// Note that numConfsPrunedByE + numConfsPrunedByS + numConfsEvaluated
		//  should equal the total number of conformations
	float KSepsilon = 0.03f;
		// the accuracy for computing energies for K*
		// a value of 0.03 means the energies computed
		//  will allow for a calculation of K*_approx
		//  that's within 3% of the true K*
	BigDecimal partial_q = new BigDecimal(0.0);
		// the partially computed partition function (updated as we go)
	BigDecimal partial_p = new BigDecimal(0.0);
		// the bound on the partition function of the pruned conformations
	BigDecimal initial_q = new BigDecimal(0.0);
		// used in mutation search as an initial partial_q if we're
		//  bootstrapping the search
	
	// Note that the mutation search functions in this class are relatively
	//  messy as a result of changing the algorithms multiple times. They
	//  could be rewritten to be much tighter and more elegant.
	
	// A dummy constructor, not very useful
	RotamerSearch() {
	}
	
	// sysStrand must be a valid stand number in theMolec
	RotamerSearch(Molecule theMolec, int systemStrand, boolean hE, boolean hV, boolean hS, boolean addH,
			boolean conRes, float eps, float stericThresh, boolean ddDielect, double dielectC,
			boolean doDihedral, boolean doSolv, double solvScFactor, double vdwMult, RotamerLibrary rlP) {
		
		rl = rlP;
		
		hElect = hE;
		hVDW = hV;
		hSteric = hS;
		addHydrogens = addH;
		connectResidues = conRes;
		KSepsilon = eps;
		overlapThresh = stericThresh;
		vdwMultiplier = vdwMult;
		distDepDielect = ddDielect;
		dielectConst = dielectC;
		doDihedE = doDihedral;
		doSolvationE = doSolv;
		solvScale = solvScFactor;
		
		bestEMin = stericE;
		
		m=theMolec;
		a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
		simpMin = new SimpleMinimizer();
		bbMin = new BBMinimizer();
		sysStrNum = systemStrand;
		curAANum = new int[m.strand[sysStrNum].numberOfResidues];
		
		sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
	}

	// the constructor if you also have a ligand
	RotamerSearch(Molecule theMolec, int systemStrand, int ligandStrand, boolean isLigAA,
			boolean hE, boolean hV, boolean hS, boolean addH,
			boolean conRes, float eps, float stericThresh, boolean ddDielect, double dielectC,
			boolean doDihedral, boolean doSolv, double solvScFactor, double vdwMult, RotamerLibrary rlP) {
		
		rl = rlP;
		
		hElect = hE;
		hVDW = hV;
		hSteric = hS;
		addHydrogens = addH;
		connectResidues = conRes;
		KSepsilon = eps;
		overlapThresh = stericThresh;
		vdwMultiplier = vdwMult;
		distDepDielect = ddDielect;
		dielectConst = dielectC;
		doDihedE = doDihedral;
		doSolvationE = doSolv;
		solvScale = solvScFactor;
		
		bestEMin = stericE;
		
		m=theMolec;
		a96ff = new Amber96ext(m, distDepDielect, dielectConst, doSolvationE, solvScale, vdwMultiplier);
		simpMin = new SimpleMinimizer();
		bbMin = new BBMinimizer();
		sysStrNum = systemStrand;
		ligStrNum = ligandStrand;
		m.strand[ligStrNum].residue[0].flexible = true;
		curAANum = new int[m.strand[sysStrNum].numberOfResidues];
		
		sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);

		haveLigand = true;
		if (isLigAA){			
			ligROT = new StrandRotamers(rl,m.strand[ligStrNum]);
			ligROT.setAllowable(0,m.strand[ligStrNum].residue[0].name);
			
			ligandReady = true;
			isTheLigandAA = true;
		}
		else {
			// the ligROT gets instantiated later (see setupGenericRotamer)
			ligROT = null; // temp line - delete when the above is ready
			ligandReady = false; // the GenericRotamer object is not ready yet
			isTheLigandAA = false;	// a generic rotamer is not an amino acid
		}
	}
	
	// This function adds the AA type named name to the list
	//  of allowable types for residue number resNum in
	//  the system strand (resNum is strand based numbering)
	public void setAllowable(int resNum, String name) {
		sysLR.setAllowable(resNum,name);
		m.residue[resNum].flexible = true;
	}
	
	public void setSplitFlags(int size){
		splitFlags = new boolean[size][size];
		for (int i=0; i<splitFlags.length; i++){
			for (int j=0; j<splitFlags.length; j++){
				splitFlags[i][j] = false;
			}
		}
	}
	
	public void setSplitFlags(boolean spFlags[][]){
		for (int i=0; i<spFlags.length; i++){
			for (int j=0; j<spFlags.length; j++){
				splitFlags[i][j] = spFlags[i][j];
			}
		}
	}
	
	//Return the split flags
	//	This version makes sure that the returned matrix is not associated with the RotamerSearch instance variable
	public boolean[][] getSplitFlags(boolean savedSpF[][]){
		savedSpF = new boolean[splitFlags.length][splitFlags.length];
		for (int i=0; i<splitFlags.length; i++){
			for (int j=0; j<splitFlags.length; j++){
				savedSpF[i][j] = splitFlags[i][j];
			}
		}
		return savedSpF;
	}
	
	public boolean[][] getSplitFlags(){
		return splitFlags;
	}
	
	public BigDecimal getBestE(){
		return new BigDecimal(bestEMin);
	}
	
	// This function clears the list of allowable AA types
	//  for residue number resNum in the system strand
	public void clearAllowable(int resNum) {
		sysLR.clearAllowable(resNum);
	}

	// Refreshes the system strand
	public void refreshSystemStrand(){
		sysLR = new StrandRotamers(rl,m.strand[sysStrNum]);
	}

	// This function computes one energy 
	private float calcTotalSnapshotEnergy(){
		
		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,-1); //compute the energy
	
		return (float)energyTerms[0]; //the total energy is in energyTerms[0]
	}
	
	// This function computes the energy of the given AAnum with all other residues for which the energy flags are set
	// NOTE: This assumes that Amber is setup and the bbMin is setup (so backbone flexibility is used)
	private float calcPartEnergy(int AAnum){
		
		bbMin.setupPartialAmber();
		double energyTerms[] = a96ff.calculateTotalEnergy(m.actualCoordinates,AAnum); //compute the energy
	
		return (float)energyTerms[0]; //the total energy is in energyTerms[0]
	}

//// BEGIN CHECK_STERICS CODE SECTION
	// This function is similar to RS_CheckStericsWithH
	// This version checks all residues against the target residue rather
	//  than just checking residues up to the specified one in the strand
	private boolean RS_CheckAllStericsWithH(int strandNum, int resNum) {
	
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;

		for(int i=0;i<res.numberOfAtoms;i++) {
			for(int q=0;q<m.numberOfStrands;q++) {
				resToCheck = m.strand[q].numberOfResidues;
				for(int w=0;w<resToCheck;w++) {
					if(!((q==strandNum) && (w==resNum))) {
						for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
							tmpAtm = m.strand[q].residue[w].atom[t];
							if ((res.atom[i].distance(tmpAtm)) < (((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
								if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
									return false;
								}
							}
						}
					}
				}
			}
		}
	
		// If you got here then everything passed
		return true;
	}

	
	// This function is similar to RS_CheckSterics
	// This version checks all residues against the target residue rather
	//  than just checking residues up to the specified one in the strand
	private boolean RS_CheckAllSterics(int strandNum, int resNum) {
	
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;
		// It's important to divide below by 100.0 rather than
		//  just 100 as java assumes precision or does some
		//  strange rounding that causes an incorrect results
		//  when 100 is used
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			if (!(res.atom[i].elementType.equalsIgnoreCase("H"))) {
				for(int q=0;q<m.numberOfStrands;q++) {
					resToCheck = m.strand[q].numberOfResidues;
					for(int w=0;w<resToCheck;w++) {
						if(!((q==strandNum) && (w==resNum))) {
							for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
								tmpAtm = m.strand[q].residue[w].atom[t];
								if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
									if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
										if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
											//System.out.println("Clash: " + tmpAtm.moleculeAtomNumber + " with residue " + resNum + " atom " + i);
											//writeCurrentMolecule("what.pdb");
											//System.exit(0);
											return false;
										}
									}
								}
							}
						}
					}
				}
			}		
		}
	
		// If you got here then everything passed
		return true;
	}

//// END CHECK_STERICS CODE SECTION

//// BEGIN HELPER FUNCTION SECTION
	
	// Loads the max pairwise energy matrix
	public void loadPairwiseEnergyMatricesMax(String allRotamerPairsEnergyName) {

		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(allRotamerPairsEnergyName));
			arpMatrixMax = (float [][])in.readObject();
			in.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while reading precomputed max energy matrices");
			System.exit(0);
		}
	}
	
	// Loads the min pairwise energy matrix
	public void loadPairwiseEnergyMatrices(String allRotamerPairsEnergyName) {

		try{
			ObjectInputStream in = new ObjectInputStream(new FileInputStream(allRotamerPairsEnergyName));
			arpMatrix = (float [][])in.readObject();
			in.close();
		}
		catch (Exception e){
			System.out.println(e.toString());
			System.out.println("ERROR: An exception occurred while reading precomputed min energy matrices");
			System.exit(0);
		}		
	}
	
	//Adds the reference energies to the intra-energies (column 0) in arpMatrix[][];
	//If doMinimize is true, then arpMatrixMax is also updated appropriately
	public void addEref(float eRef[], boolean doMinimize, boolean ligPresent, int numInAS){
		
		int ind = 1; //skip the entry [0][0], since this is the fixed template energy
		for (int i=0; i<numInAS; i++){
			for (int j=0; j<eRef.length; j++){
				int numRot = rl.getNumRotForAAtype(j);
				if (numRot==0) //ALA or GLY
					numRot = 1;
				for (int k=0; k<numRot; k++){
					arpMatrix[ind][0] -= eRef[j];
					if (doMinimize)
						arpMatrixMax[ind][0] -= eRef[j];
					ind++;
				}
			}			
		}
		
		if (ligPresent){
			int ligType = ligROT.getIndexOfNthAllowable(0,0);
			int numRot = rl.getNumRotForAAtype(ligType);
			if (numRot==0)
				numRot = 1;
			
			for (int k=0; k<numRot; k++){
				arpMatrix[ind][0] -= eRef[ligType];
				if (doMinimize)
					arpMatrixMax[ind][0] -= eRef[ligType];
				ind++;
			}
		}
	}


	// Computes the best energy (lower bound) using the arpMatrix
	// This energy is rotamer based, that is it computes the best energy
	//  for the current rotamer assignment of each amino-acid	
	private float computeBestRotEnergyBound(int numTotalRotamers, int rotamerIndexOffset[],
		boolean ligPresent) {

		float bestE = 0.0f;
		int indx1 = 0;
		int indx2 = 0;
		
		// Add shell-shell energy
		bestE = arpMatrix[0][0];
	
		if (ligPresent) {
			indx1 = ASAANums.length * numTotalRotamers + curLigRotNum + 1;
			// Ligand shell energy
			bestE += arpMatrix[0][indx1];
			// Ligand intra-residue energy
			bestE += arpMatrix[indx1][0];
			// Ligand pairwise energies
			for(int j=0;j<ASAANums.length;j++) {
				indx2 = j * numTotalRotamers + rotamerIndexOffset[ASAANums[j]] + curASRotNum[j] + 1;
				bestE += arpMatrix[indx1][indx2];
			}
		}
	
		for(int i=0;i<ASAANums.length;i++) {
			indx1 = i * numTotalRotamers + rotamerIndexOffset[ASAANums[i]] + curASRotNum[i] + 1;
			// Add the residue i shell energy
			bestE += arpMatrix[0][indx1];
			// Add the residue i intra-residue energy
			bestE += arpMatrix[indx1][0];
			// Add the pairwise energies
			for(int j=i+1;j<ASAANums.length;j++) {
				indx2 = j * numTotalRotamers + rotamerIndexOffset[ASAANums[j]] + curASRotNum[j] + 1;
				bestE += arpMatrix[indx1][indx2];
			}
		}

		return bestE;
	}

//// END HELPER FUNCTION SECTION

	
//// BEGIN MASTER MUTATION SEARCH SECTION

	// This function is the model for a mutation search
	// This is the function used by the master node to generate a list of
	//  mutations that it wishes to consider.
	// Utilizes a number of helper functions
	public int simpleMasterMutationSearch(int residueMap[], int numInAS,
		String residuesAllowed[], int numResAllowed, int theCurConfNum,
		OneMutation mutArray[], float minVol, float maxVol) {

		curConfNum = theCurConfNum;
		
		masterMutationSearchHelper(0, numInAS, residueMap, residuesAllowed,
			numResAllowed, mutArray, minVol, maxVol);
		
		return curConfNum;
	}

	// This function is similar to mutationSearchHelper
	//  the only difference is that we only compute volumes and an amino acid
	//  level energy approximation. I could have modified that function, but
	//  decided not to so as to keep that function fast (ie. this way the
	//  execution of a bunch of conditionals is saved in the normal search)
	public void masterMutationSearchHelper(int depth, int maxDepth,
		int residueMap[], String residuesAllowed[], int numResAllowed,
		OneMutation mutArray[], float minVol, float maxVol) {
	
		if (depth >= maxDepth) {
			// If we've arrived here then we're ready to
			//  compute a volume and approximate energy
			if(debug){
				System.out.print(".");
			}
			curVolume = sysLR.getCurAllowableVolume();
			if ((curVolume > minVol) && (curVolume < maxVol)) {
				// Add mutation to mutation array
				OneMutation tMut = new OneMutation();
				assignAANums(residueMap,residuesAllowed,numResAllowed);
				tMut.score = new BigDecimal("0.0");  // Added when aap removed 6/23/03
				tMut.resTypes = new String[maxDepth];
				for(int q=0;q<maxDepth;q++) {
					tMut.resTypes[q] = rl.getAAName(curAANum[residueMap[q]]);
				}
				tMut.vol = curVolume;
				if (curConfNum >= mutArray.length) {
					// If there's no space left, make space in mutArray
					OneMutation newArray[] = new OneMutation[mutArray.length + 5000];
					System.arraycopy(mutArray, 0, newArray, 0, mutArray.length);
					mutArray = newArray;
				}
				mutArray[curConfNum] = tMut;
				curConfNum++;
			}
			return;
		}

		// Check with allowed AAs
		for(int q=0;q<sysLR.getNumAllowable(residueMap[depth]);q++) {
			curAANum[residueMap[depth]] = sysLR.getIndexOfNthAllowable(residueMap[depth],q);
			// Change the residue type
			sysLR.changeResidueType(m,residueMap[depth],rl.getAAName(curAANum[residueMap[depth]]),addHydrogens,connectResidues);
			masterMutationSearchHelper(depth+1,maxDepth,residueMap,residuesAllowed,numResAllowed,mutArray,minVol,maxVol);
		}
	}
	
	// Assigns elements of the ASAANums[] array
	private void assignAANums(int residueMap[], String AAList[],
		int numAAAllowed) {
	
		boolean done = false;
		int tmpIndx = 0;
		String tmpStr = null;
		
		ASAANums = new int[residueMap.length];
		
		for(int i=0;i<residueMap.length;i++){
			done = false;
			tmpIndx = 0;
			tmpStr = m.strand[sysStrNum].residue[residueMap[i]].name;
			if (tmpStr.equalsIgnoreCase("hip") ||
				tmpStr.equalsIgnoreCase("hie") ||
				tmpStr.equalsIgnoreCase("hid") ||
				tmpStr.equalsIgnoreCase("his"))
				tmpStr = new String("hip");
			while ((!done) && (tmpIndx < numAAAllowed)) {
				if (tmpStr.equalsIgnoreCase(AAList[tmpIndx])) {
					ASAANums[i] = tmpIndx;
					done = true;
				}
				else
					tmpIndx++;
			}
		}
	}

//// END MASTER MUTATION SEARCH SECTION
///////////////////////////////////////////////////////////////////////////////////	

///////////////////////////////////////////////////////////////////////////////////
//// BEGIN PAIRWISE MUTATION _ALL_ ROTAMER SEARCH SECTION - ENERGY PRECOMPUTATION
	
	// Turns all residues off except the side-chain for the residue specified by
	//  molResNum. Computes the energy of the current system.
	//  It does not do minimization.
	// *Make sure that the forcefield types for the atoms in the
	//  residue of interest have been computed prior to calling
	//  this function as they are not computed in this function
	private float computeEnergyOfOnlyRes(int molResNum) {
		
		boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		float curEnergy = 0.0f;
			
		// Save the energy eval flag, clear them at the same time
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			m.residue[i].setEnergyEval(false, false);
		}

		m.residue[molResNum].setEnergyEval(true, false);
		
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		curEnergy = calcTotalSnapshotEnergy();
			
		// Restore the energy eval and flexibility flags
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
		}
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		
		return curEnergy;
	}
	
	//Turns all template residues on and the side-chains for all flexible residues off;
	//Computes the template energy for the current template; does not do minimization
	// *Make sure that the forcefield types for the atoms in the
	//  residue of interest have been computed prior to calling
	//  this function as they are not computed in this function
	private float computeEnergyOfOnlyTemplate(int residueMap[]) {
		
		boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		float curEnergy = 0.0f;
			
		// Save the energy eval flag, clear them at the same time
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			m.residue[i].setEnergyEval(true, true);
		}

		//Clear the energy evaluation flags for the AS residues and the ligand (if present), so that
		//		the energies only between template residues are computed
		for (int i=0; i<residueMap.length; i++){
			m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(false, true);
		}
		if (haveLigand)
			m.strand[ligStrNum].residue[0].setEnergyEval(false, false);
		
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		curEnergy = calcTotalSnapshotEnergy();
			
		// Restore the energy eval and flexibility flags
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
		}
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		
		return curEnergy;
	}

	// This function helps to compute the min/max pairwise energy matrix, two
	//  residues are allowed to mutate (as specified in residueMutatable)
	//  a steric check is done for each rotamer pair if the steric check
	//  passes then the min/max energy is computed for the pair and is
	//  saved. ALL rotamer pairs that pass the steric threshold are
	//  saved. If a pair doesn't pass the threshold then an energy of
	//  10^38 is assigned.
	// To compute the pairwise interactions of each rotameric position with
	//  the ligand only one residue should be "allowed" in residueMutatable[]
	//  and the ligPresent should be set to true.
	// A shellRun computes the energy of the "allowed" residue with all other
	//  residues that are NOT in residueMap
	// Utilizes a number of helper functions
	public float[][] simplePairwiseMutationAllRotamerSearch(int residueMap[], int numInAS, 
			String residuesAllowed[], int numResAllowed, int rotamerIndexOffset[], int numTotalRotamers,
			boolean searchDoMinimize, boolean ligPresent, boolean shellRun, boolean intraRun,
			int residueMutatable[], float retEMatrix[][], int matrixSizeRows, boolean minimizeBB, boolean templateOnly) {
		
		doMinimization = searchDoMinimize;
		computeEVEnergy = true;
		
		if (ligPresent && !ligandReady){
			System.out.println("Warning: Ligand rotamers not setup, RotamerSearch aborting");
			return null;
		}

		// Prepare Amber
		if(computeEVEnergy){
			// Amber should already be loaded
			// First turn off energy evaluation for all residues
			//   since we're only computing pairwise energies
			//   we only want specific residues on (will be turned on later)
			// If we're doing a shell run then turn the shell on
			if (shellRun) {
				for(int i=0;i<m.numberOfResidues;i++){
					m.residue[i].setEnergyEval(true, true);
					m.residue[i].flexible = false;
				}
				if (ligStrNum != -1) {
					for(int i=0;i<m.strand[ligStrNum].numberOfResidues;i++)
						m.strand[ligStrNum].residue[i].setEnergyEval(false, false);
				}
				for(int i=0;i<numInAS;i++)
					m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(false, true);
			}
			else {
				for(int i=0;i<m.numberOfResidues;i++){
					m.residue[i].setEnergyEval(false, false);
					m.residue[i].flexible = false;
				}
			}
			
			if(ligPresent) {
				a96ff.setLigandNum(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
			}
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null) {
						System.out.println("Error: simpMin not allocated and you are attempting minimization");
						return(null);
					}
					bbMin = null;
				}
				else { //backbone minimization
					if (bbMin == null) {
						System.out.println("Error: bbMin not allocated and you are attempting backbone minimization");
						return(null);
					}
					simpMin = null;
				}
			}
			
			if (shellRun) { //compute the template energies				
				if (!minimizeBB) {//side-chain minimization, so the template is fixed					
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);				
					float minE = calcTotalSnapshotEnergy();
					retEMatrix[0][0] = minE;
					retEMatrix[0][matrixSizeRows+0] = minE;
				}
			}
		}
		
		if ((minimizeBB)&&(templateOnly)){ //the template energies for backbone minimization are computed only once
			for (int i=0; i<residueMap.length; i++){
				String a = rl.getThreeLetAACode(sysLR.getCurRotType(residueMap[i]));
				if ((!a.equalsIgnoreCase("GLY"))&&(!a.equalsIgnoreCase("PRO"))){
					sysLR.changeResidueType(m,residueMap[i],"ala",addHydrogens,connectResidues);
				}
				m.strand[sysStrNum].residue[residueMap[i]].flexible = false;
				m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(false, true);
			}
			a96ff.calculateTypesWithTemplates();
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			bbMin.initialize(m, a96ff, residueMap, sysStrNum); //no ligand for template energy computation					
			pairwiseRotamerEnergyBackboneHelper(-1,-1,-1,-1,residueMap,-1,-1,	shellRun, retEMatrix, matrixSizeRows, true);					
			bbMin = new BBMinimizer(); //reset the backbone minimizer
			return retEMatrix;
		}

		if (intraRun) {
			// Compute the intra-energies for each amino-acid/rotamer
			float intraEnergyArrayMin[] = new float[matrixSizeRows];
			float intraEnergyArrayMax[] = new float[matrixSizeRows];
			for(int i=0;i<intraEnergyArrayMin.length;i++){
				intraEnergyArrayMin[i] = 0.0f;
				intraEnergyArrayMax[i] = 0.0f;
			}
			computeIntraRotEnergies(numInAS,residueMap,residuesAllowed,
				numResAllowed,rotamerIndexOffset,numTotalRotamers,
				intraEnergyArrayMin, intraEnergyArrayMax,ligPresent,residueMutatable,minimizeBB);
			// Start copying at 1 so that [0][0] doesn't get clobbered
			for(int i=1;i<matrixSizeRows;i++){
				retEMatrix[i][0] = intraEnergyArrayMin[i];
				if (retEMatrix[0].length>matrixSizeRows)
					retEMatrix[i][matrixSizeRows+0] = intraEnergyArrayMax[i];
			}
			
			return retEMatrix;
		}
		
		// Initialize curAANum array, we have to do this because we only recurse through the 9 core residues of the active site
		//  and not all 40 residues, so we use a residueMap and we have to do some preinitialization.
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++)
			curAANum[i] = -1;
		
		// Note: In this search we only search over the key active site  residues rather than all residues in the molecule, 
		//		thus maxDepth should be numInAS
		pairwiseEnergyComp (numInAS,residueMap,	residuesAllowed, numResAllowed, rotamerIndexOffset, numTotalRotamers,
				ligPresent, residueMutatable, retEMatrix, matrixSizeRows, shellRun, minimizeBB);

		return retEMatrix;
	}

	
	private void computeIntraRotEnergies(int numInAS, int residueMap[], String residuesAllowed[], int numResAllowed,
		int rotamerIndexOffset[], int numTotalRotamers, float intraEnergyArrayMin[],float intraEnergyArrayMax[],
		boolean ligPresent, int residueMutatable[], boolean minimizeBB) {
	
		//Save the energy eval and flexibilty flags, clear them at the same time
		boolean savedEnergyEvalSC[] = new boolean[m.numberOfResidues];
		boolean savedEnergyEvalBB[] = new boolean[m.numberOfResidues];
		boolean savedFlexible[] = new boolean[m.numberOfResidues];
		for(int i=0;i<m.numberOfResidues;i++){
			savedEnergyEvalSC[i] = m.residue[i].getEnergyEvalSC();
			savedEnergyEvalBB[i] = m.residue[i].getEnergyEvalBB();
			savedFlexible[i] = m.residue[i].flexible;
			m.residue[i].setEnergyEval(false, false);
			m.residue[i].flexible = false;
		}
		
		// Go through each active site residue, each AA type they could be and all their rotamers, 
		// 	saving the computed energies to the appropriate place.
		for(int i=0;i<numInAS;i++) {
			
			if (residueMutatable[i] == 1){ //only compute for those residues that are mutatable
				
				System.out.println();
				
				for(int j=0;j<numResAllowed;j++){
					
					System.out.print(".");
					
					// Apply mutation
					sysLR.changeResidueType(m,residueMap[i],residuesAllowed[j],addHydrogens,connectResidues);
					curAANum[residueMap[i]] = rl.getAARotamerIndex(residuesAllowed[j]);
					m.strand[sysStrNum].residue[residueMap[i]].flexible = true;
					m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(true, false);
	
					// Setup Amber, Setup Minimizer
					if (computeEVEnergy){
						a96ff.calculateTypesWithTemplates();
						a96ff.initializeCalculation();
						a96ff.setNBEval(hElect,hVDW);
						if (doMinimization){
							if (!minimizeBB) //side-chain minimization
								simpMin.initialize(m,sysStrNum,a96ff,sysLR,curAANum,doDihedE,rl);
							else //backbone minimization
								bbMin.initialize(m, a96ff, residueMap, sysStrNum);
						}
					}
	
					
					boolean done = false;
					int curRot = 0;
					int totRotForCur = rl.getNumRotForAAtype(curAANum[residueMap[i]]);
					
					while ( (!done) && ((curRot<totRotForCur) || (totRotForCur==0)) ){
					
						if (totRotForCur!=0) // loop through rotamers
							sysLR.applyRotamer(m, residueMap[i], curRot);
						else // no rotamers, so only the current AA state has to be computed and we exit the while loop
							done = true;
						
						computeIntraRotEnergiesHelper(i, j, curRot, totRotForCur, 
								intraEnergyArrayMin, intraEnergyArrayMax, rotamerIndexOffset, 
								residueMap, numTotalRotamers, numInAS, false, minimizeBB);
						
						curRot++;
					}				
				}
				m.strand[sysStrNum].residue[residueMap[i]].flexible = false;
				m.strand[sysStrNum].residue[residueMap[i]].setEnergyEval(false, false);
			}
		}
		
		if (ligPresent) { //go through all of the ligand rotamers and save the computed intra energies		
			// Apply mutation
			curLigNum = ligROT.getIndexOfNthAllowable(0,0);
			ligROT.changeResidueType(m,0,rl.getAAName(curLigNum),addHydrogens);
			m.strand[ligStrNum].residue[0].flexible = true;
			m.strand[ligStrNum].residue[0].setEnergyEval(true,true);

			// Setup Amber, Setup Minimizer
			if (computeEVEnergy){
				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);
				if (doMinimization){
					if (!minimizeBB) //side-chain minimization
						simpMin.initialize(m,sysStrNum,ligStrNum,a96ff,sysLR,ligROT,curAANum,curLigNum,doDihedE,rl);
					else //backbone minimization
						bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
				}
			}

			boolean done = false;
			int curRot = 0;
			int totRotForCur = rl.getNumRotForAAtype(curLigNum);
			
			while ( (!done) && ((curRot<totRotForCur) || (totRotForCur==0)) ){
				
				if (totRotForCur!=0) // loop through rotamers
					ligROT.applyRotamer(m, 0, curRot);
				else // no rotamers, so only the current AA state has to be computed and we exit the while loop
					done = true;
						
				computeIntraRotEnergiesHelper(-1, -1, curRot, totRotForCur,	intraEnergyArrayMin, 
						intraEnergyArrayMax, rotamerIndexOffset, residueMap, numTotalRotamers, numInAS, true, minimizeBB);
				
				curRot++;
			}
		}
		
		// Restore the energy eval and flexibility flags
		for(int i=0;i<m.numberOfResidues;i++){
			m.residue[i].setEnergyEval(savedEnergyEvalSC[i], savedEnergyEvalBB[i]);
			m.residue[i].flexible = savedFlexible[i];
		}
	}
	
	//Computes the intra-residue min and max energies for a given rotamer
	private void computeIntraRotEnergiesHelper(int curPos, int curAA, int curRot, int totRotForCur, 
			float intraEnergyArrayMin[], float intraEnergyArrayMax[], int rotamerIndexOffset[], 
			int residueMap[], int numTotalRotamers, int numInAS, boolean isLig, boolean minimizeBB){
		
		float curEnergy = 0.0f;	
		float beginE = 0.0f;
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
		
		//Compute min and max energies
		if (doMinimization){								
			
			//Iniitialize Amber for the current residue
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			
			if (!minimizeBB){ //side-chain minimization
				
				//first, compute the initial energy at the new position
				beginE = calcTotalSnapshotEnergy();
								
				//minimize, starting at the initial position
				simpMin.minimize(numMinSteps,isLig);
				curEnergy = calcTotalSnapshotEnergy();
				if (doDihedE) //add dihedral energies
					curEnergy += simpMin.computeDihedEnergy();
				
				//Compare to the min and max energies found so far and update, if necessary
				if (beginE<curEnergy)
					curEnergy = beginE;
				float lE = Math.min(beginE, curEnergy);
				float hE = Math.max(beginE, curEnergy);
				minEnergy = Math.min(minEnergy,lE);									
				maxEnergy = Math.max(maxEnergy,hE);

				m.updateCoordinates();//restore the actualCoordinates array to the initial values
			}
			else { //backbone minimization, so only rotate the backbone O and HN		
				
				int at[] = new int[5];
				int numAtoms = 0;
				Residue r1 = null;
				if (curPos>=0) //AS residue
					r1 = m.residue[residueMap[curPos]];
				else //ligand
					r1 = m.residue[m.strand[ligStrNum].residue[0].moleculeResidueNumber];
				
				numAtoms = r1.numberOfAtoms;
				
				//get the atoms
				for (int i=0; i<numAtoms; i++){
					if (r1.atom[i].name.equalsIgnoreCase("CA"))
						at[0] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("C"))
						at[1] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("O"))
						at[2] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("N"))
						at[3] = r1.atom[i].moleculeAtomNumber;
					else if (r1.atom[i].name.equalsIgnoreCase("H"))
						at[4] = r1.atom[i].moleculeAtomNumber;
				}
				
				//get the CA-C bond
				float dx = m.actualCoordinates[at[1]*3] - m.actualCoordinates[at[0]*3];
				float dy = m.actualCoordinates[at[1]*3+1] - m.actualCoordinates[at[0]*3+1];
				float dz = m.actualCoordinates[at[1]*3+2] - m.actualCoordinates[at[0]*3+2];
				
				//get the N-CA bond
				float dxH = m.actualCoordinates[at[0]*3] - m.actualCoordinates[at[3]*3];
				float dyH = m.actualCoordinates[at[0]*3+1] - m.actualCoordinates[at[3]*3+1];
				float dzH = m.actualCoordinates[at[0]*3+2] - m.actualCoordinates[at[3]*3+2];
				
				//get the center of rotation for O (the actualCoordinates[] of C)
				double center[] = new double[3];
				center[0] = m.actualCoordinates[at[1]*3];
				center[1] = m.actualCoordinates[at[1]*3+1];
				center[2] = m.actualCoordinates[at[1]*3+2];
				
				//get the center of rotation for H (the actualCoordinates[] of N)
				double centerH[] = new double[3];
				centerH[0] = m.actualCoordinates[at[3]*3];
				centerH[1] = m.actualCoordinates[at[3]*3+1];
				centerH[2] = m.actualCoordinates[at[3]*3+2];
				
				float rotForInitPos = bbMin.getMaxDihedRot(); //get the max phi/psi rotation
				
				//Do the sampling and minimization
				for (int curSample=0; curSample<samplesBB; curSample++){
				
					//randomly generate the rotation angle
					float rotChange[] = new float[2];
					Random randNum = new Random();
					
					if (curSample!=0) {
						for (int i=0; i<2; i++)
							rotChange[i] = (randNum.nextFloat()-0.5f)*rotForInitPos*2.0f;
					}
					else {
						for (int i=0; i<2; i++)
							rotChange[i] = 0.0f;
					}
					
					//Compute the energy corresponding to the new positions
					if (curSample!=0){ //do not apply a change for the initial position
						m.rotateAtom(at[2], dx, dy, dz, center[0], center[1], center[2], rotChange[0], false);
						m.rotateAtom(at[4], dxH, dyH, dzH, centerH[0], centerH[1], centerH[2], rotChange[1], false);
					}
					
					//compute the initial energy at the new position
					curEnergy = calcTotalSnapshotEnergy();
					
					//Compare to the min and max energies found so far and update, if necessary;
					//For intra-energies with BB flexibility, the initial point is taken as the max energy,
					//		since the O and H positions are only sampled without minimization
					minEnergy = Math.min(minEnergy,curEnergy);
					if (curSample==0)
						maxEnergy = Math.max(maxEnergy,curEnergy);
	
					m.updateCoordinates();//restore the actualCoordinates array to the initial values
				}
			}
		}
		else if (computeEVEnergy){
			a96ff.initializeCalculation();
			a96ff.setNBEval(hElect,hVDW);
			curEnergy = calcTotalSnapshotEnergy();
			m.updateCoordinates();
			minEnergy = curEnergy;
			maxEnergy = curEnergy;
		}

		// Store result
		int indx1;
		if (!isLig){ //AS residue
			if (totRotForCur==0) //no rotamers
				indx1 = curPos*numTotalRotamers + rotamerIndexOffset[curAA] + 1;
			else //the index of the current rotamer
				indx1 = curPos*numTotalRotamers + rotamerIndexOffset[curAA] + curRot + 1;
		}
		else { //ligand
			if (totRotForCur==0) //no rotamers
				indx1 = numInAS * numTotalRotamers + 1;
			else //the index of the current rotamer
				indx1 = numInAS * numTotalRotamers + curRot + 1;
		}
		
		intraEnergyArrayMin[indx1] = minEnergy;
		intraEnergyArrayMax[indx1] = maxEnergy;
	}
	
	// Sets up the computation for residue-to-template (SHL-AS, LIG-SHL) and rotamer-rotamer (AS-AS, LIG-AS) energies;
	// Searches among amino acid types for all mutatable residues (as determined by residueMutatable);
	//		Sets all non-mutatable residues in residueMap to either Gly or Ala (and leaves all Pro)
	public void pairwiseEnergyComp (final int maxDepth,
		int residueMap[], String residuesAllowed[],
		int numResAllowed, int rotamerIndexOffset[], int numTotalRotamers,
		boolean ligPresent, int residueMutatable[],
		float retEMatrix[][], int matrixSizeRows, boolean shellRun, boolean minimizeBB) {
		
		
		if (ligPresent) { //there is a ligand and its energies should be computed
			if(ligROT.getNumAllowable(0)==0) {  // this shouldn't happen
				System.out.println("ERROR: Ligand has no allowables but you are using a ligand?");
				System.exit(1);
			}
			else {
				// Because this is a ligand there can be only one allowable type, ie. it can not mutate;
				// Change the residue type to this one type
				curLigNum = ligROT.getIndexOfNthAllowable(0,0);
				ligROT.changeResidueType(m,0,rl.getAAName(curLigNum),addHydrogens);
				m.strand[ligStrNum].residue[0].flexible = true;
				m.strand[ligStrNum].residue[0].setEnergyEval(true, true);
			}
		}
		
		//Check with all allowed amino acids for all of the mutatable positions in residueMutatable[];
		//	all non-mutatable positions in residueMutatable (these are necessarily in residueMap[]) are set ot Gly or Ala (or remain Pro)
		
		int numMut = 0; //should be either 1 (shellRun) or 2 (rotamer-rotamer energies)
		int mutDepth[] = new int[2];
		for (int i=0; i<mutDepth.length; i++)
			mutDepth[i] = -1;
		
		for (int depth=0; depth<maxDepth; depth++){	
			if (residueMutatable[depth]!=0){
				mutDepth[numMut] = depth;
				numMut++;
			}
		}
		
		if ( ( (!ligPresent) && ( ((shellRun) && (numMut!=1)) || ((!shellRun) && (numMut!=2)) ) )
				|| ( (ligPresent) && ( ((shellRun) && (numMut!=0)) || ((!shellRun) && (numMut!=1)) ) ) ) {
			System.out.println("ERROR: incorrect number of mutatable positions in residueMutatable[]");
			System.exit(1);
		}
			
		for (int depth=0; depth<maxDepth; depth++){				
			if (residueMutatable[depth]==0){ //not a mutatable residue, but in residueMap[]
				
				// Make this residue a "gly", if it is not already GLY or PRO;
				//		If minimizeBB, then change to ALA, if not already GLY or PRO
				String a = rl.getThreeLetAACode(sysLR.getCurRotType(residueMap[depth]));
				if ((!a.equalsIgnoreCase("GLY"))&&(!a.equalsIgnoreCase("PRO"))){
					
					String tmpAA = "gly";
					if (minimizeBB)
						tmpAA = "ala";
				
					sysLR.changeResidueType(m,residueMap[depth],tmpAA,addHydrogens,connectResidues);
				}
				
				curAANum[residueMap[depth]] = -1;
				m.strand[sysStrNum].residue[residueMap[depth]].flexible = false;
				m.strand[sysStrNum].residue[residueMap[depth]].setEnergyEval(false, shellRun);
			}
		}
		
		pairwiseEnergyCompAllMutatedResHelper(maxDepth,	residueMap,	residuesAllowed, numResAllowed, rotamerIndexOffset,
				numTotalRotamers, -1, -1, ligPresent, residueMutatable, retEMatrix, matrixSizeRows, shellRun, minimizeBB, numMut, mutDepth, 0);
	}
	
	// Helper to pairwiseEnergyComp();
	// Searches among amino acid types for all mutatable residues (as determined by residueMutatable[]);
	//		the non-mutatable residues in residueMutatable[] have already been set to either Gly or Ala (all Pro remain)
	private void pairwiseEnergyCompAllMutatedResHelper(int maxDepth,
			int residueMap[], String residuesAllowed[],
			int numResAllowed, int rotamerIndexOffset[], int numTotalRotamers,
			int res1Num, int res2Num,boolean ligPresent, int residueMutatable[],
			float retEMatrix[][], int matrixSizeRows, boolean shellRun, boolean minimizeBB, final int numMut, int mutDepth[], int curMut){
		
		if (curMut>=numMut){// If we've arrived here, then we have assigned a full amino acid sequence and we are ready to generate all rotamer combinations
			
			if(debug){
				System.out.print(".");
			}
			if (computeEVEnergy){
				a96ff.calculateTypesWithTemplates();
				a96ff.initializeCalculation();
				a96ff.setNBEval(hElect,hVDW);
				if (doMinimization){
					if (!minimizeBB){ //side-chain minimization
						if(ligPresent)
							simpMin.initialize(m,sysStrNum,ligStrNum,a96ff,sysLR,ligROT,curAANum,curLigNum,doDihedE,rl);
						else
							simpMin.initialize(m,sysStrNum,a96ff,sysLR,curAANum,doDihedE,rl);
					}
					else { //backbone minimization
						if (ligPresent)
							bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
						else
							bbMin.initialize(m, a96ff, residueMap, sysStrNum);
					}
				}
			}
			int res1 = -1, res2 = -1;
			for(int i=0;i<maxDepth;i++){
				if (residueMutatable[i]==1) {
					if (res1 == -1)
						res1 = i;
					else
						res2 = i;
				}
			}
			
			if (ligPresent){
				// The first residue is the ligand (so -2 is passed), the second is numerically res1
				pairwiseMutationAllRotamerSearch(true, maxDepth, residueMap, -2, res1, -2,
					res1Num,-1,-1,rotamerIndexOffset,numTotalRotamers,
					ligPresent, retEMatrix, matrixSizeRows, minimizeBB, numMut, mutDepth, 0);
			}
			else{
				pairwiseMutationAllRotamerSearch(false, maxDepth, residueMap, res1, res2, res1Num,
					res2Num,-1,-1,rotamerIndexOffset,numTotalRotamers,
					ligPresent,retEMatrix, matrixSizeRows, minimizeBB, numMut, mutDepth, 0);
			}
			
			return;
		}
		else {
			boolean isFirstAA = (res1Num == -1);
			int depth = mutDepth[curMut];
			for(int q=0;q<sysLR.getNumAllowable(residueMap[depth]);q++) {
				
				int AAindex = sysLR.getIndexOfNthAllowable(residueMap[depth],q);
				// Change the residue type
				if (isFirstAA)
					res1Num = AAindex;
				else
					res2Num = AAindex;
				curAANum[residueMap[depth]] = AAindex;
				
				// Apply mutation
				sysLR.changeResidueType(m,residueMap[depth],rl.getAAName(curAANum[residueMap[depth]]),addHydrogens,connectResidues);
				
				m.strand[sysStrNum].residue[residueMap[depth]].flexible = true;
				m.strand[sysStrNum].residue[residueMap[depth]].setEnergyEval(true, shellRun);
				pairwiseEnergyCompAllMutatedResHelper(maxDepth,
					residueMap,	residuesAllowed, numResAllowed, rotamerIndexOffset,
					numTotalRotamers,res1Num,res2Num, ligPresent, residueMutatable, 
					retEMatrix, matrixSizeRows, shellRun, minimizeBB, numMut, mutDepth, curMut+1);
			}
		}
	}

	// Called by pairwiseEnergyCompAllMutatedResHelper();
	// Computes the energies among all rotamer combinations for residues res1 and res2 (if res2==1, then computes the res1-to-template energies)
	private void pairwiseMutationAllRotamerSearch(boolean isLigRun, int maxDepth, int residueMap[],
		int res1, int res2, int res1AANum, int res2AANum, int res1RotNum, int res2RotNum,
		int rotamerIndexOffset[], int numTotalRotamers, boolean ligPresent, float retEMatrix[][], 
		int matrixSizeRows, boolean minimizeBB, final int numMut, int mutDepth[], int curMut) {

		// If we're at the ligand depth
		if (isLigRun) {
			isLigRun = false;
			int numLigRot = rl.getNumRotForAAtype(curLigNum);
			if (numLigRot==0) { //ligand is ALA or GLY
				pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
						res1AANum, res2AANum, 0, res2RotNum, rotamerIndexOffset,
						numTotalRotamers, ligPresent, retEMatrix, matrixSizeRows, minimizeBB, numMut, mutDepth, curMut);
			}
			else {
				for(int w=0;w<rl.getNumRotForAAtype(curLigNum);w++){
					ligROT.applyRotamer(m, 0, w);
					pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
						res1AANum, res2AANum, w, res2RotNum, rotamerIndexOffset,
						numTotalRotamers, ligPresent, retEMatrix, matrixSizeRows, minimizeBB, numMut, mutDepth, curMut);
				}
			}
		}
		else {
			if (curMut>=numMut) {
				// Index jumps to the right rotamer block with res1*numTotalRotamers
				//  then jumps to the right offset within this block with rotamerIndexOffset[]
				//  then the rotamer number is added res1AANum and finally 1 is added because
				//  the 0 position is reserved for interactions made with the backbone
				int res1Strand = 0;
				int res1Resnum = 0;
				int indx1 = -1, indx2 = -1;
				boolean shellRun = false;
				
				if (res1 == -2) {
					// res1 is the ligand
					indx1 = maxDepth * numTotalRotamers + res1RotNum + 1;
					res1Strand = ligStrNum;
					res1Resnum = 0;
				}
				else {
					indx1 = res1 * numTotalRotamers + rotamerIndexOffset[res1AANum] + res1RotNum + 1;
					res1Strand = sysStrNum;
					res1Resnum = residueMap[res1];	
				}

				if (res2 == -1) {
					// then we're doing a shell search
					indx2 = 0;
					shellRun = true;
				}	
				else {
					indx2 = res2 * numTotalRotamers + rotamerIndexOffset[res2AANum] + res2RotNum + 1;
				}
				
				// If we've gotten here then check the sterics to make sure
				//  we're sterically allowable
				boolean stericallyGood = false;
				if (hSteric) {
					if (RS_CheckAllStericsWithH(res1Strand,res1Resnum)) {
						if (res2 != -1) {
							if (RS_CheckAllStericsWithH(sysStrNum,residueMap[res2]))
								stericallyGood = true;
						} else 
							stericallyGood = true;
					}
				}
				else {
					if (RS_CheckAllSterics(res1Strand,res1Resnum)) {
						if (res2 != -1) {
							if (RS_CheckAllSterics(sysStrNum,residueMap[res2]))
								stericallyGood = true;
						} else 
							stericallyGood = true;
					}
				}
				
				if ((stericallyGood)) { //good steric found or doing backbone minimization
					
					// After minimization do a m.updateCoordinates() to resync the actualCoordinates which were changed
					//  in the minimization procedure									
					
					if (doMinimization){
						if (!minimizeBB) //side-chain minimization
							pairwiseRotamerEnergySidechainHelper (res1,res2,res1Strand,res1Resnum,residueMap,indx1,indx2,
									shellRun, retEMatrix, matrixSizeRows);
						else //backbone minimization
							pairwiseRotamerEnergyBackboneHelper (res1,res2,res1Strand,res1Resnum,residueMap,indx1,indx2,
									shellRun, retEMatrix, matrixSizeRows, false);
					}
					else if (computeEVEnergy){
						
						a96ff.initializeCalculation();
						a96ff.setNBEval(hElect,hVDW);
						float curEnergy = 0.0f;	
						curEnergy = calcTotalSnapshotEnergy();
						
						//Remove the intra-residue energies here
						float tmpE = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber);
						curEnergy -= tmpE;
						if (res2 != -1) {
							tmpE = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[res2]].moleculeResidueNumber);
							curEnergy -= tmpE;
						}
						else {
							// If this is a shell run then subtract off the shell-to-shell
							//  interaction energies so that they are not over counted
							if (!minimizeBB)
								curEnergy -= retEMatrix[0][0];//same as mutationEnergiesMin[0][0]
							else
								curEnergy -= computeEnergyOfOnlyTemplate(residueMap);
						}
						
						//Store the computed energy
						if (shellRun) {
							retEMatrix[0][indx1] = curEnergy;
							retEMatrix[0][matrixSizeRows+indx1] = curEnergy;
						}
						else {
							retEMatrix[indx1][indx2] = curEnergy;
							retEMatrix[indx2][indx1] = curEnergy;
							retEMatrix[indx1][matrixSizeRows+indx2] = curEnergy;
							retEMatrix[indx2][matrixSizeRows+indx1] = curEnergy;
						}
					}
					else {
						System.out.println("This should not happen. No energy evaluation specified");
						System.exit(1);
					}					
															
					//restore to the coordinates before the energy computation
					m.updateCoordinates();
				}
				else {
					if (shellRun) {
						retEMatrix[0][indx1] = stericE;
						retEMatrix[0][matrixSizeRows+indx1] = stericE;
					}
					else {
						retEMatrix[indx1][indx2] = stericE;
						retEMatrix[indx2][indx1] = stericE;
						retEMatrix[indx1][matrixSizeRows+indx2] = stericE;
						retEMatrix[indx2][matrixSizeRows+indx1] = stericE;
					}
				}

				return;
			}

			else {
				// Next check with different rotamers
				int depth = mutDepth[curMut];
				int numRot = rl.getNumRotForAAtype(curAANum[residueMap[depth]]);
				
				if (numRot == 0) {//If there are no rotamers for this AA then allow the default conformation; this will only happen with Ala and Gly
					if (depth==res1)
						pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
							res1AANum, res2AANum, 0, res2RotNum, rotamerIndexOffset,
							numTotalRotamers, ligPresent, retEMatrix, matrixSizeRows, minimizeBB,
							numMut, mutDepth, curMut+1);
					else
						pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
							res1AANum, res2AANum, res1RotNum, 0, rotamerIndexOffset,
							numTotalRotamers, ligPresent, retEMatrix, matrixSizeRows, minimizeBB,
							numMut, mutDepth, curMut+1);
				}
				else {		
					for(int w=0;w<numRot;w++) {					
						sysLR.applyRotamer(m, residueMap[depth], w);
						if (depth==res1)
							pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
								res1AANum, res2AANum, w, res2RotNum, rotamerIndexOffset,
								numTotalRotamers,ligPresent, retEMatrix, matrixSizeRows, minimizeBB, 
								numMut, mutDepth, curMut+1);
						else
							pairwiseMutationAllRotamerSearch(false,maxDepth,residueMap,res1,res2,
								res1AANum, res2AANum, res1RotNum, w, rotamerIndexOffset,
								numTotalRotamers, ligPresent, retEMatrix, matrixSizeRows, minimizeBB, 
								numMut, mutDepth, curMut+1);
					}
				}
			}
		}
	}	
	
	//This method computes the minimized energy for a given pair of rotamers (or rot-to-template)
	//Called by pairwiseMutationAllRotamerSearch(.)
	private void pairwiseRotamerEnergySidechainHelper (int res1,int res2,int res1Strand,int res1Resnum,int residueMap[],
			int indx1,int indx2,boolean shellRun,float retEMatrix[][],int matrixSizeRows){
		
		//Initialize Amber for the current pair
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
		
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
		
		float curEnergy = 0.0f;
		float beginE = 0.0f;
		
		/////////////////////////////
		//formally making sure that when minimization is performed,
		//		the minimized value is never greater than the initial value
		beginE = calcTotalSnapshotEnergy();
		
		//Remove the intra-residue energies here
		float tmpE1 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber);
		beginE -= tmpE1;
		if (res2 != -1) {
			tmpE1 = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[res2]].moleculeResidueNumber);
			beginE -= tmpE1;
		}
		else {
			// If this is a shell run then subtract off the shell-to-shell
			//  interaction energies so that they are not over counted
			beginE -= retEMatrix[0][0];
		}
		/////////////////////////////
		
		//Minimize
		if ((res1 == -2) && (shellRun))
			simpMin.minimize(numMinSteps,true); //ligand only
		else
			simpMin.minimize(numMinSteps,false);
		
		curEnergy = calcTotalSnapshotEnergy(); //compute the energy after minimization
		
		//Remove the intra-residue energies here
		float tmpE2 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber);
		curEnergy -= tmpE2;
		if (res2 != -1) {
			tmpE2 = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[res2]].moleculeResidueNumber);
			curEnergy -= tmpE2;
		}
		else {
			// If this is a shell run then subtract off the shell-to-shell
			//  interaction energies so that they are not over counted
			curEnergy -= retEMatrix[0][0];//same as mutationEnergiesMin[0][0]
		}
		/////////////////////////////

		//Compare to the min and max energies found so far and update, if necessary
		if (beginE<curEnergy)
			curEnergy = beginE;
		float lE = Math.min(beginE, curEnergy);
		float hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);		
								
		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		
		if (shellRun) {
			retEMatrix[0][indx1] = minEnergy;
			retEMatrix[0][matrixSizeRows+indx1] = maxEnergy;
		}
		else {
			retEMatrix[indx1][indx2] = minEnergy;
			retEMatrix[indx2][indx1] = minEnergy;
			retEMatrix[indx1][matrixSizeRows+indx2] = maxEnergy;
			retEMatrix[indx2][matrixSizeRows+indx1] = maxEnergy;
		}
		
		return;
	}
	
	//This method helps perform random sampling for a given pair of rotamers for backbone minimization
	//Called by pairwiseMutationRotamerSearch(.)
	private void pairwiseRotamerEnergyBackboneHelper (int res1,int res2,int res1Strand,int res1Resnum,int residueMap[],
			int indx1,int indx2,boolean shellRun,float retEMatrix[][],int matrixSizeRows, boolean templateOnly){		

		//Initialize Amber for the current pair
		a96ff.initializeCalculation();
		a96ff.setNBEval(hElect,hVDW);
	
		float minEnergy = (float)Math.pow(10,30);
		float maxEnergy = -(float)Math.pow(10,30);
			
		float beginE;
		if ((shellRun)&&(!templateOnly)) //res-to-template with backbone flexibility
			beginE = calcPartEnergy(0);
		else 
			beginE = calcTotalSnapshotEnergy();
		
		//Remove the intra-residue energies here
		if (!templateOnly){
			float tmpE1 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber);
			beginE -= tmpE1;
			if (res2 != -1) {
				tmpE1 = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[res2]].moleculeResidueNumber);
				beginE -= tmpE1;
			}
		}
		
		bbMin.minimizeFull(true); //minimize
		
		float curEnergy;
		if ((shellRun)&&(!templateOnly)) //res-to-template with backbone flexibility
			curEnergy = calcPartEnergy(0);
		else 
			curEnergy = calcTotalSnapshotEnergy();
		
		//Remove the intra-residue energies here
		if (!templateOnly){
			float tmpE1 = computeEnergyOfOnlyRes(m.strand[res1Strand].residue[res1Resnum].moleculeResidueNumber);
			curEnergy -= tmpE1;
			if (res2 != -1) {
				tmpE1 = computeEnergyOfOnlyRes(m.strand[sysStrNum].residue[residueMap[res2]].moleculeResidueNumber);
				curEnergy -= tmpE1;
			}
		}
		
		if (beginE<curEnergy)
			curEnergy = beginE;
		float lE = Math.min(beginE, curEnergy);
		float hE = Math.max(beginE, curEnergy);
		minEnergy = Math.min(minEnergy,lE);
		maxEnergy = Math.max(maxEnergy,hE);	
								
		m.updateCoordinates();//restore the actualCoordinates array to the initial values
		/////////////////////////////////////////////
		
		if (templateOnly){
			retEMatrix[0][0] = minEnergy;
			retEMatrix[0][matrixSizeRows] = maxEnergy;
		}
		else if (shellRun) {
			retEMatrix[0][indx1] = minEnergy;
			retEMatrix[0][matrixSizeRows+indx1] = maxEnergy;
		}
		else {
			retEMatrix[indx1][indx2] = minEnergy;
			retEMatrix[indx2][indx1] = minEnergy;
			retEMatrix[indx1][matrixSizeRows+indx2] = maxEnergy;
			retEMatrix[indx2][matrixSizeRows+indx1] = maxEnergy;
		}
		
		return;
	}
//// END PAIRWISE MUTATION _ALL_ ROTAMER SEARCH SECTION - ENERGY PRECOMPUTATION
///////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
//// BEGIN ROTAMER SEARCH SECTION
	
	//For the given mutation sequence,
	//Count the total number of conformations, the number of confs pruned by MinDEE, the
	//	number of remaining conformations; the number of conformations above each level
	//	(from level i+1 to the top level, which is the ligand); the number of rotamers for each flexible residue;
	//	the number of rotamers (total and non-pruned by MinDEE) for the flexible residues only;
	//The MinDEE matrix eliminatedRotAtRes[] should already be computed;
	//For each pruned rotamer at curLevel, count the number of *new* conformations that 
	//	are pruned as a result, then decrease the number of rotamers for curLevel to
	//	make sure we do not overcount the pruned conformations (the next pruned rotamer
	//	should not count the conformations that include the first rotamer, and therefore
	//	have already been pruned)
	private int countConfs(int numInAS, int numRes, int numTotalRotamers, int rotamerIndexOffset[], boolean ligPresent, 
			int residueMap[], int numRotForRes[], int numRotForResNonPruned[]){
		
		int numTotalRotRed = 0;
		numConfsTotal = BigInteger.valueOf(1);
		numConfsAboveLevel = new BigInteger[numRes];

		int curNumRot;
		
		//Store the number of rotamers for each AA in the current mutation
		//The ligand (if present) is in the last level
		for (int curLevel=0; curLevel<numRes; curLevel++){
			if ((ligPresent)&&(curLevel==(numRes-1))) //the ligand level
				curNumRot = rl.getNumRotForAAtype(curLigNum);
			else
				curNumRot = rl.getNumRotForAAtype(curAANum[residueMap[curLevel]]);
			if (curNumRot==0) //GLY or ALA
				curNumRot = 1;
			numRotForRes[curLevel] = curNumRot;
			numRotForResNonPruned[curLevel] = numRotForRes[curLevel];
			numTotalRotRed += numRotForRes[curLevel];
			numConfsTotal = numConfsTotal.multiply(BigInteger.valueOf(numRotForRes[curLevel]));
		}
		
		BigInteger numPruned = BigInteger.ZERO;
		int numPrunedThisLevel;
		
		//Count the number of rotamers pruned by MinDEE
		for (int curLevel=0; curLevel<numRes; curLevel++){
			numPrunedThisLevel = 0;
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				int curIndex;
				if ((ligPresent)&&(curLevel==(numRes-1))) //the ligand level
					curIndex = numInAS*numTotalRotamers + curRot;
				else
					curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[curAANum[residueMap[curLevel]]] + curRot;
				
				if (eliminatedRotAtRes[curIndex]){
					numPruned = numPruned.add(compPrunedConfsByRot(numRotForResNonPruned,numRes,curLevel));
					numPrunedThisLevel++;
				}
			}
			numRotForResNonPruned[curLevel] -= numPrunedThisLevel;
		}
		
		//Count the number of conformations below each level (flexible residue)
		numConfsAboveLevel[numRes-1] = new BigInteger("1"); //the last level
		if (numRes>1){
			for (int curLevel=numRes-2; curLevel>=0; curLevel--){
				numConfsAboveLevel[curLevel] = numConfsAboveLevel[curLevel+1].multiply(BigInteger.valueOf(numRotForResNonPruned[curLevel+1]));
			}
		}
		
		numConfsPrunedByMinDEE = numPruned; //set the number of confs pruned by MinDEE
		numConfsLeft = numConfsTotal.subtract(numConfsPrunedByMinDEE);
		numConfsPrunedByS = BigInteger.valueOf(0);
		
		return numTotalRotRed;
	}
	
	//Computes the number of conformations pruned by a rotamer (eliminated by MinDEE) at curLevel:
	//	multiply the number of rotamers for each level different from curLevel
	private BigInteger compPrunedConfsByRot(int numRotForRes[], int numRes, int curLevel){
		
		BigInteger numPruned = BigInteger.ONE;
		
		for (int i=0; i<numRes; i++){
			if (i!=curLevel)
				numPruned = numPruned.multiply(BigInteger.valueOf(numRotForRes[i]));
		}
		
		return numPruned;
	}
	
	//Computes the number of conformations that are pruned by MinDEE due to steric clashes:
	//	for each level, count the number of rotamers that are not pruned by MinDEE due to steric clash, then 
	//	multiply for all levels; the result is the number of conformations that do not include
	//	pruned rotamers; return (totalNumConfs - this number) as the number of conformations pruned by MinDEE
	private BigInteger countPrunedByMinDEESteric(int numInAS, int numRes, int numTotalRotamers, int rotamerIndexOffset[], boolean ligPresent, 
			int residueMap[], int numRotForRes[], boolean prunedIsSteric[]){
		
		BigInteger numConfNotPruned = BigInteger.ONE;
		int numNonPruned[] = new int[numRes];
		for (int i=0; i<numRotForRes.length; i++)
			numNonPruned[i] = numRotForRes[i];
		
		for (int curLevel=0; curLevel<numRes; curLevel++){
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				int curIndex;
				if ((ligPresent)&&(curLevel==(numRes-1))) //the ligand level
					curIndex = numInAS*numTotalRotamers + curRot;
				else
					curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[curAANum[residueMap[curLevel]]] + curRot;
				if ((eliminatedRotAtRes[curIndex])&&(prunedIsSteric[curIndex])){
					numNonPruned[curLevel]--;
				}
			}
		}
		
		for (int curLevel=0; curLevel<numRes; curLevel++)
			numConfNotPruned = numConfNotPruned.multiply(BigInteger.valueOf(numNonPruned[curLevel]));
		
		return (numConfsTotal.subtract(numConfNotPruned));
	}
	
	//Sets-up the repeat mutation search if the desired accuracy has not been achived, although
	//	all non-pruned conformations have been generated by A*;
	//Updates eliminatedRotAtRes[] to not prune a subset of the pruned rotamers, such that the
	//	number of conformations not pruned by this is at least numPrunedConfsToAllow
	private void setupRepeatRun(BigInteger numPrunedConfsToAllow, int numRotForResNonPruned[], 
			int numInAS){
		
		BigInteger totalNumConfsUnpruned = new BigInteger("0");
		for (int i=0; i<eliminatedRotAtRes.length; i++){
			if ((eliminatedRotAtRes[i])&&(!prunedIsSteric[i])){ //pruned non-steric rotamer
				
				int curLevel = (int)Math.floor(i/152); //the residue number for the cur rot index
				numRotForResNonPruned[curLevel]++;
				eliminatedRotAtRes[i] = false;
				
				BigInteger numConfsAtLevelUnpruned = new BigInteger("1"); //count the unpruned confs by unpruning the cur rot index
				for (int j=0; j<numInAS; j++){
					if (j!=curLevel){
						numConfsAtLevelUnpruned = numConfsAtLevelUnpruned.multiply(BigInteger.valueOf(numRotForResNonPruned[j]));
					}
				}
				
				totalNumConfsUnpruned = totalNumConfsUnpruned.add(numConfsAtLevelUnpruned);
				if (totalNumConfsUnpruned.compareTo(numPrunedConfsToAllow)>=0) // num pruned confs reduced by the necessary number
					break;
			}
		}
	}
	
	// This function performs a rotamer search to compute
	//  a partition function for the slave node.
	// Utilizes a number of helper functions
	public void slaveDoRotamerSearch(boolean searchComputeEVEnergy,
		boolean searchDoMinimization, int numInAS, int numAAAllowed,
		int numTotalRotamers, int rotamerIndexOffset[], String AAList[],
		int residueMap[], boolean usingInitialBest, BigDecimal initialBest,
		CommucObj cObj, boolean minimizeBB, boolean saveConfs, String fName) {

		// A rotamer search is performed. For each residue,
		//  every allowable rotamer is tried in combination
		//  with every other allowable rotamer

		computeEVEnergy = searchComputeEVEnergy;
		doMinimization = searchDoMinimization;
		ASAANums = new int[numInAS];
		curASRotNum = new int[numInAS];
		int curResToASMap[] = new int[m.strand[sysStrNum].numberOfResidues];
			// This map maps the system strand residues back to the AS numbering
			// So something like 8 -> 0, 10 -> 1, 11 -> 2, ...
		curLigRotNum = 0;
		numConfsPrunedByE = BigInteger.ZERO;
		numConfsPrunedByS = BigInteger.ZERO;
		numConfsEvaluated = BigInteger.ZERO;
		numConfsPrunedByMinDEE = BigInteger.ZERO;
		allPruned = false;
		//confEnergies = new Vector(2048,1024);
		bestEnergyFound = 9999999.0f;
		bestEMin = 9999999.0f;
		bestEUnMin = 9999999.0f;
		if (usingInitialBest)
			initial_q = initialBest.multiply(new BigDecimal((double)(1-KSepsilon)));
		else
			initial_q = new BigDecimal(0.0);
		partial_q = new BigDecimal(0.0);
		partial_p = new BigDecimal(0.0);
		
		for (int i=0; i<numInAS; i++) //the AS residues are flexible - this is used by simpMin to set up the minimizer
			m.residue[residueMap[i]].flexible = true;
		if (haveLigand) //the ligand is also flexible
			m.strand[ligStrNum].residue[0].flexible = true;
		

		if (searchDoMinimization && !searchComputeEVEnergy){
			System.out.println("Warning: In order to do minimization computeEVEnergy must be true");
			return;
		}

		if (haveLigand && !ligandReady){
			System.out.println("Warning: Ligand rotamers not setup, RotamerSearch aborting");
			return;
		}

		// Prepare Amber
		if(searchComputeEVEnergy){
			// Amber should already be loaded
			if(haveLigand) {
				a96ff.setLigandNum(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
			}
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null) {
						System.out.println("Warning: Attempting minimization run but simpMin not allocated, RotamerSearch aborting");
						return;
					}
					bbMin = null;
				}
				else { //backbone minimization
					if (bbMin == null) {
						System.out.println("Warning: Attempting minimization run but bbMin not allocated, RotamerSearch aborting");
						return;
					}
					simpMin = null;
				}
			}
		}

		// Make sure the allRotamerPairsEnergyName matrices exist
		if (arpMatrix == null){
			System.out.println("Warning: allRotamerPairsEnergy matrix not loaded");
			return;
		}
		
		if (eliminatedRotAtRes == null){
			System.out.println("Warning: MinDEE matrix not computed");
			return;
		}

		// Setup the residue number to AS number map
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			curResToASMap[i] = -1;
		}
		for(int i=0;i<residueMap.length;i++){
			curResToASMap[residueMap[i]] = i;
		}

		if (haveLigand)
			slaveMutationRotamerSearch(-1, numInAS, residueMap, AAList, numAAAllowed,
				numTotalRotamers, rotamerIndexOffset,	curResToASMap, haveLigand, minimizeBB, saveConfs, fName);
		else
			slaveMutationRotamerSearch(0, numInAS, residueMap, AAList, numAAAllowed,
				numTotalRotamers, rotamerIndexOffset, curResToASMap, haveLigand, minimizeBB, saveConfs, fName);
	
		// Store results to the communication object
		if (cObj!=null){
			if (haveLigand) {
				cObj.EL_searchNumConfsTotal = numConfsTotal.intValue();
				if (allPruned) { //all of the conformations were pruned by MinDEE, as there is a residue with no remaining rotamers
					cObj.EL_searchNumConfsPrunedByS = 0;
					cObj.EL_searchNumPrunedMinDEE = numConfsTotal.intValue();
					cObj.EL_searchNumConfsEvaluated = 0;
					cObj.EL_searchNumConfsLeft = 0;
				}
				else {
					cObj.EL_searchNumConfsPrunedByS = numConfsPrunedByS.intValue();
					cObj.EL_searchNumPrunedMinDEE = numConfsPrunedByMinDEE.intValue();
					cObj.EL_searchNumConfsEvaluated = numConfsEvaluated.intValue();
					cObj.EL_searchNumConfsLeft = numConfsLeft.intValue();
				}
				cObj.EL_searchBestEnergyFound = bestEnergyFound;
			} else {
				cObj.E_searchNumConfsTotal = numConfsTotal.intValue();
				if (allPruned) { //all of the conformations were pruned by MinDEE, as there is a residue with no remaining rotamers
					cObj.E_searchNumConfsPrunedByS = 0;
					cObj.E_searchNumPrunedMinDEE = numConfsTotal.intValue();
					cObj.E_searchNumConfsEvaluated = 0;
					cObj.E_searchNumConfsLeft = 0;
				}
				else {
					cObj.E_searchNumConfsPrunedByS = numConfsPrunedByS.intValue();
					cObj.E_searchNumPrunedMinDEE = numConfsPrunedByMinDEE.intValue();
					cObj.E_searchNumConfsEvaluated = numConfsEvaluated.intValue();
					cObj.E_searchNumConfsLeft = numConfsLeft.intValue();
				}
				cObj.E_searchBestEnergyFound = bestEnergyFound;
			}
	
			if (haveLigand)
				cObj.EL_allPruned = allPruned;
			else
				cObj.E_allPruned = allPruned;
			
			// Compute q_X
			if (haveLigand) {
				cObj.q_EL = partial_q;
				cObj.bestBoundE = (double)bestEUnMin;
				cObj.bestBoundEMin = (double)bestEMin;
			} else {
				cObj.q_E = partial_q;
				cObj.bestUnBoundE = (double)bestEUnMin;
				cObj.bestUnBoundEMin = (double)bestEMin;
			}
		}
		else {
			if (haveLigand)
				System.out.println("Statistics (bound):");
			else
				System.out.println("Statistics (unbound):");
			System.out.println("Best Energy:  " + bestEnergyFound);
			System.out.println("partial_q: " + partial_q);
			System.out.println("partial_p: " + partial_p);
			System.out.println("NumConfsTotal:      	" + numConfsTotal);
			System.out.println("NumConfsPrunedByMinDEE: " + numConfsPrunedByMinDEE);
			System.out.println("NumConfsPrunedByS:  	" + numConfsPrunedByS);
			System.out.println("NumConfsEvaluated:  	" + numConfsEvaluated);
			System.out.println("NumConfsLeft:       	" + numConfsLeft);
		}
		
	}


	// This function does a mutation search then for each allowable mutation
	//  a simple rotamer search is performed
	private void slaveMutationRotamerSearch(int depth, int maxDepth, int 
		residueMap[], String AAList[], int numAAAllowed, int numTotalRotamers,
		int rotamerIndexOffset[], int curResToASMap[], boolean ligPresent, boolean minimizeBB,
		boolean saveConfs, String fName) {
	
		// If we're at the ligand depth
		if (depth == -1) {
			if(ligROT.getNumAllowable(0)==0) {
				curLigNum = -1;  // this shouldn't happen
				System.out.println("Error: Ligand has no allowables but you are using a ligand?");
				slaveMutationRotamerSearch(depth+1,maxDepth,residueMap,AAList,numAAAllowed,
					numTotalRotamers,rotamerIndexOffset,curResToASMap,ligPresent,minimizeBB,saveConfs,fName);
			}
			else {
				// Because this is a ligand there can be only one
				//  allowable type, ie. it can not mutate
				// Change the residue type to this one type
				curLigNum = ligROT.getIndexOfNthAllowable(0,0);
				ligROT.changeResidueType(m,0,rl.getAAName(curLigNum),addHydrogens);
				slaveMutationRotamerSearch(depth+1,maxDepth,residueMap,AAList,numAAAllowed,
					numTotalRotamers,rotamerIndexOffset,curResToASMap,ligPresent,minimizeBB,saveConfs,fName);
			}
		}
		else {
			if (depth >= maxDepth) {
				// If we've arrived here then we're ready to
				//  do a rotamerSearch
				if(debug)
					System.out.println("One Mutation Conf Found and is being tested");
				if (computeEVEnergy){
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					if (doMinimization){
						if (!minimizeBB){ //side-chain minimization
							if(ligPresent)
								simpMin.initialize(m,sysStrNum,ligStrNum,a96ff,sysLR,ligROT,curAANum,curLigNum,doDihedE,rl);
							else
								simpMin.initialize(m,sysStrNum,a96ff,sysLR,curAANum,doDihedE,rl);
						}
						else { //backbone minimization
							if (ligPresent)
								bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
							else
								bbMin.initialize(m, a96ff, residueMap, sysStrNum);
						}
					}
				}

				// Determine the following
				// -AANums of all AS residues (so that we can index properly into rotamerIndexOffset, others)
				//   residue i (i=0..8) has 3 letter code AAList(AANums[i])
				// -total number of conformations for this mutation -> assign to numConfsTotal, numConfsLeft
				// -number of conformations below each level (ie. if we prune one of the 3rd AAs rotamers
				//   because of steric overlap, how many confs are we pruning)
				assignAANums(residueMap,AAList,numAAAllowed);
				//int numConfs = computeTotalNumConfs(residueMap,curResToASMap,ligPresent);
				
				  // also computes num of confs at each level
				//numConfsLeft = (numConfsTotal-numConfsPrunedByMinDEE); //MinDEE should have already been applied
				
				//Perform A* search: compute the partial partition function q*
				slaveRotamerSearchAStar(ligPresent, maxDepth, numTotalRotamers, rotamerIndexOffset, curResToASMap,
						residueMap, minimizeBB, saveConfs, fName);
				
				return;
			}

			// If there are no allowables then test with 'native' form of the enzyme
			if (sysLR.getNumAllowable(residueMap[depth]) == 0) {
				curAANum[residueMap[depth]] = -1;
				slaveMutationRotamerSearch(depth+1,maxDepth,residueMap,AAList,numAAAllowed,
					numTotalRotamers,rotamerIndexOffset,curResToASMap,ligPresent,minimizeBB,saveConfs,fName);
			}
			else {
				// Otherwise check with different AAs
				for(int q=0;q<sysLR.getNumAllowable(residueMap[depth]);q++) {
					curAANum[residueMap[depth]] = sysLR.getIndexOfNthAllowable(residueMap[depth],q);
					// Change the residue type
					sysLR.changeResidueType(m,residueMap[depth],rl.getAAName(curAANum[residueMap[depth]]),addHydrogens,connectResidues);
					slaveMutationRotamerSearch(depth+1,maxDepth,residueMap,AAList,numAAAllowed,
						numTotalRotamers,rotamerIndexOffset,curResToASMap,ligPresent,minimizeBB,saveConfs,fName);
				}
			}
		}
	}

	// Calls AStar repeatedly while the returned conformations still have energy below the threshold;
	//		computes the partial partition function q*;
	// Called by slaveMutationRotamerSearch(.)
	private void slaveRotamerSearchAStar(boolean ligPresent, int numInAS, int numTotalRotamers,
			int rotamerIndexOffset[], int curResToASMap[], int residueMap[], boolean minimizeBB,
			boolean saveConfs, String fName){

		/*///////////////////////////////////////////////////////////////////////
		PrintStream logPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream("/net/grad/shaqbuzz/i+"+numRotSamples+".txt");
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
			numRotSamples++;
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		//for (int i=0;i<9;i++){logPS.println(curAANum[residueMap[i]]);}logPS.flush();
		///////////////////////////////////////////////////////////////////////*/
		
		ExpFunction ef = new ExpFunction();
		
		
		int treeLevels; // total num levels in the conformation tree
		if (ligPresent) //determine num tree levels: if ligPresent, then numInAS+1
			treeLevels = numInAS+1;
		else
			treeLevels = numInAS;
		
		int numRotForRes[] = new int[treeLevels]; //the number of rotamers for each flexible residue (AS+lig) during a mutation search
		int numRotForResNonPruned[] = new int[treeLevels]; //the number of non-pruned (by MinDEE) rotamers for each flexible residue
		int numTotalRotRed = 0;		//the total number of rotamers for the flexible residues only (during a mutation search)
		int numTotalRotRedNonPruned = 0; //the total num of non-pruned rotamers for the flexible residues
		int indicesEMatrix[] = null;//original (in the non-reduced matrices) indices of non-pruned rot to be included
		boolean eliminatedRotAtPosRed[] = null; //reduced MinDEE matrix
		float arpMatrixRed[][] = null; //reduced min energy matrix
		
		
		//Count the total number of conformations, the number of conformations pruned by MinDEE,
		//	the remaining conformations; the total num rotamers for the flexible residues, the num
		//	rotamers (total and non-pruned by MinDEE) for each fllexible residue;
		numTotalRotRed = countConfs(numInAS, treeLevels, numTotalRotamers, rotamerIndexOffset, ligPresent, 
				residueMap, numRotForRes, numRotForResNonPruned);
		
		BigInteger k_const = numConfsPrunedByMinDEE;
	
		BigInteger numConfsPrunedMinDEESteric = countPrunedByMinDEESteric(numInAS, treeLevels, numTotalRotamers, rotamerIndexOffset, ligPresent, 
				residueMap, numRotForRes, prunedIsSteric);
		
		k_const = k_const.subtract(numConfsPrunedMinDEESteric); //only the non-steric prunings are used in the computation of k_const
		
		//Bound the contribution of the conformations pruned by MinDEE
		BigDecimal pStar = ef.exp(-Ec_const/constRT).multiply(new BigDecimal(k_const));
		
		final double ro = (double)KSepsilon /(double)(1-KSepsilon);
		
		System.out.println("k_const: "+k_const+" pStar: "+pStar+" numConfsPrunedMinDEESteric: "+numConfsPrunedMinDEESteric);
		
		//Count the number of non-pruned rotamers
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if (numRotForResNonPruned[curLevel]==0){ //no non-pruned rotamers for curLevel, so no possible conformations
				allPruned = true;
				
				if (!k_const.equals(BigInteger.ZERO)){ //some non-sterics pruned but accuracy not achieved, so search must be repeated
					BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
					
					BigDecimal f = psi.divide(ef.exp(-Ec_const/constRT),4); //rounding is ROUND_HALF_UP
					BigInteger l_const = k_const.subtract( BigInteger.valueOf( (long)Math.ceil(f.doubleValue()) ) );
					setupRepeatRun(l_const, numRotForResNonPruned, numInAS); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
					repeatSearch = true;
				}
				
				return;
			}
			else
				numTotalRotRedNonPruned += numRotForResNonPruned[curLevel];
		}
		
		//logPS.println(Ec_const+" "+numConfsPrunedByMinDEE+" "+numConfsPrunedMinDEESteric+" "+k_const+" "+pStar+" "+numTotalRotRedNonPruned);
		
		indicesEMatrix = new int[numTotalRotRedNonPruned];
		eliminatedRotAtPosRed = new boolean[numTotalRotRed];
		arpMatrixRed = new float[numTotalRotRedNonPruned+1][numTotalRotRedNonPruned+1];//include the intra-energies in the last column
																		//and the shell-residue energies in the last row
		
		//Reduce the matrices
		reduceMatrices(eliminatedRotAtPosRed, arpMatrixRed,	indicesEMatrix, numRotForRes, numRotForResNonPruned, 
				treeLevels, numTotalRotRedNonPruned, numInAS, numTotalRotamers, rotamerIndexOffset, residueMap,
				ligPresent);
		
		//Set-up the A* search
		stericF = null;
		if (!minimizeBB) {//do not do a steric check if backbone minimization
			if (ligPresent)
				stericF = new StericCheck(curAANum,curResToASMap,residueMap,eliminatedRotAtPosRed,numRotForRes,
						m,overlapThresh,hSteric,numConfsLeft,numConfsAboveLevel,sysStrNum,sysLR,ligStrNum,ligROT,curLigNum,rl);
			else
				stericF = new StericCheck(curAANum,curResToASMap,residueMap,eliminatedRotAtPosRed,numRotForRes,
						m,overlapThresh,hSteric,numConfsLeft,numConfsAboveLevel,sysStrNum,sysLR,rl);
		}
		
		MSAStar AStarSearch = new MSAStar(treeLevels,numRotForResNonPruned,arpMatrixRed,stericF);
		
		curConf = new int[treeLevels]; //the rotamer sequence				
		boolean run1 = true;

		while (numConfsLeft.compareTo(BigInteger.ZERO)==1){

			curConf = AStarSearch.doAStar(run1,numInAS,null,eliminatedRotAtPosRed,sysLR,null,numRotForRes,residueMap,true,rl); //the current rotamer sequence); //the current rotamer sequence
			run1 = false;

			//Update the number of remaining conformations
			if (stericF!=null){
				numConfsLeft = stericF.getNumConfsLeft();
				numConfsPrunedByS = stericF.getNumConfsPrunedByS();
			}
			
			for (int curRotCheck=0; curRotCheck<treeLevels; curRotCheck++){//check if the conformation is valid
				if (curConf[curRotCheck]==-1){ // no valid conformations remaining
					if (partial_q.multiply(new BigDecimal(ro)).compareTo(pStar)<0){ //approximation accuracy not achieved
					
						if (!k_const.equals(BigInteger.ZERO)){ //some non-sterics pruned but accuracy not achieved, so search must be repeated
							BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
							
							BigDecimal f = psi.divide(ef.exp(-Ec_const/constRT),4); //rounding is ROUND_HALF_UP
							BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
							setupRepeatRun(l_const, numRotForResNonPruned, numInAS); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
							repeatSearch = true;
						}
					}
					
					return;
				}
			}
			
			//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
			//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
			//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
			//	of the rot to be applied for res 3 is 5)
			int conf[] = new int[curConf.length]; //the conformation with the actual rotamer numbers
			conf = getActualConf(curConf,eliminatedRotAtPosRed,treeLevels,numRotForRes,conf);
			
			//Apply the rotamers of the current conformation
			int curAS = 0;	
			for (int curLevel=0; curLevel<m.strand[sysStrNum].numberOfResidues; curLevel++){
				if (curResToASMap[curLevel]!=-1){//make a change only to the AS residues: use the native type for the other residues
										
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0){//not GLY or ALA
						sysLR.applyRotamer(m, curLevel, conf[curAS]);
						curASRotNum[curResToASMap[curLevel]] = conf[curAS];
					}
					else { //GLY or ALA
						curASRotNum[curResToASMap[curLevel]] = 0;
					}
					curAS++; //prepare the next AS residue
				}
			}		
			if (ligPresent){ //apply the ligand rotamer
				if (rl.getNumRotForAAtype(curLigNum)!=0){//not GLY or ALA
					ligROT.applyRotamer(m, 0, conf[treeLevels-1]);//the ligand level
					curLigRotNum = conf[treeLevels-1];
				}
				else { //GLY or ALA
					curLigRotNum = 0;
				}
			}
			
			for(int i=0;i<treeLevels;i++)System.out.print(conf[i]+" ");System.out.println();
			
			//Check the energy of the conformation and compute the score if necessary
			double minELowerBound = (double)(computeBestRotEnergyBound(numTotalRotamers,rotamerIndexOffset,ligPresent));
			//double psi = Math.max(initial_q,partial_q);
			BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
			
			//double curThreshold = -constRT * (Math.log(psi)+Math.log(ro/numConfsLeft));
			double curThreshold;
			BigDecimal diff_qp = psi.subtract(pStar);
			if (diff_qp.compareTo(new BigDecimal(0.0))<0) //the contribution of the pruned confs is bigger than ro*partial_q, so the search cannot be halted
				curThreshold = stericE;
			else
				curThreshold = -constRT * (ef.log(diff_qp).doubleValue()-Math.log(numConfsLeft.doubleValue()));

			System.out.println("conf: "+numConfsEvaluated.add(BigInteger.ONE)+" minELowerBound: "+minELowerBound+" curThreshold: "+curThreshold);
			System.out.println("pStar(double): "+pStar.doubleValue()+" qStar(double): "+partial_q.doubleValue()+" rho*qStar(double): "+ro*partial_q.doubleValue());
			
			if (minELowerBound > curThreshold) 
				return;
			
			else { //the energy of the new conformation is still below the threshold
				
				// After minimization do a m.updateCoordinates() to
				//  resync the actualCoordinates which were changed
				//  in the minimization procedure
				numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);
				numConfsLeft = numConfsLeft.subtract(BigInteger.ONE);
				if (stericF!=null)
					stericF.setNumConfsLeft(numConfsLeft); //update the number of conformations still to examine
				curVolume = 0.0f;
				float energy = 0.0f;
				if (doMinimization){
					if (!minimizeBB){ //side-chain minimization
						energy = calcTotalSnapshotEnergy();
						if (energy < bestEUnMin) {
							bestEUnMin = energy;
						}						
						simpMin.minimize(numMinSteps,false);
						energy = calcTotalSnapshotEnergy();
						if (doDihedE) //add dihedral energies
							energy += simpMin.computeDihedEnergy();
					}
					else { //backbone minimization
						energy = calcTotalSnapshotEnergy();
						if (energy < bestEUnMin) {
							bestEUnMin = energy;
						}						
						bbMin.minimizeFull(false);
						energy = calcTotalSnapshotEnergy();
					}
					
					if (saveConfs)
						saveMolecule(m,(fName+numConfsEvaluated),energy);
					
					m.updateCoordinates();
					if (energy < bestEMin) {
						bestEMin = energy;
					}
					if (energy < bestEnergyFound) {
						bestEnergyFound = energy;
					}
				}
				else if (computeEVEnergy){
					energy = calcTotalSnapshotEnergy();
					m.updateCoordinates();
					if (energy < bestEnergyFound) {
						bestEnergyFound = energy;
					}
					if (energy < bestEUnMin) {
						bestEUnMin = energy;
					}
				}
				partial_q = partial_q.add(ef.exp(-((double)(energy)) / constRT));
				
				System.out.println("energy: "+energy);
			}
		}
		if ((numConfsLeft.equals(BigInteger.ZERO))&&(!k_const.equals(BigInteger.ZERO))){ //no conformations remaining, non-sterics pruned			
			if (partial_q.multiply(new BigDecimal(ro)).compareTo(pStar)<0){ //approximation accuracy not achieved, so repeat search
				BigDecimal psi = initial_q.max(partial_q.multiply(new BigDecimal(ro)));
				
				BigDecimal f = psi.divide(ef.exp(-Ec_const/constRT),4); //rounding is ROUND_HALF_UP				
				BigInteger l_const = k_const.subtract(BigInteger.valueOf((long)Math.ceil(f.doubleValue())));
				setupRepeatRun(l_const, numRotForResNonPruned, numInAS); //accuracy not achieved, so repeat the search with reduced num pruned conf MinDEE
				repeatSearch = true;
			}
		}
	}
	
	//Reduce the original MinDEE and min energy matrices to contain only the entries for the
	//	rotamers of the current mutation sequence: the min energy matrix contains only the non-pruned entries,
	//	whereas the MinDEE matrix contains both the pruned and non-pruned entries for the current mutation sequence;
	//The matrices eliminatedRotAtPosRed[], arpMatrixRed[][], indicesEMatrix[] are changed and returned to 
	//	slaveRotamerSearchAStar()
	private void reduceMatrices(boolean eliminatedRotAtPosRed[], float arpMatrixRed[][], int indicesEMatrix[], 
			int numRotForRes[], int numRotForResNonPruned[], int treeLevels, 
			int numTotalRotRedNonPruned, int numInAS, int numTotalRotamers,
			int rotamerIndexOffset[], int residueMap[], boolean ligPresent){
		
		//Do the mapping from the original MinDEE and min energy matrices to the reduced ones
		int curIndexRed = 0;//index into the reduced matrices
		int pruningIndex = 0;//index into the reduced MinDEE matrix
		int curIndex;//index into the original minDEE matrix (curIndex+1 is in the original min energy matrix)
		
		//The entries for the ligand rotamers are after the the ones for the AS rotamers in the matrices
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				if ((ligPresent)&&(curLevel==(treeLevels-1))) //the ligand level
					curIndex = numInAS*numTotalRotamers + curRot;						
				else
					curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[curAANum[residueMap[curLevel]]] + curRot;
				
				if (!eliminatedRotAtRes[curIndex]){ //not pruned, so add its index
					indicesEMatrix[curIndexRed] = 1+curIndex;//add 1 for the shell/intra row/column in the original matrix
					curIndexRed++;
				}
				
				eliminatedRotAtPosRed[pruningIndex] = eliminatedRotAtRes[curIndex];
				//logPS.println(pruningIndex+" "+curIndex+" "+eliminatedRotAtRes[curIndex]);logPS.flush();
				pruningIndex++;
			}
		}

		//Reduce the min energy matrix
		for (int curRot1=0; curRot1<numTotalRotRedNonPruned; curRot1++){
			for (int curRot2=0; curRot2<numTotalRotRedNonPruned; curRot2++){
				arpMatrixRed[curRot1][curRot2] = arpMatrix[indicesEMatrix[curRot1]][indicesEMatrix[curRot2]];//pairwise
			}
			arpMatrixRed[curRot1][numTotalRotRedNonPruned] = arpMatrix[indicesEMatrix[curRot1]][0];//store intra-energies in the last column
			arpMatrixRed[numTotalRotRedNonPruned][curRot1] = arpMatrix[0][indicesEMatrix[curRot1]];//store shell-residue E in the last row
		}
	}
	
	//Get the actual numbers for the rotamers of a conformation that is returned by A* by including the information
	//		for the pruned rotamers
	private int [] getActualConf(int curConf[], boolean eliminatedRotAtPosRed[], 
			int treeLevels, int numRotForRes[], int conf[]){
		
		int curPruningInd = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			int curRotInd = 0;
			for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){
				if (!eliminatedRotAtPosRed[curPruningInd]){
					if (curRotInd==curConf[curLevel])
						conf[curLevel] = curRot;
					curRotInd++;
				}
				curPruningInd++;
			}
		}
		return conf;
	}
	
//// END ROTAMER SEARCH SECTION
/////////////////////////////////////////////////////////////////////////////

	
/*
 * 
 * DEE section
 * 
 */		
	//Compute the two interval terms in the summation of the MinDEE criteria
	public void doCompMinDEEIntervals(int numResInActiveSite, int numTotalRotamers,
			int numLigRotamers, int residueMap[], int rotIndexOffset[], boolean prunedRotAtRes[],
			boolean scaleInt, float maxIntScale){
		
		System.out.print("Computing MinDEE interval terms..");	
		
		float dist[][] = null;
		if (scaleInt)
			dist = doCompDistSC(numResInActiveSite,residueMap);
		
		MinDEEIntervals compInt = new MinDEEIntervals(arpMatrix, arpMatrixMax, numResInActiveSite, numTotalRotamers,
				numLigRotamers, rotIndexOffset,	residueMap, sysLR, prunedRotAtRes, scaleInt, dist, maxIntScale, rl);
		
		compInt.compMinDEEIntervals();
		indIntMinDEE = compInt.getIndInt();
		pairIntMinDEE = compInt.getPairInt();
		
		System.out.println("done.");
		
		if (debug){
			System.out.println();
			System.out.print(" ind: ");
			for (int curPos=0; curPos<indIntMinDEE.length; curPos++){
				System.out.print(indIntMinDEE[curPos]+" ");
			}
			System.out.println();
			System.out.print(" pair: ");
			for (int curPos=0; curPos<pairIntMinDEE.length; curPos++){
				System.out.print(pairIntMinDEE[curPos]+" ");
			}
			System.out.println();
		}
	}
	
	//Compute the distance between the side-chains for each active site residue pair (the ligand is not considered here);
	//	Returns the minimum distance between a pair of atoms in the two side-chains, for each side-chain pair
	private float[][] doCompDistSC(int numResInActiveSite, int residueMap[]){
		
		float dist[][] = new float[numResInActiveSite][numResInActiveSite];
		
		for (int i=0; i<numResInActiveSite; i++){
			Residue r1 = m.strand[sysStrNum].residue[residueMap[i]];
			for (int j=i+1; j<numResInActiveSite; j++){
				Residue r2 = m.strand[sysStrNum].residue[residueMap[j]];
				dist[i][j] = r1.getDist(r2,false);
				dist[j][i] = dist[i][j];
			}
		}
		return dist;
	}
	
	//Prune all rotamers that are incompatible with the template (intra E + res-to-template E >= stericE)
	public boolean[] DoPruneStericTemplate(int numInAS, int totalNumRotamers, int numLigRotamers, int residueMap[],
			int rotamerIndexOffset[], boolean prunedRotAtRes[], double stericE){
		
		eliminatedRotAtRes = prunedRotAtRes;
		
		if (prunedIsSteric==null){ //no pruning runs done yet
			prunedIsSteric = new boolean[prunedRotAtRes.length];
			for (int i=0; i<prunedIsSteric.length; i++)
				prunedIsSteric[i] = false;
		}
		
		int numAAtypes[] = new int[numInAS];
		for (int i=0; i<numAAtypes.length; i++) //the number of AAs allowed for each AS residue
			numAAtypes[i] = sysLR.getNumAllowable(residueMap[i]);
		
		int numPruned = 0;
		
		//Compute for the AS residues first
		for (int curPos=0; curPos<numInAS; curPos++){
			
			for (int AA=0; AA<numAAtypes[curPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[curPos],AA);
				
				//find how many rotamers are allowed for the current AA type at the given residue;
				//note that ala and gly have 0 possible rotamers
				int numRotForCurAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForCurAAatPos==0)	//ala or gly
					numRotForCurAAatPos = 1;
				
				for(int curRot=0; curRot<numRotForCurAAatPos; curRot++){
					
					int index_r = 1 + curPos*totalNumRotamers + rotamerIndexOffset[curAA] + curRot;
					
					if ((arpMatrix[index_r][0]+arpMatrix[0][index_r])>=stericE){
						
						eliminatedRotAtRes[index_r-1] = true;
						
						prunedIsSteric[index_r-1] = true;
						numPruned++;
					}
				}
			}
		}
		
		//If there is a ligand, compute for the lig rotamers as well
		if (numLigRotamers!=0){
			for (int curRot=0; curRot<numLigRotamers; curRot++){
				
				int index_r = 1 + numInAS*totalNumRotamers+curRot;
				
				if ((arpMatrix[index_r][0]+arpMatrix[0][index_r])>=stericE){
					
					eliminatedRotAtRes[index_r-1] = true;
					
					prunedIsSteric[index_r-1] = true;
					numPruned++;
				}
			}
		}
		
		System.out.println("Number of rotamers pruned due to incompatibility with the template: "+numPruned);		
		
		return eliminatedRotAtRes;
	}
	
	//Marks all rotamer pairs for which the min energy matrix entry is greater than cutoff as having a steric clash;
	//		If the max energy matrix exists, the corresponding entries are marked in the same way;
	//Note that row 0 and column 0 in the energy matrix are the intra- and res-to-template energies, so we must
	//		start checking at index 1 for both rows and columns, since we are interested in the pairwise energies
	public void preprocessPairs(float cutoff){
		
		if (arpMatrix==null)
			return;
		else {
			for (int i=1; i<arpMatrix.length; i++){
				for (int j=1; j<arpMatrix.length; j++){
					
					if (arpMatrix[i][j]>cutoff){
						
						arpMatrix[i][j] = stericE;
						if (arpMatrixMax!=null)
							arpMatrixMax[i][j] = stericE;
					}
				}
			}
		}
	}
	
	//Do Bounds Pruning
	public boolean[] DoMinBounds(int numInAS, int totalNumRotamers, int numLigRotamers, int residueMap[],
			int rotamerIndexOffset[], double pruningE, boolean prunedRotAtRes[], float initEw, boolean useSF, boolean boundKS){	
		
		MSMinBounds minBoundsRun = new MSMinBounds(arpMatrix,numInAS,totalNumRotamers,
				numLigRotamers,rotamerIndexOffset,residueMap,sysLR,
				pruningE,prunedRotAtRes,splitFlags,useSF,initEw,boundKS,rl);
		
		if (!boundKS){ //use Bounds to prune new rotamers
			eliminatedRotAtRes = minBoundsRun.ComputeEliminatedRotConf();
		}
		else { //compute Ec
			minBoundsRun.ComputeEliminatedRotConf();
			prunedIsSteric = minBoundsRun.getPrunedSteric();
			Ec_const = minBoundsRun.getEc();
		}
		
		minBoundsRun = null;
		
		return eliminatedRotAtRes;
	}
	
	//Do Bounding Flags
	public void DoBoundFlags(int numInAS, int totalNumRotamers, int numLigRotamers, int residueMap[],
			int rotamerIndexOffset[], double pruningE, boolean prunedRotAtRes[], float initEw, boolean useSF){	
		
		BoundFlags bFlagsRun = new BoundFlags(arpMatrix,numInAS,totalNumRotamers,
				numLigRotamers,rotamerIndexOffset,residueMap,sysLR,
				pruningE,prunedRotAtRes,splitFlags,useSF,initEw,rl);
		
		splitFlags = bFlagsRun.ComputeEliminatedRotConf();
		
		bFlagsRun = null;
	}
	
	//Do simple Goldstein DEE
	public boolean [] DoDEEGoldstein(int numResInActiveSite, int numTotalRotamers,
			int numLigRotamers, int residueMap[], int rotIndexOffset[], float initEw, boolean prunedRotAtRes[],
			boolean doMinimize, boolean useSF, boolean minimizeBB){
		
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		DEEGoldstein DEERun = new DEEGoldstein(arpMatrix, arpMatrixMax, numResInActiveSite, numTotalRotamers,
				numLigRotamers, rotIndexOffset,	residueMap, initEw, sysLR, prunedRotAtRes, doMinimize, 
				indIntMinDEE, pairIntMinDEE, splitFlags, useSF, minimizeBB, rl);
		
		eliminatedRotAtRes = DEERun.ComputeEliminatedRotConf();
		
		DEERun = null;
		
		return eliminatedRotAtRes;
	}
	
	//SplitDEE (conformational splitting) with 1 or 2 split positions
	public boolean [] DoDEEConfSplitting(int numResInActiveSite, int numTotalRotamers, int numLigRotamers, 
			int residueMap[], int rotIndexOffset[], float initEw, boolean prunedRotAtRes[], boolean resInMut[],
			int majorSplitPos[], int curDepth, boolean doMinimize, boolean useSF, int numSplits, 
			boolean distrDEE, boolean minimizeBB){
		
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		if (numSplits==1){ //1 split position
			DEESplit1f DEERunConfSplitting = new DEESplit1f(arpMatrix, arpMatrixMax, 
					numResInActiveSite, numTotalRotamers, numLigRotamers, rotIndexOffset, residueMap, initEw, 
					sysLR, prunedRotAtRes, resInMut, majorSplitPos, curDepth, doMinimize, indIntMinDEE, pairIntMinDEE, 
					splitFlags, useSF, distrDEE, minimizeBB, rl);
			
			eliminatedRotAtRes = DEERunConfSplitting.ComputeEliminatedRotConf();
			splitFlags = DEERunConfSplitting.getSplitFlags();
			
			DEERunConfSplitting = null;
		}
		else{ //2 split positions
			DEESplit2f DEERunConfSplitting = new DEESplit2f(arpMatrix, arpMatrixMax, 
					numResInActiveSite, numTotalRotamers, numLigRotamers, rotIndexOffset, residueMap, initEw, 
					sysLR, prunedRotAtRes, resInMut, majorSplitPos, curDepth, doMinimize, indIntMinDEE, pairIntMinDEE, 
					splitFlags, useSF, distrDEE, minimizeBB, rl);
			
			eliminatedRotAtRes = DEERunConfSplitting.ComputeEliminatedRotConf();
			splitFlags = DEERunConfSplitting.getSplitFlags();
			
			DEERunConfSplitting = null;
		}
		
		return eliminatedRotAtRes;
	}
	
	//Do Goldstein DEE pairs
	public void DoDEEPairs(int numResInActiveSite, int numTotalRotamers, int numLigRotamers, int residueMap[], 
			int rotIndexOffset[], float initEw, boolean prunedRotAtRes[], boolean resInPair[], boolean doMinimize, 
			boolean useSF, boolean magicBullet, boolean distrDEE, boolean minimizeBB, boolean scaleInt, float maxScale){
		
		//arpmatrix has a row/column for the backbone energies, so we just need
		//the # remaining entries, which is length-1 (only the AS and ligand rotamers)
		DEEGoldsteinPairs DEERunPairs = new DEEGoldsteinPairs(arpMatrix, arpMatrixMax, numResInActiveSite, numTotalRotamers,
				numLigRotamers, rotIndexOffset,	residueMap, initEw, sysLR, prunedRotAtRes, resInPair, doMinimize, 
				splitFlags, useSF, magicBullet, distrDEE, minimizeBB, scaleInt, maxScale, rl);
		
		DEERunPairs.ComputeEliminatedRotConf();		
		splitFlags = DEERunPairs.getSplitFlags();
		
		DEERunPairs = null;
	}
	/*
	 * 
	 * End of DEE section
	 * 
	 */	
	
//////////////////////////////////////////////////////////////////////////
//	Compute Min GMEC section
//////////////////////////////////////////////////////////////////////////
	public void doAStarGMEC(String fileName, boolean searchComputeEVEnergy, 
			boolean searchDoMinimization,int numInAS, int numTotalRotamers, int rotamerIndexOffset[], 
			int residueMap[], String resDefault[],int numMut, float Ew, double bestScore, 
			CommucObj cObj, boolean approxMinGMEC, float lambda, boolean minimizeBB, boolean useEref, float eRef[]) {
	
		// A rotamer search is performed. For each residue,
		//  every allowable rotamer is tried in combination
		//  with every other allowable rotamer
		// If we're doing a mutation search then residues
		//  are allowed to mutate
	
		numConfsEvaluated = BigInteger.ZERO;
		computeEVEnergy = searchComputeEVEnergy;
		doMinimization = searchDoMinimization;
		ASAANums = new int[numInAS];
		curASRotNum = new int[numInAS];
		int curResToASMap[] = new int[m.strand[sysStrNum].numberOfResidues];
			// This map maps the system strand residues back to the AS numbering
			// So something like 8 -> 0, 10 -> 1, 11 -> 2, ...
		curLigRotNum = 0;
		
		for (int i=0; i<numInAS; i++) //the AS residues are flexible - this is used by simpMin to set up the minimizer
			m.residue[residueMap[i]].flexible = true;
		if (haveLigand) //the ligand is also flexible
			m.strand[ligStrNum].residue[0].flexible = true;
		
	
		if (searchDoMinimization && !searchComputeEVEnergy){
			System.out.println("Warning: In order to do minimization computeEVEnergy must be true");
			return;
		}
	
		if (haveLigand && !ligandReady){
			System.out.println("Warning: Ligand rotamers not setup, RotamerSearch aborting");
			return;
		}
	
		// Prepare Amber
		if(searchComputeEVEnergy){
			// Amber should already be loaded
			if(haveLigand) {
				a96ff.setLigandNum(m.strand[ligStrNum].residue[0].moleculeResidueNumber);
			}
			if (doMinimization){
				if (!minimizeBB){ //side-chain minimization
					if (simpMin == null) {
						System.out.println("Warning: Attempting minimization run but simpMin not allocated, RotamerSearch aborting");
						return;
					}
					bbMin = null;
				}
				else { //backbone minimization
					if (bbMin == null) {
						System.out.println("Warning: Attempting minimization run but bbMin not allocated, RotamerSearch aborting");
						return;
					}
					simpMin = null;
				}
			}
		}
	
		// Make sure the allRotamerPairsEnergyName matrices exist
		if (arpMatrix == null) {
			System.out.println("Warning: allRotamerPairsEnergy matrix not loaded");
			return;
		}
	
		// Setup the residue number to AS number map
		for(int i=0;i<m.strand[sysStrNum].numberOfResidues;i++){
			curResToASMap[i] = -1;
		}
		for(int i=0;i<residueMap.length;i++){
			curResToASMap[residueMap[i]] = i;
		}
		
		doAStarGMECHelper(haveLigand, numInAS, numTotalRotamers, rotamerIndexOffset, curResToASMap,
				residueMap, fileName, resDefault, numMut, Ew, bestScore, cObj, approxMinGMEC, lambda, minimizeBB, useEref, eRef);
	}

	// Calls AStar repeatedly while the returned conformations still have energy below the threshold;
	//		computes the partial partition function q*;
	// Called by mutationRotamerSearch(.)
	private void doAStarGMECHelper(boolean ligPresent, int numInAS, int numTotalRotamers,
			int rotamerIndexOffset[], int curResToASMap[], int residueMap[], String fileName, 
			String resDefault[], int numMut, float Ew, double bestScore, CommucObj cObj, 
			boolean approxMinGMEC, float lambda, boolean minimizeBB, boolean useEref, float eRef[]){
		
		/*/////////////////////////////////////////////////////////////////////
		PrintStream debugPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream("debug.txt");
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			debugPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}
		///////////////////////////////////////////////////////////////////////*/
		
		boolean outputFile = (fileName!=null); //output to file
	
		int treeLevels; // total num levels in the conformation tree
		if (ligPresent) //determine num tree levels: if ligPresent, then numInAS+1
			treeLevels = numInAS+1;
		else
			treeLevels = numInAS;
		
		int numRotForRes[] = new int[treeLevels]; //the number of rotamers for each flexible residue (AS+lig) during a mutation search
		int numRotForResNonPruned[] = new int[treeLevels]; //the number of non-pruned (by MinDEE) rotamers for each flexible residue
		int numTotalRotRed = 0;		//the total number of rotamers for the flexible residues only (during a mutation search)
		int numTotalRotRedNonPruned = 0; //the total num of non-pruned rotamers for the flexible residues
		int indicesEMatrix[] = null;//original (in the non-reduced matrices) indices of non-pruned rot to be included
		boolean eliminatedRotAtPosRed[] = null; //reduced MinDEE matrix
		float arpMatrixRed[][] = null; //reduced min energy matrix
		
		int numTotalConf = 1;
		int curNumRot = 0;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				curNumRot = rl.getNumRotForAAtype(ligROT.getIndexOfNthAllowable(0,0));
				if (curNumRot==0) //GLY or ALA
					curNumRot = 1;
			}
			else { //AS residue				
				curNumRot = 0;
				for (int i=0; i<sysLR.getNumAllowable(residueMap[curLevel]); i++){ //add the rot for all allowable AA at this residue
					int newRot = rl.getNumRotForAAtype(sysLR.getIndexOfNthAllowable(residueMap[curLevel],i));
					if (newRot==0) //GLY or ALA
						newRot = 1;
					curNumRot += newRot; 
				}
			}
			numRotForRes[curLevel] = curNumRot;
			numRotForResNonPruned[curLevel] = numRotForRes[curLevel];
			numTotalRotRed += numRotForRes[curLevel];
			
			numTotalConf *= numRotForRes[curLevel];
		}
		
		int numPrunedThisLevel;
		
		for (int curLevel=0; curLevel<treeLevels; curLevel++){ //for each residue
			int curIndex;		System.out.println("curLevel "+curLevel);
			numPrunedThisLevel = 0;
			if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){ //for all rotamers for the given AA
					curIndex = numInAS*numTotalRotamers + curRot;
					if (eliminatedRotAtRes[curIndex]){
						numPrunedThisLevel++;
					}
				}
			}
			else { //AS residue
				for (int curAA=0; curAA<sysLR.getNumAllowable(residueMap[curLevel]); curAA++){ //for all allowed AA's
					int index = sysLR.getIndexOfNthAllowable(residueMap[curLevel],curAA);
					int newRot = rl.getNumRotForAAtype(index);
					if (newRot==0)
						newRot = 1;			System.out.print(newRot+" ");
					for (int curRot=0; curRot<newRot; curRot++){ //for all rotamers for the given AA
						curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[index]+curRot;
						System.out.print(curIndex+" "+eliminatedRotAtRes[curIndex]+" ");
						if (eliminatedRotAtRes[curIndex]){
							numPrunedThisLevel++;
						}
					}System.out.println();
				}
			}
			numRotForResNonPruned[curLevel] -= numPrunedThisLevel;
		}
		
		int numConfNonPrunedMinDEE = 1;
		for (int curLevel=0; curLevel<treeLevels; curLevel++){
			if (numRotForResNonPruned[curLevel]==0){ //no non-pruned rotamers for curLevel, so no possible conformations
				System.out.println("MinDEE has pruned all possible conformations: position "+curLevel); //at least the minGMEC should be remaining
				if (cObj!=null) {//output num conf info to cObj
					cObj.EL_searchNumConfsTotal = numTotalConf;
					cObj.EL_searchNumPrunedMinDEE = numTotalConf;
					cObj.EL_searchNumConfsEvaluated = 0; //no confs evaluated
					cObj.bestBoundEMin = stericE;
				}
				
				return;
			}
			else {
				numTotalRotRedNonPruned += numRotForResNonPruned[curLevel];
				
				numConfNonPrunedMinDEE *= numRotForResNonPruned[curLevel];
			}
		}
		
		if (cObj!=null) {//output num conf info to cObj
			cObj.EL_searchNumConfsTotal = numTotalConf;
			cObj.EL_searchNumPrunedMinDEE = (numTotalConf - numConfNonPrunedMinDEE);
		}
		
		for(int i=0;i<treeLevels;i++)System.out.print(numRotForRes[i]+" ");System.out.println();
		System.out.println(numTotalRotRed);
		for(int i=0;i<treeLevels;i++)System.out.print(numRotForResNonPruned[i]+" ");System.out.println(numTotalRotRedNonPruned);
		
		indicesEMatrix = new int[numTotalRotRedNonPruned];
		eliminatedRotAtPosRed = new boolean[numTotalRotRed];
		arpMatrixRed = new float[numTotalRotRedNonPruned+1][numTotalRotRedNonPruned+1];//include the intra-energies in the last column
																		//and the shell-residue energies in the last row
		
		int curIndexRed = 0;//index into the reduced matrices
		int pruningIndex = 0;//index into the reduced MinDEE matrix
		int curIndex;//index into the original minDEE matrix (curIndex+1 is in the original min energy matrix)
		
		//The entries for the ligand rotamers are after the the ones for the AS rotamers in the matrices
		for (int curLevel=0; curLevel<treeLevels; curLevel++){			
			if ((ligPresent)&&(curLevel==(treeLevels-1))){ //the ligand level
				for (int curRot=0; curRot<numRotForRes[curLevel]; curRot++){ //for all rotamers for the given AA
					curIndex = numInAS*numTotalRotamers + curRot;
					if (!eliminatedRotAtRes[curIndex]){ //not pruned, so add its index
						indicesEMatrix[curIndexRed] = 1+curIndex;//add 1 for the shell/intra row/column in the original matrix
						curIndexRed++;
					}
					eliminatedRotAtPosRed[pruningIndex] = eliminatedRotAtRes[curIndex];
					//logPS.println(pruningIndex+" "+curIndex+" "+eliminatedRotAtRes[curIndex]);logPS.flush();
					pruningIndex++;
				}
			}
			else { //AS residue
				for (int curAA=0; curAA<sysLR.getNumAllowable(residueMap[curLevel]); curAA++){ //for all allowed AA's
					int index = sysLR.getIndexOfNthAllowable(residueMap[curLevel],curAA);
					int newRot = rl.getNumRotForAAtype(index);
					if (newRot==0)
						newRot = 1;
					for (int curRot=0; curRot<newRot; curRot++){ //for all rotamers for the given AA
						curIndex = curLevel*numTotalRotamers + rotamerIndexOffset[index]+curRot;
						
						if (!eliminatedRotAtRes[curIndex]){ //not pruned, so add its index
							indicesEMatrix[curIndexRed] = 1+curIndex;//add 1 for the shell/intra row/column in the original matrix
							curIndexRed++;
						}
						eliminatedRotAtPosRed[pruningIndex] = eliminatedRotAtRes[curIndex];
						//logPS.println(pruningIndex+" "+curIndex+" "+eliminatedRotAtRes[curIndex]);logPS.flush();
						pruningIndex++;
					}
				}
			}
		}
		int AAdefault[] = new int[resDefault.length];
		for (int i=0; i<AAdefault.length; i++){
			AAdefault[i] = rl.getAARotamerIndex(resDefault[i]);
		}
		
		
		//for (int i=0; i<AAdefault.length; i++)debugPS.print(AAdefault[i]+" ");debugPS.println();debugPS.flush();	
		
		for (int i=0; i<AAdefault.length; i++)System.out.print(AAdefault[i]+" ");System.out.println();
		//for(int i=0;i<rotIndexOffsetRed.length;i++)System.out.print(rotIndexOffsetRed[i]+" ");System.out.println();
		//for(int i=0;i<AAindMap.length;i++)System.out.print(AAindMap[i]+" ");System.out.println();
		System.out.println("pruneIndex "+pruningIndex);
		for (int i=0;i<curIndexRed;i++)System.out.print(indicesEMatrix[i]+" ");System.out.println("curIndexRed "+curIndexRed);
	
		//Reduce the min energy matrix
		for (int curRot1=0; curRot1<numTotalRotRedNonPruned; curRot1++){
			for (int curRot2=0; curRot2<numTotalRotRedNonPruned; curRot2++){
				arpMatrixRed[curRot1][curRot2] = arpMatrix[indicesEMatrix[curRot1]][indicesEMatrix[curRot2]];//pairwise
			}
			arpMatrixRed[curRot1][numTotalRotRedNonPruned] = arpMatrix[indicesEMatrix[curRot1]][0];//store intra-energies in the last column
			arpMatrixRed[numTotalRotRedNonPruned][curRot1] = arpMatrix[0][indicesEMatrix[curRot1]];//store shell-residue E in the last row
		}
		
		//Declaring the logPS output here prevents opening an empty file and returning
		//	(for example, if all conformations are pruned by MinDEE above)
		PrintStream logPS = null; //the output file for conf info
		if (outputFile){
			try {			
				FileOutputStream fileOutputStream = new FileOutputStream(fileName,true); //append file if more than 1 partition
				BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
				logPS = new PrintStream( bufferedOutputStream );
			}
			catch (Exception ex) {
				System.out.println("ERROR: An exception occured while opening log file");
			}
		}
		
		
		//Set-up the A* search					
		MSAStar MSAStarSearch = new MSAStar(treeLevels,numRotForResNonPruned,arpMatrixRed,null);
		
		curConf = new int[treeLevels]; //the rotamer sequence				
		boolean run1 = true;
		int numConfsOutput = 0;//num confs output to file
		float lowestBound = stericE;
		
		if (!doMinimization)
			approxMinGMEC = false; //reset approxMinGMEC, since it is valid only for MinDEE, and not for traditional DEE

		float curMinE = stericE;
		curMinE = (float)Math.min(curMinE,bestScore); //the best energy for the already computed sequences
		bestEMin = (float)Math.min(bestEMin, bestScore);
	
		while (true){
			
			//debugPS.println("curMinE: "+curMinE);debugPS.flush();
			
			if (cObj!=null){				
				/*if (numConfsEvaluated>=cObj.confSeq.length){
					CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[2*cObj.confSeq.length];
					System.arraycopy(cObj.confSeq,0,tempConf,0,cObj.confSeq.length);
					cObj.confSeq = tempConf;
				}
				cObj.confSeq[numConfsEvaluated] = cObj.new ConfInfo(treeLevels);*/
				
				cObj.EL_searchNumConfsEvaluated = numConfsEvaluated.intValue();
			}
			
			//clear the values from the previous run
			for (int i=0; i<numInAS; i++){
				curAANum[i] = -1;
				curASRotNum[i] = -1;
			}
			curLigNum = -1;
			curLigRotNum = -1;
			
			
			curConf = MSAStarSearch.doAStar(run1,numMut,AAdefault,eliminatedRotAtPosRed,sysLR,resDefault,numRotForRes,residueMap,false,rl); //the current rotamer sequence
			
			System.out.println("confNum: "+(numConfsEvaluated.add(BigInteger.ONE)));
			System.out.print("curConf: ");for(int i=0;i<treeLevels;i++)System.out.print(curConf[i]+" ");System.out.println();
			
			//debugPS.println("confNum: "+(numConfsEvaluated+1));
			//debugPS.print("curConf: ");for(int i=0;i<treeLevels;i++)debugPS.print(curConf[i]+" ");debugPS.println();
			//debugPS.flush();
			
			
			
			for (int curRotCheck=0; curRotCheck<treeLevels; curRotCheck++){//check if the conformation is valid
				if (curConf[curRotCheck]==-1){ // no valid conformations remaining
					/*CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[numConfsEvaluated];
					System.arraycopy(cObj.confSeq,0,tempConf,0,numConfsEvaluated);
					cObj.confSeq = tempConf;*/
					
					if (cObj!=null)
						cObj.bestScore = new BigDecimal(curMinE); //update the best score so far to supply to the next partition
					
					MSAStarSearch = null;
					
					if (outputFile){
						logPS.flush(); //there may still be results to output
					}
					return;
				}
			}
			
			//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
			//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
			//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
			//	of the rot to be applied for res 3 is 5)
			int conf[] = new int[curConf.length]; //the conformation with the actual rotamer numbers
			conf = getActualConf(curConf,eliminatedRotAtPosRed,treeLevels,numRotForRes,conf);
			
			System.out.print("actualConf: ");for(int i=0;i<treeLevels;i++)System.out.print(conf[i]+" ");System.out.println();
			
			//debugPS.print("actualConf: ");for(int i=0;i<treeLevels;i++)debugPS.print(conf[i]+" ");debugPS.println();debugPS.flush();
			
			
			//Extract the AA numbers for the current conformation and appply the corresponding AA types
			for (int i=0; i<treeLevels; i++){
				if ((ligPresent)&&(i==(treeLevels-1))){ //the ligand level: only one type is allowed
					curLigNum = ligROT.getIndexOfNthAllowable(0,0);
					ligROT.changeResidueType(m,0,rl.getAAName(curLigNum),addHydrogens);
				}
				else { // AS residue
					curAANum[residueMap[i]] = getAAIndex(conf[i],i,resDefault,residueMap);
					sysLR.changeResidueType(m,residueMap[i],rl.getAAName(curAANum[residueMap[i]]),addHydrogens,connectResidues);
				}
			}
			
			for (int i=0; i<numInAS; i++)
				ASAANums[i] = curAANum[residueMap[i]];
			
			System.out.print("curAANum: ");for(int i=0;i<numInAS;i++)System.out.print(ASAANums[i]+" ");System.out.println(curLigNum);
			
			
			//debugPS.print("curAANum: ");for(int i=0;i<numInAS;i++)debugPS.print(ASAANums[i]+" ");debugPS.println(curLigNum);debugPS.flush();
			
			
			//Extract and apply the rotamers of the current conformation
			int curAS = 0;	
			int curRot = 0;
			for (int curLevel=0; curLevel<m.strand[sysStrNum].numberOfResidues; curLevel++){
				if (curResToASMap[curLevel]!=-1){//make a change only to the AS residues: use the native type for the other residues
										
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0){//not GLY or ALA
						curRot = conf[curAS] - getRotSum(conf[curResToASMap[curLevel]],curResToASMap[curLevel],resDefault,residueMap);
						sysLR.applyRotamer(m, curLevel, curRot);
						curASRotNum[curResToASMap[curLevel]] = curRot;
					}
					else { //GLY or ALA
						curASRotNum[curResToASMap[curLevel]] = 0;
					}
					curAS++; //prepare the next AS residue
				}
			}		
			if (ligPresent){ //apply the ligand rotamer
				if (rl.getNumRotForAAtype(curLigNum)!=0){//not GLY or ALA
					ligROT.applyRotamer(m, 0, conf[treeLevels-1]);//the ligand level
					curLigRotNum = conf[treeLevels-1];
				}
				else { //GLY or ALA
					curLigRotNum = 0;
				}
			}	
			
			System.out.print("curRot: ");for(int i=0;i<numInAS;i++)System.out.print(curASRotNum[i]+" ");System.out.println(curLigRotNum);
			
			//debugPS.print("curRot: ");for(int i=0;i<numInAS;i++)debugPS.print(curASRotNum[i]+" ");debugPS.println(curLigRotNum);debugPS.flush();
			
			
			
			// After minimization do a m.updateCoordinates() to
			//  resync the actualCoordinates which were changed
			//  in the minimization procedure
			
			//////////////////////////////////////////////////////////////////////////////////////////
			// Energy computation
			float unMinE = 0.0f;
			float minE = 0.0f;
				
			float minELowerBound = computeBestRotEnergyBound(numTotalRotamers,rotamerIndexOffset,ligPresent);
			
			//debugPS.println("minELowerBound: "+minELowerBound);debugPS.flush();
			
			
			if (run1) //this is the first extracted conformation, so it has the lowest energy lower bound, so store it
				lowestBound = minELowerBound;
			
			
			boolean done = false;
			
			if ((minELowerBound>(curMinE+Ew)) && (!run1)){ //we already have all confs within Ew of the minGMEC
				done = true;
			}
			else if (approxMinGMEC){ //running the heuristic halting condition
				
				if ((minELowerBound>(lowestBound+lambda)) && (!run1)) //compare the current bound to the lowest bound
					done = true;
			}
			
			
			if (done){ //we already have all required conformations
				
				if (cObj!=null)
					cObj.bestScore = new BigDecimal(curMinE); //update the best score so far to supply to the next partition
				
				MSAStarSearch = null;
				
				/*CommucObj.ConfInfo tempConf[] = new CommucObj.ConfInfo[numConfsEvaluated];				
				System.arraycopy(cObj.confSeq,0,tempConf,0,numConfsEvaluated);				
				cObj.confSeq = tempConf;*/
				
				if (outputFile){
					logPS.flush();
				}
				return; //stop the search
			}
			
			else{
				if (computeEVEnergy){
					a96ff.calculateTypesWithTemplates();
					a96ff.initializeCalculation();
					a96ff.setNBEval(hElect,hVDW);
					if (doMinimization){
						if (!minimizeBB){
							if(ligPresent){
								simpMin.initialize(m,sysStrNum,ligStrNum,a96ff,sysLR,ligROT,curAANum,curLigNum,doDihedE,rl);
							}
							else
								simpMin.initialize(m,sysStrNum,a96ff,sysLR,curAANum,doDihedE,rl);
						}
						else { //backbone minimization
							if (ligPresent)
								bbMin.initialize(m, a96ff, residueMap, sysStrNum, ligStrNum);
							else
								bbMin.initialize(m, a96ff, residueMap, sysStrNum);
						}
					}
				}
				if (doMinimization){
					if (!minimizeBB) {//side-chain minimization
						unMinE = calcTotalSnapshotEnergy();
						simpMin.minimize(numMinSteps,false);
						minE = calcTotalSnapshotEnergy();
						if (doDihedE) //add dihedral energies
							minE += simpMin.computeDihedEnergy();
					}
					else {//backbone minimization
						unMinE = calcTotalSnapshotEnergy();
						bbMin.minimizeFull(false);
						minE = calcTotalSnapshotEnergy();
					}
					m.updateCoordinates();
				}
				else if (computeEVEnergy){ //no minimization, so traditional DEE
					minE = calcTotalSnapshotEnergy(); //the sum of the precomputed energy terms
					m.updateCoordinates();
				}
				
				float totEref = 0.0f;
				if (useEref)
					totEref = getTotSeqEref(eRef,numInAS,residueMap,ligPresent);
				unMinE -= totEref;
				minE -= totEref;
				
				bestEMin = Math.min(bestEMin,minE); //for updating from within KStar (getBestE(.))
				curMinE = Math.min(curMinE,minE); //for the halting condition
				////////////////////////////////////////////////////////////////////////////////////////////
				
				
				System.out.println(minELowerBound+" "+minE+" "+curMinE);				
				
				//Since we only need to save the information for conformations whose energy is within Ew of
				//	the minGMEC, we can compare the minimized energy for the current conformation to the
				//	lowest minimized energy in the search so far and output only the conformations that
				//	are within Ew. This optimization is important, since writing to the output file is
				//	relatively very exomputationally expensive; moreover, the output file can become very
				//	big, so the optimization reduces the storage requirement. This approach is most beneficial
				//	when low minimized energies are returned early in thse serach, so there will be only a small
				//	number of extra output conformations.
				if (outputFile){//Output to file
					
					if ((approxMinGMEC)||(minE<=(curMinE+Ew))){ //heuristic stopping condition or minE within Ew of the current lowest energy
						
						numConfsOutput++;
						logPS.print(numConfsOutput+" ");
						for (int i=0; i<treeLevels; i++){
							if ((ligPresent)&&(i==(treeLevels-1)))
								logPS.print(rl.getAAName(curLigNum)+" ");
							else
								logPS.print(rl.getAAName(curAANum[residueMap[i]])+" ");
						}
						for (int i=0; i<treeLevels; i++){
							if ((ligPresent)&&(i==(treeLevels-1)))
								logPS.print(curLigRotNum+" ");
							else
								logPS.print(curASRotNum[i]+" ");
						}
						//logPS.println();
						logPS.print("unMinE: "+unMinE+" ");
						logPS.print("minE: "+minE+" ");
						if (doMinimization)
							logPS.println("minBound: "+minELowerBound);
						//if (numConfsOutput%100==0) //flush only every 100 confs, since output is relatively *very* computationally expensive
							logPS.flush();
					}
				}
			}
			
			numConfsEvaluated = numConfsEvaluated.add(BigInteger.ONE);
			
			run1 = false;
		}
	}
	
	//Returns the reference energy for the current amino acid sequence assignment (for the mutatable positions only)
	private float getTotSeqEref(float eRef[], int numInAS, int residueMap[], boolean ligPresent){
		
		float totEref = 0.0f;
		for (int i=0; i<numInAS; i++)
			totEref += eRef[curAANum[residueMap[i]]];
		
		if (ligPresent)
			totEref += eRef[curLigNum];
		
		return totEref;
	}
	
	private int getAAIndex(int rotIndex, int curRes, String resDefault[], int residueMap[]){
		
		int rotSum = 0;
		for (int i=0; i<sysLR.getNumAllowable(residueMap[curRes]); i++){
			int curRot = rl.getNumRotamers(rl.getAAName(sysLR.getIndexOfNthAllowable(residueMap[curRes],i)));
			if (curRot==0) //GLY or ALA
				curRot = 1;
			rotSum += curRot;
			if (rotSum>rotIndex)
				return (sysLR.getIndexOfNthAllowable(residueMap[curRes],i));
		}
		return -1;
	}
	
	private int getRotSum(int rotIndex, int curRes, String resDefault[], int residueMap[]){
		
		int rotSum = 0;
		for (int i=0; i<sysLR.getNumAllowable(residueMap[curRes]); i++){
			int curRot = rl.getNumRotamers(rl.getAAName(sysLR.getIndexOfNthAllowable(residueMap[curRes],i)));
			if (curRot==0) //GLY or ALA
				curRot = 1;
			if ((rotSum+curRot)>rotIndex)
				return rotSum;
			else
				rotSum += curRot;
		}
		return -1;
	}
//////////////////////////////////////////////////////////////////////////
//	End of Compute Min GMEC section
//////////////////////////////////////////////////////////////////////////
	
	
//////////////////////////////////////////////////////////////////////////
//	Generate Backbones section
//////////////////////////////////////////////////////////////////////////
	public void doGenBackbones(String runName, int numInAS, int residueMap[], double theta, double alpha,
			int numSamples, boolean systematicSampling){
	
		//setup the log file to store information about the generated backbones
		PrintStream logPS = null;
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream((runName+".bb.all"));
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
			System.exit(1);
		}
		
		Backbone bb = new Backbone();
		
		//the initial phi and psi angles for the residues with flexible backbones
		double initFi[] = new double[numInAS];
		double initPsi[] = new double[numInAS];
		
		for (int i=0; i<numInAS; i++){
			initFi[i] = bb.getFiPsi(m, sysStrNum, residueMap[i], 0);
			initPsi[i] = bb.getFiPsi(m, sysStrNum, residueMap[i], 1);
			
			System.out.println("AS residue: "+i+" ("+m.strand[sysStrNum].residue[residueMap[i]].getResNumber()+") phi: "+initFi[i]+" psi: "+initPsi[i]);
			
			if ((initFi[i]==0.0)||(initPsi[i]==0.0)){
				System.out.println("ERROR: Residue "+residueMap[i]+" does not have a valid (phi,psi) pair.");
				System.exit(1);
			}
		}
		
		if (systematicSampling){ //systematic sampling (generates a large number of backbones)
			//the number of steps (from -theta to +theta)
			int numSteps = 1;
			if (alpha!=0.0)
				numSteps = 2*(int)(theta/alpha) + 1; //count the initial phi/psi values as well
			
			doGenBackbonesHelper(runName, numInAS, residueMap, theta, alpha, numSteps, bb, initFi, initPsi, logPS, 0);
		}
		else { //random sampling (generates a small number of backbones)
			outputBB(logPS, numInAS, runName, residueMap, bb); //output the initial backbone first
			doGenBackbonesRandHelper(runName, numInAS, residueMap, theta, numSamples, bb, logPS);
		}
	}
	
	//Generates up to numSamples backbones by applying random phi/psi changes within theta degrees 
	//		This version is used to generate a small number of backbones (the systematic version below
	//		will generate a very large number of backbones, even for a very small number of steps)
	private void doGenBackbonesRandHelper(String runName, int numInAS, int residueMap[], double theta,
			int numSamples, Backbone bb, PrintStream logPS){
		
		Random randNum = new Random();
		
		for (int curSample=0; curSample<numSamples; curSample++){ //generate up to numSamples backbones
			
			float curFiChange[] = new float[numInAS]; //the phi changes for each residue
			float curPsiChange[] = new float[numInAS]; //the psi changes for each residue
			
			for(int curRes=0; curRes<numInAS; curRes++){ //apply the random (phi,psi) changes for each residue
			
				//Get the random (phi,psi) change
				curFiChange[curRes] = 2*(randNum.nextFloat()-0.5f)*(float)theta; //phi change within (-theta,theta)
				curPsiChange[curRes] = 2*(randNum.nextFloat()-0.5f)*(float)theta; //psi change within (-theta,theta)
				
				//Apply the (phi,psi) change
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],curFiChange[curRes],0);
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],curPsiChange[curRes],1);
			}
			
			//we have a full backbone conformation, so output
			outputBB(logPS, numInAS, runName, residueMap, bb);
			
			for(int curRes=0; curRes<numInAS; curRes++){ //restore the original (phi,psi) for each residue
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],-curFiChange[curRes],0);
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],-curPsiChange[curRes],1);
			}
		}
	}
	
	//Recursively generates backbones by applying phi/psi changes within theta degrees 
	//		at a step of alpha degrees to each residue with a flexible backbone (systematic sampling)
	private void doGenBackbonesHelper(String runName, int numInAS, int residueMap[], double theta, double alpha, 
			int numSteps, Backbone bb, double initFi[], double initPsi[], PrintStream logPS, int curRes){
		
		if (curRes==numInAS){//we have a full backbone conformation, so output
			outputBB(logPS, numInAS, runName, residueMap, bb);
		}
		else { //we only have a partial backbone, so apply changes to the current residue
			
			if (curRes==0) 
				System.out.println("Starting..");
			
			//First, move the phi/psi angle to -theta, then apply the changes to +theta, at alpha steps
			
			//move phi to -(theta+alpha), so that the first move in the *for* statement below goes to -theta
			bb.applyFiPsi(m,sysStrNum,residueMap[curRes],-(theta+alpha),0);
			
			//apply phi changes up to +theta
			for (int curStepFi=0; curStepFi<numSteps; curStepFi++){ //apply each alpha step up to a displacement of +theta from the initial phi
				
				if (curRes==0)
					System.out.print("*");
				
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],alpha,0);
				
				//move psi to -(theta+alpha), so that the first move in the *for* statement below goes to -theta
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],-(theta+alpha),1);
				
				for (int curStepPsi=0; curStepPsi<numSteps; curStepPsi++){ //apply each alpha step up to a displacement of +theta from the initial psi
				
					if (curRes==0)
						System.out.print(".");
					
					bb.applyFiPsi(m,sysStrNum,residueMap[curRes],alpha,1);
					if (checkStericsBBonly(sysStrNum,residueMap[curRes])) {//passed steric test, so move to the next residue
						doGenBackbonesHelper(runName,numInAS,residueMap,theta,alpha,numSteps,bb,initFi,initPsi,logPS,curRes+1);
					}
				}
				
				//restore initial psi
				bb.applyFiPsi(m,sysStrNum,residueMap[curRes],-theta,1);
				
				if (curRes==0)
					System.out.println();
			}
			
			//restore initial phi
			bb.applyFiPsi(m,sysStrNum,residueMap[curRes],-theta,0);
			
			if (curRes==0)
				System.out.println("done");
		}
	}
	
	//Checks if the current backbone is sterically allowed and outputs the pdb and log information
	private void outputBB(PrintStream logPS, int numInAS, String runName, int residueMap[], Backbone bb){
		
		//Check all residues for sterics against the whole strand since all residues
		//		are already assigned (up to this point, we have checked for sterics only
		//		against the residues up to a given AS residue);
		//The backbone movement may have actually caused unallowed sterics between rigid
		//		parts of the molecule, so we check for this possibility
		boolean stericAllowed = true;
		for (int i=0; i<m.strand[sysStrNum].numberOfResidues; i++){
			if (!checkAllStericsBBonly(sysStrNum,i)){
				stericAllowed = false;
				break;
			}
		}
		
		if (stericAllowed){ //all sterics are allowed
		
			String fileName = (runName+System.currentTimeMillis()+".pdb");
			
			saveMolecule(m,fileName,0.0f);//save the molecule
			
			//output the molecule file name and all (phi,psi) pairs for the residues with flexible backbones
			logPS.print(fileName+" ");
			for (int i=0; i<numInAS; i++){
				logPS.print("( "+bb.getFiPsi(m, sysStrNum, residueMap[i], 0)+" , "+bb.getFiPsi(m, sysStrNum, residueMap[i], 1)+" ) ");
			}
			logPS.println();
			logPS.flush();
		}
	}
	
	// This function checks the sterics of the given conformation;
	//  it checks backbone atoms of the given residue resNum against all
	//  backbone atoms that are in residues 0..resNum-1 for the given strand only
	//  If any two atoms overlap by more than overlapThresh then
	//  false is returned
	// strandNum is the number of the strand containing resNum
	// Hydrogens are NOT used in checking sterics in this function
	private boolean checkStericsBBonly(int strandNum, int resNum) {
	
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			
			if (isBBatom(res.atom[i])){ //backbone atom
				
				resToCheck = resNum;
				for(int w=0;w<resToCheck;w++) {
					for(int t=0;t<m.strand[strandNum].residue[w].numberOfAtoms;t++) {
						tmpAtm = m.strand[strandNum].residue[w].atom[t];
						
						if (isBBatom(tmpAtm)){
							if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
								if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
									if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
										return false;
									}
								}
							}
						}
					}
				}	
			}
		}
	
		return true; //if we are here, then everything passed the steric test
	}
	
	//This function is similar to checkStericsBBonly(), but it checks all residues in strandNum against resNum, 
	//		instead of just checking the residues up to resNum
	private boolean checkAllStericsBBonly(int strandNum, int resNum) {
		
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;
		
		for(int i=0;i<res.numberOfAtoms;i++) {
			
			if (isBBatom(res.atom[i])){ //backbone atom
				
				resToCheck = m.strand[strandNum].numberOfResidues;
				for(int w=0;w<resToCheck;w++) {
					if (w!=resNum){
						for(int t=0;t<m.strand[strandNum].residue[w].numberOfAtoms;t++) {
							tmpAtm = m.strand[strandNum].residue[w].atom[t];
							
							if (isBBatom(tmpAtm)){
								if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
									if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
										if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
											return false;
										}
									}
								}
							}
						}
					}
				}	
			}
		}
	
		return true; //if we are here, then everything passed the steric test
	}
	
	//Determines if the given atom is a backbone atom
	private boolean isBBatom(Atom at){
		
		return ((at.name.equalsIgnoreCase("N"))||(at.name.equalsIgnoreCase("CA"))
				||(at.name.equalsIgnoreCase("C"))||(at.name.equalsIgnoreCase("O")));
	}
//////////////////////////////////////////////////////////////////////////
//	End of Generate Backbones section
//////////////////////////////////////////////////////////////////////////
	
	
//////////////////////////////////////////////////////////////////////////
	
	/*//Determines how many residue positions in the system strand (pos is strand-relative numbering)
	//		are within dist from residue position pos
	public boolean [] getProxAS(int pos, float dist, boolean as[]){
		
		Residue r1 = m.strand[sysStrNum].residue[pos];
		for (int i=0; i<m.strand[sysStrNum].numberOfResidues; i++){
			if (i!=pos){
				Residue r2 = m.strand[sysStrNum].residue[i];
				if (r1.getDistSC(r2)<=dist)
					as[i] = true;
				else
					as[i] = false;
			}
		}
		
		return as;
	}*/
	
	
//////////////////////////////////////////////////////////////////////////////////
public void saveMolecule(Molecule m, String fname, float energy){
	boolean printSegID = false;
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
	
}//end of RotamerSearch class

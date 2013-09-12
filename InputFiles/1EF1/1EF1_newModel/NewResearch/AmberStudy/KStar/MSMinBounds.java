///////////////////////////////////////////////////////////////////////////////////////////////
//	MSMinBounds.java
//
//	Version:           0.3
//
//
//	  authors:
// 	  initials    name                 organization                email
//	---------   -----------------    ------------------------    ----------------------------
//	  ISG		 Ivelin Georgiev	  Duke University			  ivelin.georgiev@duke.edu
//
///////////////////////////////////////////////////////////////////////////////////////////////

/**
* Written by Ivelin Georgiev (2004-2006)
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

/* Performs two separate operations:
 * 
 * 1) Applies the Bounds/MinBounds pruning criteria: computes a lower bound on the energy of all
 * 		conformations that contain a given rotamer i_r, for each rotamer (boundKStar is false);
 * 
 * 2) Compute Ec, a lower bound on the energy of all conformations that contain a pruned rotamer, and
 * 		prunedIsSteric[], all conformations that are pruned due to unallowed sterics (boundKStar is true)
 * 
 */
public class MSMinBounds {
	
	private class MSRotBounds implements RyanComparable{
		int index;
		double Ec; //min bound for rot index
		
		public int compareTo(Object otherObject) {
			MSRotBounds otherBound = (MSRotBounds) otherObject;
			if (Ec > otherBound.Ec) return -1;
			if (Ec < otherBound.Ec) return 1;
			return 0;
		}
	}
	//two pairwise energy matrices: one for the min energies,
	//and one for the max; the definitions of the elements
	//[0][0], [i][0], [0][i], and [i][j] are the same as with K*
	private float pairwiseMinEnergyMatrix [][] = null;
	
	//eliminated rotamers at position i, for all positions
	private boolean eliminatedRotAtPos [] = null;
	
	//number of residues under consideration
	private int numSiteResidues;
	
	//for each residue, number of possible amino acids
	private int numTotalRot;
	
	//number of possible rotamers for the ligand
	int numLigRot;
	
	//offset of the given rotamer in the total rotamer set (?152?)
	int rotIndOffset[];
	
	//the number of AA tyes allowed for each AS residue
	int numAAtypes[] = null;
	
	//number of rotamers for the current AA type at the given residue
	//int numRotForAAtypeAtRes[];
	
	//this value depends on the particular value specified in the pairwise energy matrices;
	//		in KSParser, this value is 10^38;
	//entries with this particular value will not be examined, as they are not allowed;
	//note that when computing E intervals, if a steric is not allowed, (maxE-minE)=0,
	//		so no comparison with stericE is necessary there
	private float bigE = (float)Math.pow(10,38);
	
	//steric energy that determines incompatibility of a rotamer with the template
	float stericE = bigE;
	
	//size of the pairwise energy matrix
	private int PEMsize;
	
	//private double curEw = 0.0f;	//the max allowable difference from the GMEC (checkSum<=curEw should not be pruned)
	
	//PrintStream logPS = null;
	
	//stores the min bound for each rotamer
	MSRotBounds indBounds[] = null;
	
	//the minimum lower energy bound for all pruned conformations
	double Ec = bigE;
	
	//the lowest energy bound involving the current rotamer
	double curEc;
	
	double minEmin = bigE; //used in Ec computation
	
	//the number of rotamers for the current mutation sequence only
	int numRotForMut = 0;
	
	//determines if a rotamer index is a part of the current mutation sequence
	boolean rotInMutInd[];
	
	//the rotamer library
	RotamerLibrary rl = null;
	
	//The system rotamer handler
	StrandRotamers sysLR = null;
	
	//The mapping from AS position to actual residue numbers
	int residueMap[] = null;
	
	//determines if Ec and prunedIsSteric are computed for the ensemble-based bound to the total
	//		contribution of all pruned conformations
	boolean boundKStar = false;
	
	//flag that a rotamer is pruned because of a steric clash, and not because of energy difference
	boolean prunedIsSteric[] = null;
	
	double pruningE = bigE; //the lower-bound energy cutoff for pruning
	
	double Ew = 0.0; //the E window allowed from the best energy
	
	//the precomputed single and pair interval terms
	double indInt[] = null;
	double pairInt[] = null;
	
	//split flags for all rotamer pairs
	boolean splitFlags[][] = null;
	
	//determines if split flags are used
	boolean useFlags = false;
	
	//the percent of pruned non-steric rotamers
	//final double lambda = 0.3;//0.08;

	//constructor
	MSMinBounds(float arpMatrix[][], int numResInActiveSite, int numTotalRotamers,int numLigRotamers,	
			int rotamerIndexOffset[], int resMap[], StrandRotamers systemLRot, double pruneE, 
			boolean prunedRotAtRes[], boolean spFlags[][], boolean useSF, float initEw, boolean boundKS, RotamerLibrary rlP) {
		

		/*///////////////////////////////////////////////////////////////////////////////////
		try {			
			FileOutputStream fileOutputStream = new FileOutputStream("/net/grad/shaqbuzz/j"+numRuns+".txt");
			BufferedOutputStream bufferedOutputStream = new BufferedOutputStream( fileOutputStream );
			logPS = new PrintStream( bufferedOutputStream );
		}
		catch (Exception ex) {
			System.out.println("ERROR: An exception occured while opening log file");
		}numRuns++;
		///////////////////////////////////////////////////////////////////////////////////*/
		
		//size of the pairwise energy matrix: (km+1)x(km+1)
		PEMsize = arpMatrix.length;
		
		splitFlags = spFlags;
		pairwiseMinEnergyMatrix = arpMatrix;
		rotIndOffset = rotamerIndexOffset;
		eliminatedRotAtPos = prunedRotAtRes;
		residueMap = resMap;
		sysLR = systemLRot;
		rl = rlP;
		useFlags = useSF;
		boundKStar = boundKS;
		
		numSiteResidues = numResInActiveSite;		// tested with 9 AS
		numTotalRot = numTotalRotamers;				// ?152?
		numLigRot = numLigRotamers;					// 0 if no ligand
		pruningE = pruneE;
		Ew = initEw;
		
		//subtract 1 from PEMsize, as the backbone row/column in the
		//	pairwise energy matrix is not a part of the eliminated
		//	rotamer matrix
		rotInMutInd = new boolean[PEMsize-1];
		indBounds = new MSRotBounds[PEMsize-1];
		
		for (int i=0; i<(PEMsize-1); i++){
			rotInMutInd[i] = false; //true only if the given rotamer index is used in the current mutation sequence
			
			indBounds[i] = new MSRotBounds();
			indBounds[i].index = i;
			indBounds[i].Ec = bigE;
		}
		
		numAAtypes = new int[numSiteResidues];
		for (int i=0; i<numAAtypes.length; i++) //the number of AAs allowed for each AS residue
			numAAtypes[i] = sysLR.getNumAllowable(residueMap[i]);
		
		prunedIsSteric = new boolean[eliminatedRotAtPos.length];
		for (int i=0; i<prunedIsSteric.length; i++)
			prunedIsSteric[i] = false;
		
		Ec = bigE;
		curEc = 0.0;
	}
	
	public double getEc(){
		return Ec;
	}
	
	public boolean [] getPrunedSteric(){
		return prunedIsSteric;
	}
	
	//Precompute the terms for indInt[] and pairInt[]
	private void precomputeInt() {
		
		int numRes = numSiteResidues;		
		if (numLigRot!=0) //ligand is present
			numRes++;
		
		indInt = new double[numRes];
		pairInt = new double[numRes];
		
		for (int curPos=0; curPos<numRes; curPos++){
		
			curEc = 0.0;
			SumMaxIndInt(curPos);
			indInt[curPos] = curEc;			//formula term 3
			
			curEc = 0.0;
			SumSumMaxPairInt(curPos);
			pairInt[curPos] = curEc;			//formula term 4
			
			curEc = 0.0;
		}
	}

	//Compute the conformations that can be eliminated
	//Return a boolean matrix in which an element is true if
	//the corresponding r at i can be eliminated, and false otherwise
	public boolean[] ComputeEliminatedRotConf (){
		
		//precompute the terms for indInt[] and pairInt[]
		precomputeInt();
		
		boolean done = false;
		int prunedCurRun = 0;
		int numRuns = 1;
		
		while (!done) {
			
			prunedCurRun = 0;
			
			System.out.println("Current run: "+numRuns);
			
			int numRotForCurAAatPos;
			
			//Compute for the AS residues first
			for (int curPos=0; curPos<numSiteResidues; curPos++){
				
				System.out.print("Starting AS residue "+curPos);
								
				for (int AA=0; AA<numAAtypes[curPos]; AA++){
					
					System.out.print(".");
					
					int curAA = sysLR.getIndexOfNthAllowable(residueMap[curPos],AA);
				
					//find how many rotamers are allowed for the current AA type at the given residue;
					//note that ala and gly have 0 possible rotamers
					numRotForCurAAatPos = rl.getNumRotForAAtype(curAA);
					if (numRotForCurAAatPos==0)	//ala or gly
						numRotForCurAAatPos = 1;
					
					for(int curRot=0; curRot<numRotForCurAAatPos; curRot++){
						
						numRotForMut++;
						
						int index = curPos*numTotalRot + rotIndOffset[curAA] + curRot;
						rotInMutInd[index] = true; //rot index is in cur mut sequence
						
						if ((boundKStar)||(!eliminatedRotAtPos[index])){ //boundKStar or not already pruned
						
							//logPS.println((curPos*numTotalRot + rotIndOffset[curAA] + curRot));logPS.flush();
							curEc = pairwiseMinEnergyMatrix[0][0]; //initialize to Et' for each rotamer
							
							CanEliminate(curPos, curAA, curRot, numRotForCurAAatPos);
							indBounds[curPos*numTotalRot + rotIndOffset[curAA] + curRot].Ec = Math.min(indBounds[curPos*numTotalRot + rotIndOffset[curAA] + curRot].Ec, curEc); //update the lowest energy bound if necessary
						}
					}
				}
				System.out.println("done");
			}
			
			//If there is a ligand, compute MinDEE for the lig rotamers as well
			if (numLigRot!=0){
				System.out.print("Starting ligand run");
				System.out.print("..");
				for (int curRot=0; curRot<numLigRot; curRot++){
						
					numRotForMut++;
					
					int index = numSiteResidues*numTotalRot+curRot;
					rotInMutInd[index] = true; //rot index is in cur mut sequence
					
					if ((boundKStar)||(!eliminatedRotAtPos[index])){ //boundKStar or not already pruned
					
						curEc = pairwiseMinEnergyMatrix[0][0]; //initialize to Et' for each rotamer
						
						CanEliminateLig(curRot);
						indBounds[numSiteResidues*numTotalRot+curRot].Ec = Math.min(indBounds[numSiteResidues*numTotalRot+curRot].Ec, curEc); //update the lowest energy bound if necessary
					}
				}
				System.out.println("done");
			}
			
			RyanQuickSort rqs = new RyanQuickSort();
			rqs.Sort(indBounds);
			
			//Determine the pruned rotamers: all rotamers whose lower energy bound is above the given
			//	cutoff value are pruned
			double minE = bigE;
			double maxE = -bigE;
			double minEc = bigE;
			int numRot = 0;
			int numStericFromBounds = 0;
			int numTotalSteric = 0;
			for (int i=0; i<indBounds.length; i++){
				if (rotInMutInd[indBounds[i].index]){ //rot index is in current mutation
					numRot++;
					
					if (!boundKStar) {//pruning is performed
						
						if ((!eliminatedRotAtPos[indBounds[i].index])){ //not already pruned
						
							if (indBounds[i].Ec>pruningE+Ew){ //higher than the cutoff energy, so prune
								eliminatedRotAtPos[indBounds[i].index] = true;
								prunedCurRun++;
							
								if (indBounds[i].Ec >= stericE){ //pruned due to unallowed steric (from MinBounds only)
									numStericFromBounds++;
								}
								
								minE = Math.min(minE,indBounds[i].Ec);
								maxE = Math.max(maxE,indBounds[i].Ec);					
							}
						}
					}
					else {//Ec and prunedIsSteric[] are computed
					
						if (eliminatedRotAtPos[indBounds[i].index]) //check the Ec among all pruned rotamers
							minEc = Math.min(minEc,indBounds[i].Ec);
						
						if (indBounds[i].Ec >= stericE){ //pruned due to unallowed steric (from all criteria used)
							prunedIsSteric[indBounds[i].index] = true;
							numTotalSteric++;
						}
					}
				}
			}
			
			if (boundKStar)
				Ec = Math.min(Ec,minEc); //the minimum pruned lower bound (among all pruned rotamers, for boundsKStar)	
			
			
			if ((boundKStar)||(prunedCurRun==0)) //boundKStar or no rotamers pruned this run, so done
				done = true;
			else 
				numRuns++;
			
			
			System.out.println("Number of rotamers for the current sequence: "+numRot);
			
			if (!boundKStar)
				System.out.println("Number of rotamers pruned this run: "+prunedCurRun);
			System.out.println();
			
			System.out.println("minE: "+minE+" maxE: "+maxE+" pruningE: "+pruningE+" Ew: "+Ew+" Ec: "+Ec);
			if (!boundKStar)
				System.out.println("Num rotamers pruned due to unallowed sterics (from Bounds): "+numStericFromBounds);
			else
				System.out.println("Num rotamers pruned due to unallowed sterics (from all criteria): "+numTotalSteric);
			System.out.println();				
			
			/*if (boundKStar){
				for (int i=0; i<indBounds.length; i++){
					if (rotInMutInd[indBounds[i].index]){
						System.out.print("rank: "+i+" rotIndex: "+indBounds[i].index+" minBound: "+indBounds[i].Ec);
						System.out.println(" pruned: "+eliminatedRotAtPos[indBounds[i].index]);
					}
				}
			}*/
		}
		
		return eliminatedRotAtPos;
	}
	
	//Called only by ComputeEliminatedRotConf(.)
	private void CanEliminate (int posNum, int AANumAtPos, int rotNumAtPos, int numRotForCurAAatPos){
		
		double minIndVoxelE;
		double minShellResE;
		
		int index_r; //already checked if previously pruned
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r = 1 + posNum*numTotalRot + rotIndOffset[AANumAtPos] + rotNumAtPos;
		minIndVoxelE = pairwiseMinEnergyMatrix[index_r][0];					//formula term 1
		minShellResE = pairwiseMinEnergyMatrix[0][index_r];
		
		curEc += minIndVoxelE + minShellResE;//System.out.println(++count+" "+curEc);
		
		if (curEc>=stericE) //rotamer incompatible with template, so prune
			return;
		
		curEc += indInt[posNum];							//formula term 3
		curEc += pairInt[posNum];						//formula term 5	
		SumMinDiffPVE(posNum, AANumAtPos, rotNumAtPos);	//formula term 4
	}
	
	////////////////////////////////////////////////////////////////////////
	
	//Called only by CanEliminate(.)
	private void SumMinDiffPVE (int atPos, int withAA, int withRot1){		
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){
			
			if (curPos != atPos)
				
				IndMinDiffPVE(atPos, withAA, withRot1, curPos);
		}
		
		if (numLigRot!=0){ //there is a ligand
			//add the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			LigandIndMinMinPVE(atPos, withAA, withRot1);
		}
	}
	
	//Called by SumMinMinPVE(.)
	private void IndMinDiffPVE (int firstPos, int firstAA, int firstRot1, int secondPos){

		double curEmin;
		
		int index1, index2;
		int numRotForAAatPos;
		
		//r at i
		index1 = 1 + firstPos*numTotalRot + rotIndOffset[firstAA] + firstRot1;
		
		if ((boundKStar)||(!eliminatedRotAtPos[index1-1])){ //boundsKStar or not pruned
		
			//find the minimum E among all the rotamers (all the rotamers for the given AA assignment)
			//for the given residue
			for (int AA=0; AA<numAAtypes[secondPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[secondPos],AA);;
				
				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){			
					
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = 1 + secondPos*numTotalRot + rotIndOffset[curAA] + curRot;	
					
					if ((boundKStar)||(!eliminatedRotAtPos[index2-1])){ //boundsKStar or not pruned
						
						if ((boundKStar)||((!useFlags)||(!splitFlags[index1-1][index2-1]))){ //boundKStar or (not using split flags or not flagged)
						
							curEmin = pairwiseMinEnergyMatrix[index1][index2];
							minEmin = Math.min(minEmin,curEmin);
						}
					}
				}
			}		
			curEc += minEmin;//System.out.println(++count+" "+curEc+" "+minEmin);
		}
		minEmin = bigE; //re-initialize
	}
	
	//Called by SumMinMinPVE(.)
	private void LigandIndMinMinPVE(int firstPos, int firstAA, int firstRot1){
		
		double curEmin;
		
		int index1, index2;
		
		//r at i
		index1 = 1 + firstPos*numTotalRot + rotIndOffset[firstAA] + firstRot1;	
		
		if ((boundKStar)||(!eliminatedRotAtPos[index1-1])){ //boundsKStar or not pruned

			for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){
				
				//s at j (the ligand residue)
				index2 = 1 + numSiteResidues*numTotalRot + curLigPos;	
				
				if ((boundKStar)||(!eliminatedRotAtPos[index2-1])){ //boundsKStar or not pruned
					
					if ((boundKStar)||((!useFlags)||(!splitFlags[index1-1][index2-1]))){ //boundKStar or (not using split flags or not flagged)
					
						curEmin = pairwiseMinEnergyMatrix[index1][index2];
						minEmin = Math.min(minEmin,curEmin);
					}
				}
			}		
			curEc += minEmin;//System.out.println(++count+" "+curEc);
		}
		minEmin = bigE;
	}
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////		
	//Called only by CanEliminate(.)
	private void SumMaxIndInt (int withoutPos){
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){			
			if (curPos != withoutPos)			
				MaxIndInt(curPos);
		}
		
		if (numLigRot!=0){ //ther is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position j here for which to add
			if (withoutPos!=numSiteResidues) //if we are not currently checking ligand rotamers for pruning
				LigandMaxIndInt();
		}
	}
	
	//Called by SumMaxIndInt(.)
	private void MaxIndInt (int atPos){
		
		int numRotForAAatPos;
		
		for (int AA=0; AA<numAAtypes[atPos]; AA++){
			
			int curAA = sysLR.getIndexOfNthAllowable(residueMap[atPos],AA);
			
			numRotForAAatPos = rl.getNumRotForAAtype(curAA);
			if (numRotForAAatPos==0)	//ala or gly
				numRotForAAatPos = 1;
			
			for (int curRot=0; curRot<numRotForAAatPos; curRot++){		
				IndInt(atPos, curAA, curRot);
			}
		}
		
		curEc += minEmin;//System.out.println(++count+" "+curEc);
		minEmin = bigE;
	}

	//Called by MaxIndInt(.)
	private void IndInt (int atPos, int atAA, int atRot){
		
		//s at j
		int index1 = 1 + atPos*numTotalRot + rotIndOffset[atAA] + atRot;	
		
		if ((boundKStar)||(!eliminatedRotAtPos[index1-1])){ //boundsKStar or not pruned

			double minE = pairwiseMinEnergyMatrix [index1][0];			
			double minShell = pairwiseMinEnergyMatrix[0][index1];
			
			minEmin = Math.min(minEmin,minE+minShell);
		}
	}
	
	//Called by SumMaxIndInt(.)
	private void LigandMaxIndInt(){

		for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){			
			LigandIndInt(curLigPos);
		}
		
		curEc += minEmin;//System.out.println(++count+" "+curEc);
		minEmin = bigE;
	}

	//Called by LigandMaxIndInt(.)
	private void LigandIndInt (int ligRot){
		
		//s at j (the ligand residue)
		int index1 = 1 + numSiteResidues*numTotalRot + ligRot;
		
		if ((boundKStar)||(!eliminatedRotAtPos[index1-1])){ //boundsKStar or not pruned
	
			double minE = pairwiseMinEnergyMatrix [index1][0];;
			double minShell = pairwiseMinEnergyMatrix[0][index1];
			
			minEmin = Math.min(minEmin,minE+minShell);
		}
	}
	//////////////////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////////////////
	//Called only by CanEliminate(.)
	private void SumSumMaxPairInt(int withoutPos){
		
		//get the contribution from the active site residue rotamers
		for (int curPos1=0; curPos1<numSiteResidues; curPos1++){
			if (curPos1 != withoutPos){
				for (int curPos2=0; curPos2<curPos1; curPos2++){
					if (curPos2 != withoutPos){
						MaxPairInt(curPos1,curPos2);
					}
				}
			}
		}
		
		if (numLigRot!=0){ //there is a ligand
			//get the contribution from the ligand rotamers: there is only one ligand residue,
			//so there is only one position k here for which to add;
			//the range of j is the number of active site residues
			if (withoutPos!=numSiteResidues){ //if we are not currently checking ligand rotamers for pruning
				for (int curPos=0; curPos<numSiteResidues; curPos++){
					if (curPos != withoutPos){						
						LigandMaxPairInt(curPos);
					}
				}
			}
		}
	}
	
	//Called by SumSumMaxPairInt(.)
	private void MaxPairInt (int atPos1, int atPos2){
		
		int numRotForAAatPos1;
		
		for (int AA1=0; AA1<numAAtypes[atPos1]; AA1++){
			
			int curAA1 = sysLR.getIndexOfNthAllowable(residueMap[atPos1],AA1);
			
			numRotForAAatPos1 = rl.getNumRotForAAtype(curAA1);
			if (numRotForAAatPos1==0)	//ala or gly
				numRotForAAatPos1 = 1;
		
			for (int curRot1=0; curRot1<numRotForAAatPos1; curRot1++){
				
				int numRotForAAatPos2;
				
				for (int AA2=0; AA2<numAAtypes[atPos2]; AA2++){
					
					int curAA2 = sysLR.getIndexOfNthAllowable(residueMap[atPos2],AA2);;
					
					numRotForAAatPos2 = rl.getNumRotForAAtype(curAA2);
					if (numRotForAAatPos2==0)	//ala or gly
						numRotForAAatPos2 = 1;
					
					for (int curRot2=0; curRot2<numRotForAAatPos2; curRot2++){			
						PairInt(atPos1, curAA1, curRot1, atPos2, curAA2, curRot2);
					}
				}
			}
		}
		
		curEc += minEmin;//System.out.println(++count+" "+curEc);
		minEmin = bigE;
	}
	
	//Called by MaxPairInt(.)
	private void PairInt (int atPos1, int atAA1, int atRot1, int atPos2, int atAA2, int atRot2){
		
		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		int index1 = 1 + atPos1*numTotalRot + rotIndOffset[atAA1] + atRot1;//u at k
		int index2 = 1 + atPos2*numTotalRot + rotIndOffset[atAA2] + atRot2;//s at j
		
		if ((boundKStar)||((!eliminatedRotAtPos[index1-1])&&(!eliminatedRotAtPos[index2-1]))){ //boundsKStar or not pruned
			
			if ((boundKStar)||((!useFlags)||(!splitFlags[index1-1][index2-1]))){ //boundKStar or (not using split flags or not flagged)

				double minE = pairwiseMinEnergyMatrix[index1][index2];
				
				minEmin = Math.min(minEmin,minE);
			}
		}
	}
	
	//Called by SumSumMaxPairInt(.)
	private void LigandMaxPairInt (int atPos){

		int numRotForAAatPos;
		
		for (int AA=0; AA<numAAtypes[atPos]; AA++){
			
			int curAA = sysLR.getIndexOfNthAllowable(residueMap[atPos],AA);
			
			numRotForAAatPos = rl.getNumRotForAAtype(curAA);
			if (numRotForAAatPos==0)	//ala or gly
				numRotForAAatPos = 1;
		
			for (int curRot=0; curRot<numRotForAAatPos; curRot++){
				
				for (int curLigPos=0; curLigPos<numLigRot; curLigPos++){				
					LigandPairInt(atPos, curAA, curRot, curLigPos);
				}
			}
		}
		
		curEc += minEmin;//System.out.println(++count+" "+curEc);
		minEmin = bigE;
	}
	
	//Called by LigandMaxPairInt(.)
	private void LigandPairInt (int atPos, int atAA, int atRot, int ligRot){
		
		//There is a displacement: colum 0 and row 0 have special entries, 
		//so pairwise energies start from row 1, column 1
		int index1 = 1 + numSiteResidues*numTotalRot + ligRot;//u at k (the ligand residue)
		int index2 = 1 + atPos*numTotalRot + rotIndOffset[atAA] + atRot;//s at j
		
		if ((boundKStar)||((!eliminatedRotAtPos[index1-1])&&(!eliminatedRotAtPos[index2-1]))){ //boundsKStar or not pruned
			
			if ((boundKStar)||((!useFlags)||(!splitFlags[index1-1][index2-1]))){ //boundKStar or (not using split flags or not flagged)
		
				double minE = pairwiseMinEnergyMatrix[index1][index2];
				
				minEmin = Math.min(minEmin,minE);
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////
	
	///////////////////////////////////////////////////////////////////////////////////////
	//Same as CanEliminate(), just checks the ligand rotamers for pruning
	//Called by ComputeEliminatedRotConf()
	private void CanEliminateLig (int curLigRot){
		
		double minIndVoxelE;
		double minShellResE;
		
		int index_r; //already checked if previously pruned
		
		//In the energy matrix, column 0 gives the individual energies for each r at i;
		//skip row 0, as the individual energies start from row 1 (and are all in column 0)
		index_r = 1 + numSiteResidues*numTotalRot + curLigRot;
		minIndVoxelE = pairwiseMinEnergyMatrix[index_r][0]; 					//formula term 1
		minShellResE = pairwiseMinEnergyMatrix[0][index_r];
		
		curEc += minIndVoxelE + minShellResE;
		
		if (curEc>=stericE) //rotamer incompatible with template, so prune
			return;
		
		curEc += indInt[numSiteResidues];							//formula term 3
		curEc += pairInt[numSiteResidues];						//formula term 5
		SumMinDiffPVELig(curLigRot);							//formula term 4
	}
	
	//Same as SumMinDiffPVE(), just checks the ligand rotamers for pruning;
	//Called by CanEliminateLig()
	private void SumMinDiffPVELig (int withRot1){
		
		//get the contribution from the active site residue rotamers
		for (int curPos=0; curPos<numSiteResidues; curPos++){			
			IndMinDiffPVELig(withRot1, curPos);
		}
	}
	
	//Same as IndMinDiffPVE(), just checks the ligand rotamers for pruning
	//Called by SumMinDiffPVELig()
	private void IndMinDiffPVELig (int firstRot1, int secondPos){
		
		double curEmin;
		
		int index1, index2;
		int numRotForAAatPos;
		
		//r at i
		index1 = 1 + numSiteResidues*numTotalRot + firstRot1;
		
		if ((boundKStar)||(!eliminatedRotAtPos[index1-1])){ //boundsKStar or not pruned
		
			//find the minimum E among all the rotamers (all the rotamers for the given AA assignment)
			//for the given residue
			for (int AA=0; AA<numAAtypes[secondPos]; AA++){
				
				int curAA = sysLR.getIndexOfNthAllowable(residueMap[secondPos],AA);
				
				numRotForAAatPos = rl.getNumRotForAAtype(curAA);
				if (numRotForAAatPos==0)	//ala or gly
					numRotForAAatPos = 1;
				
				for (int curRot=0; curRot<numRotForAAatPos; curRot++){			
					
					//There is a displacement: column 0 and row 0 have special entries, 
					//so pairwise energies start from row 1, column 1
					
					//s at j
					index2 = 1 + secondPos*numTotalRot + rotIndOffset[curAA] + curRot;
					
					if ((boundKStar)||(!eliminatedRotAtPos[index2-1])){ //boundsKStar or not pruned
						
						if ((boundKStar)||((!useFlags)||(!splitFlags[index1-1][index2-1]))){ //boundKStar or (not using split flags or not flagged)
							curEmin = pairwiseMinEnergyMatrix[index1][index2];				
							minEmin = Math.min(minEmin,curEmin);
						}
					}
				}
			}			
			curEc += minEmin;
		}
		minEmin = bigE;
	}
}

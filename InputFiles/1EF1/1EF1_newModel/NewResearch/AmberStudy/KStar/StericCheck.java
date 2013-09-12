///////////////////////////////////////////////////////////////////////////////////////////////
//	StericCheck.java
//
//	Version:           0.3
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
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

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Hashtable;
import java.math.*;

public class StericCheck {	
	
	//numAS residues
	int numInAS = -1;
	
	//the cur AA at each residue
	int curAANum[] = null;
	
	//the cur AA for the ligand
	int ligNum = -1;
	
	//the mapping from residue number to AS num and from AS num to res num
	int curResToASMap[] = null;
	int residueMap[] = null;
	
	//the current molecule
	Molecule m = null;
	
	//the system and ligand strand numbers
	int sysStrNum = -1;
	int ligStrNum = -1;
	
	//the rotamer library
	RotamerLibrary rl = null;
	
	//the sys and lig rotamer handlers
	StrandRotamers sysLR = null;
	StrandRotamers ligROT = null;
	
	//the overlap threshold for the steric check
	double overlapThresh = 1.5;
	
	//determines whether hydrogens are used in the steric check
	boolean hSteric = false;
	
	//the MinDEE pruning matrix
	boolean eliminatedRot[] = null;
	
	//num rotamers for each residue
	int numRotForRes[] = null;

	//the number of conformations not examined yet
	BigInteger numConfsLeft = new BigInteger("0");
	
	//the number of conformations above each level (from level i+1 to the top level,
	//	which is the ligand)
	BigInteger numConfsAboveLevel[] = null;
	
	//the number of conformations pruned by the steric filter
	BigInteger numConfsPrunedByS = new BigInteger("0");
	
	//determines whther a ligand is present
	boolean ligPresent = false;
	
	//PrintStream logPS = null;
	
	StericCheck (int curAANumP[],int curResToASMapP[],int residueMapP[],boolean eliminatedRotAtPosRedP[],
			int numRotForResP[],Molecule mP, double overlapThreshP,	boolean hS, BigInteger numConfsLeftP,
			BigInteger numConfsAboveLevelP[], int sysStrNumP, StrandRotamers sysLRP, int ligStrNumP, 
			StrandRotamers ligROTP, int curLigNum, RotamerLibrary rlP){
		
		curAANum = curAANumP;
		curResToASMap = curResToASMapP;
		residueMap = residueMapP;
		numInAS = residueMap.length;
		m = mP;
		overlapThresh = overlapThreshP;
		hSteric = hS;
		eliminatedRot = eliminatedRotAtPosRedP;
		numRotForRes = numRotForResP;
		numConfsLeft = numConfsLeftP;
		numConfsAboveLevel = numConfsAboveLevelP;
		numConfsPrunedByS = new BigInteger("0");
		sysStrNum = sysStrNumP;
		sysLR = sysLRP;
		ligStrNum = ligStrNumP;
		ligROT = ligROTP;
		ligNum = curLigNum;
		ligPresent = true;
		rl = rlP;
		
		//logPS = logPSP;
	}
	
	StericCheck (int curAANumP[],int curResToASMapP[],int residueMapP[],boolean eliminatedRotAtPosRedP[],
			int numRotForResP[],Molecule mP, double overlapThreshP, boolean hS, BigInteger numConfsLeftP, 
			BigInteger numConfsAboveLevelP[], int sysStrNumP, StrandRotamers sysLRP, RotamerLibrary rlP){
		
		curAANum = curAANumP;
		curResToASMap = curResToASMapP;
		residueMap = residueMapP;
		numInAS = residueMap.length;
		m = mP;
		overlapThresh = overlapThreshP;
		hSteric = hS;
		eliminatedRot = eliminatedRotAtPosRedP;
		numRotForRes = numRotForResP;
		numConfsLeft = numConfsLeftP;
		numConfsAboveLevel = numConfsAboveLevelP;
		numConfsPrunedByS = new BigInteger("0");
		sysStrNum = sysStrNumP;
		sysLR = sysLRP;
		ligPresent = false;
		rl = rlP;
		
		//logPS = logPSP;
	}
	
	//Return the number of conformations not examined yet
	public BigInteger getNumConfsLeft() {
		return numConfsLeft;
	}
	
	//Return the number of conformations pruned by the steric filter
	public BigInteger getNumConfsPrunedByS(){
		return numConfsPrunedByS;
	}
	
	//Sets the number of conformations remaining to be examined
	public void setNumConfsLeft(BigInteger numCLeft){
		numConfsLeft = numCLeft;
	}
	
	//Sets up the steric check for the given partial (full) conformation: curConf[] has been
	//	assigned for levels 0...(curTopLevel-1) and curNode should be applied at curTopLevel;
	//Returns true if the partial conformation is sterically allowed
	//Called by A*
	public boolean checkAllowedSteric (int curTopLevel, int curConf[], int curNode){		
		
		//As the rotamers given to A* are only the non-pruned ones, there is a difference between the
		//	rotamer numbers returned by A* and the actual rotamer numbers for each residue (that is,
		//	A* may return rot 4 for res 3, but rot 3 for res 3 may be pruned, and so the actual number
		//	of the rot to be applied for res 3 is 5)
		int curPruningInd = 0;
		int curRotInd, compInd;
		int conf[] = new int[curTopLevel+1];
		for (int curRes=0; curRes<=curTopLevel; curRes++){//compute the actual rot numbers for levels 0...curTopLevel
			curRotInd = 0;
			for (int curRot=0; curRot<numRotForRes[curRes]; curRot++){
				if (!eliminatedRot[curPruningInd]){
					if (curRes==curTopLevel) //the rotamer (curNode) for curTopLevel is not a part of curConf[]
						compInd = curNode;
					else
						compInd = curConf[curRes];
					
					if (curRotInd==compInd)
						conf[curRes] = curRot;
					curRotInd++;
				}
				curPruningInd++;
			}
		}

		//Backup the atom coordinates, so that they can be restored after the steric check, as 
		//		applyRotamer() changes both the actualCoordinates[] and the atom coordinates,
		//		so we cannot restore the original position just using m.updateCoordinates()
		m.backupAtomCoord();
		
		//Apply the rotamers of the current partial conformation (up to level (curTopLevel-1))
		int curAS = 0;
		boolean applyLig = false;
		
		//If there is a ligand and curTopLevel is the ligand level (full conformation), apply the lig rotamer;
		//	otherwise, curTopLevel is an AS residue
		if ((ligPresent)&&(curTopLevel==numInAS)){ //apply the ligand rotamer
			if (rl.getNumRotForAAtype(ligNum)!=0){//not GLY or ALA
				ligROT.applyRotamer(m, 0, conf[curTopLevel]);//the ligand level
			}
			applyLig = true;
		}
		for (int curRes=0; curRes<m.strand[sysStrNum].numberOfResidues; curRes++){
			if (curAS<curTopLevel){ //apply for AS res 0...(curTopLevel-1)
				if (curResToASMap[curRes]!=-1){//make a change only to the AS residues: use the native type for the other residues
										
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0){//not GLY or ALA
						sysLR.applyRotamer(m, curRes, conf[curAS]);
					}
					curAS++; //prepare the next AS residue
				}
			}
			else if (!applyLig) { //we need to apply an AS rot at curTopLevel, as there is no ligand
				if (curResToASMap[curRes]!=-1){
					if (rl.getNumRotForAAtype(curAANum[residueMap[curAS]])!=0)//not GLY or ALA
						sysLR.applyRotamer(m, curRes, conf[curTopLevel]);
					break;
				}
			}
			else //we have already applied all of the rotamers for the given partial conformation
				break;
		}
		
		/*logPS.println("curTopLevel "+curTopLevel+" curNode "+curNode+" curConf ");
		for (int i=0;i<=curTopLevel;i++)logPS.print(conf[i]+" ");logPS.println();
		if (!ligPresent){
			for (int i=0;i<=curTopLevel;i++)logPS.print(sysLR.getCurRotNum(residueMap[i])+" ");logPS.println();logPS.flush();
		}
		else{
			for (int i=0;i<curTopLevel;i++)logPS.print(sysLR.getCurRotNum(residueMap[i])+" ");
			logPS.print(ligROT.getCurRotNum(0));logPS.println();logPS.flush();
		}*/
		
		boolean allowedSteric = true;
		//Do the steric checks
		if ((ligPresent)&&(curTopLevel==numInAS)){ //check the ligand (which is at the top level) against all other residues
			if (hSteric)
				allowedSteric = RS_CheckAllStericsWithH(ligStrNum,0);
			else
				allowedSteric = RS_CheckAllSterics(ligStrNum,0);
		}
		else {
			if (hSteric)
				allowedSteric = RS_CheckAllStericsWithH(sysStrNum,residueMap[curTopLevel]);
			else
				allowedSteric = RS_CheckAllSterics(sysStrNum,residueMap[curTopLevel]);
		}
		
		m.restoreAtomCoord(); //restore the atom coordinates
		m.updateCoordinates(); //restore the actualCoordinates
		
		if (!allowedSteric){ //decrease the number of remaining conformations
			numConfsLeft = numConfsLeft.subtract(numConfsAboveLevel[curTopLevel]);
			numConfsPrunedByS = numConfsPrunedByS.add(numConfsAboveLevel[curTopLevel]);
		}
		//logPS.println("allowedSteric "+allowedSteric+" confsLeft "+numConfsLeft+" confsPrunedByS "+numConfsPrunedByS);
		//logPS.println();logPS.flush();
		
		return allowedSteric;
	}

	//The two functinos below differ from the corresponding functions in RotamerSearch
	//		in the order the system and ligand residues are handled: for these functions,
	//		it is assumed that the system residues are assigned first and then the
	//		ligand residue (if present); in RotamerSearch the CheckSterics assumes that
	//		the assignement is done for the ligand first, and then for the AS residues
	
	// This function is similar to RS_CheckStericsWithH
	// This version checks all residues against the target residue rather
	//  than just checking residues up to the specified one in the strand;
	// The AS residues with res numbers greater than the resNum for the
	//	target residue are not checked, as they have not been assigned yet
	private boolean RS_CheckAllStericsWithH(int strandNum, int resNum) {
		
		//The mapping from AS res to the system strand residue numbering
		boolean curASToResMap[] = new boolean[m.strand[sysStrNum].numberOfResidues];
		for (int i=0; i<curASToResMap.length; i++)
			curASToResMap[i] = false;
		
		for (int i=0; i<residueMap.length; i++)
			curASToResMap[residueMap[i]] = true;
	
		Residue res = m.strand[strandNum].residue[resNum];
	
		Atom tmpAtm = null;
		int resToCheck = 0;

		for(int i=0;i<res.numberOfAtoms;i++) {
			for(int q=0;q<m.numberOfStrands;q++) {
				resToCheck = m.strand[q].numberOfResidues;
				for(int w=0;w<resToCheck;w++) {
					if(!((q==strandNum) && (w==resNum))) {//not the same residue
						if ((!((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)&&(curASToResMap[w])))&&
								(!((q==ligStrNum)&&(strandNum==sysStrNum)))){//not an AS residue with a bigger res number AND not the ligand
							for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
								tmpAtm = m.strand[q].residue[w].atom[t];
								if ((res.atom[i].distance(tmpAtm)) < (((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
									if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
										return false;
									}
								}
							}
						}
						else if ((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)&&(curASToResMap[w])){ //an unassigned AS residue, so check only the backbone atoms of that residue
							for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
								tmpAtm = m.strand[q].residue[w].atom[t];
								if ((tmpAtm.name.equalsIgnoreCase("C"))||(tmpAtm.name.equalsIgnoreCase("O"))||
										(tmpAtm.name.equalsIgnoreCase("N"))||(tmpAtm.name.equalsIgnoreCase("H"))||
										(tmpAtm.name.equalsIgnoreCase("CA"))||(tmpAtm.name.equalsIgnoreCase("HA"))){
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
			}
		}
	
		// If you got here then everything passed
		return true;
	}

	
	// This function is similar to RS_CheckSterics
	// This version checks all residues against the target residue rather
	//  than just checking residues up to the specified one in the strand
	// The AS residues with res numbers greater than the resNum for the
	//	target residue are not checked, as they have not been assigned yet;
	// Also, when resNum is a residue in the AS, the ligand is not checked,
	//	as it has not been assigned yet
	private boolean RS_CheckAllSterics(int strandNum, int resNum) {
		
		//The mapping from AS res to the system strand residue numbering
		boolean curASToResMap[] = new boolean[m.strand[sysStrNum].numberOfResidues];
		for (int i=0; i<curASToResMap.length; i++)
			curASToResMap[i] = false;
		
		for (int i=0; i<residueMap.length; i++)
			curASToResMap[residueMap[i]] = true;
	
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
						if(!((q==strandNum) && (w==resNum))) {//not the same residue
							if ((!((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)&&(curASToResMap[w])))&&
									(!((q==ligStrNum)&&(strandNum==sysStrNum)))){//not an AS residue with a bigger res number AND not the ligand
								for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
									tmpAtm = m.strand[q].residue[w].atom[t];
									if (!(tmpAtm.elementType.equalsIgnoreCase("H"))) {
										if ((res.atom[i].distance(tmpAtm) < ((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
											if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
												/*if((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)){logPS.println("\n Clash: "+tmpAtm.moleculeAtomNumber+" res "+w+" atom "+t+" with res "+resNum+" atom "+i);logPS.flush();
												logPS.println("type "+tmpAtm.elementType+" "+res.atom[i].elementType);logPS.flush();
												writeCurrentMolecule("what.pdb");
												System.exit(0);}*/
												return false;
											}
										}
									}
								}
							}
							else if ((strandNum==sysStrNum)&&(q==strandNum)&&(w>resNum)&&(curASToResMap[w])){ //an unassigned AS residue, so check only the backbone atoms of that residue
								for(int t=0;t<m.strand[q].residue[w].numberOfAtoms;t++) {
									tmpAtm = m.strand[q].residue[w].atom[t];
									if ((tmpAtm.name.equalsIgnoreCase("C"))||(tmpAtm.name.equalsIgnoreCase("O"))||
											(tmpAtm.name.equalsIgnoreCase("N"))||(tmpAtm.name.equalsIgnoreCase("CA"))){ //the hydrogens are not included in this check
										if ((res.atom[i].distance(tmpAtm)) < (((tmpAtm.radius + res.atom[i].radius)/100.0) - overlapThresh)){
											if (!(res.atom[i].bondedTo(tmpAtm.moleculeAtomNumber))) {
												/*logPS.println("\n Clash: "+tmpAtm.moleculeAtomNumber+" res "+w+" atom "+t+" with res "+resNum+" atom "+i);logPS.flush();
												logPS.println("type "+tmpAtm.elementType+" "+res.atom[i].elementType);logPS.flush();
												writeCurrentMolecule("what.pdb");
												System.exit(0);*/
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
		}
	
		// If you got here then everything passed
		return true;
	}
}

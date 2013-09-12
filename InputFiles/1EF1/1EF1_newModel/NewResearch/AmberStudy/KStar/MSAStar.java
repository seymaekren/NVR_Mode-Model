///////////////////////////////////////////////////////////////////////////////////////////////
//	MSAStar.java
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
* Written by Ivelin Georgiev (2004-2007)
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

/*
 * Uses A* search for single or multiple mutation sequences simultaneously to return the minimum-energy conformation; 
 * 		each consecutive run returns the next-lowest min
 * 
 * The single-mutation-sequence A* uses the steric filter object; the steric filter cannot be used with multiple-mutation-sequence A*
 * 
 * The multiple-mutation-sequence A* can allow the generation only of sequences with up to numMaxChanges mutations from the wildtype
 * 
 * The parameters from the calling program need to be modified to fit the reduced matrices here:
 * 		only the entries for the possible mutations (AA assignments) are stored in the matrices;
 * 
 * The expansion queue is returned after each call, so that the next call can start with the
 * 		saved expansion queue: this allows to return the next-lowest-energy conformation
 * 
 * 
 * If A* is not able to find a non-pruned node for a given level (all of the possible nodes for that level
 * 		have been pruned by DEE), then it returns immediately, with the node for that level set to -1
 * 
 */


public class MSAStar {

	//number of residues under consideration
	private int numTreeLevels;
	
	//number of rotamers possible for each residue (given by the assigned AA type)
	private int numNodesForLevel[] = null;
	
	//the total number of possible rotamers for the given mutation
	private int numTotalNodes;
	
	//the offset in the array index for each level
	private int nodeIndexOffset[] = null;
		
	//eliminated rotamers at residue i, for all residues
	//private boolean eliminatedNodesAtLevel [] = null;
	
	//the current sequence: the number of the corrsponding rotamer for each level assigned so far 
	private int curConf[] = null;

	//the reduced min pairwise energy matrix
	private float pairwiseMinEnergyMatrix [][] = null;
	
	//the leaf nodes visible from the expansions
	private ExpansionQueue curExpansion;
	
	int numExpanded = 0;
	
	int topL = 0;
	int numTopL = 0;
	
	//the steric check filter (*only* for single-sequence A*)
	StericCheck stericF = null;
	
	//constructor
	/*
	 * We assume that the parameters supplied (energy and DEE information) have already been modified
	 * 		to consider only the residues and rotamers for the possible mutations; i.e., the matrices are of
	 * 		reduced size
	 */
	MSAStar (int treeLevels, int numRotForRes[], float arpMatrixRed[][], StericCheck stF){
	
		numTreeLevels = treeLevels;
		
		numNodesForLevel = new int [numTreeLevels];
		nodeIndexOffset = new int [numTreeLevels];
		numTotalNodes = 0;
		
		for (int i=0; i<numTreeLevels; i++){
			nodeIndexOffset[i] = numTotalNodes;
			numNodesForLevel[i] = numRotForRes[i];
			numTotalNodes += numNodesForLevel[i];
		}
	
		//the min energy matrix: the last column contains the intra-energy for each rotamer; the last row
		//		contains the shell-residue energy for each rotamer
		pairwiseMinEnergyMatrix = new float [numTotalNodes+1][numTotalNodes+1];
		
		for (int i=0; i<numTotalNodes+1; i++){
			for (int j=0; j<numTotalNodes+1; j++){
				pairwiseMinEnergyMatrix[i][j] = arpMatrixRed[i][j];
			}			
		}
		
		//the current expansion list
		curExpansion = new ExpansionQueue();
		
		//the current conformation
		curConf = new int [numTreeLevels];
		for (int i=0; i<numTreeLevels; i++){
			curConf[i] = -1;
		}
		
		stericF = stF;
	}
	
	//Find the lowest-energy conformation and return an array with the number of the corresponding
	//		chosen rotamer for each residue;
	//		the mapping to the rotamer information is done by the calling procedure;
	//		the chosen conformation should be marked so that the next call to AStar will return the
	//			conformation with the next lowest energy value, and so on
	/*
	 * Look at the minimum value in the expansion list; determine the node level corresponding to this value;
	 * 		expand the node into the next level; update the expansion list (adding the new nodes and deleting
	 * 		the expanded node) and determine the f(n)=g(n)+h(n) scores for the new nodes
	 * 
	 * To get the next lowest conformation, the state of the expansion queue is saved after each complete run
	 * 		of A*: after the lowest-energy conformation is found, the queue is saved, and A* returns this
	 * 		conformation. To find the second conformation, A* runs on the saved queue, and this is repeated
	 * 		for all subsequent conformations
	 */
	public int[] doAStar (boolean run1, int numMaxChanges, int nodesDefault[], boolean prunedNodes[],
			StrandRotamers sysLR, String resDefault[], int numForRes[], int residueMap[], boolean singleSeq, RotamerLibrary rl){
		
		int curLevelNum = 0;
		double hScore;
		double gScore;
		double fScore;
		QueueNode expNode = null;
		QueueNode newNode = null;
		
		for (int i=0; i<numTreeLevels; i++){ //initialize for this run
			curConf[i] = -1;
		}

		if (run1) {//if this is the first run of A*, then we need to set-up the empty expansion queue
			
			//initially, we are at level zero; all nodes at that level are visible;
			//	compute their f(n) score and add to the expansion queue
			for (int curNode=0; curNode<numNodesForLevel[0]; curNode++){
				
				if ((stericF==null)||(stericF.checkAllowedSteric(0,curConf,curNode))){//do not do a steric check if backbone minimization

					curConf[0] = curNode;//this is the only node in the current conformation
					
					//compute f for the current node
					gScore = gCompute (curLevelNum, curConf);
					hScore = hCompute (curLevelNum, curConf);
					fScore = gScore + hScore;
					
					//create a queueNode with the corresponding information
					newNode = new QueueNode (curNode, 0, curConf, fScore);
					
					//insert in the expansion list
					curExpansion.insert(newNode);
				}
			}
			if (curConf[0]==-1) //no sterically allowed nodes at the first residue, so no possible conformations				
				return curConf;				
		}

		boolean done = false;		
		//start forming the conformation by expanding the lowest-valued node and updating the expansion queue
		/*
		 * While not at the last level
		 * 	For the current minimum node in the queue
		 * 		For all possible nodes at the next level
		 * 			If the steric is allowed
		 * 				Compute the f(n) scores; Add to the queue
		 * 		Delete the expanded node from the queue
		 */
		while (!done) {	
			
			if(!run1){
				for (int i=0; i<numTreeLevels; i++) //reinitialize for each consecutive node to be expanded
					curConf[i] = -1;
			}
			
			expNode = curExpansion.curFront;//get the current min node
			
			if (expNode==null){//the queue is empty
				return curConf; //so return a sequence of -1's to flag the stop of the search
			}
			
			else { //the queue is not empty
				
				printState(expNode);
				
				for (int i=0; i<=expNode.level; i++){//get the corresponding conformation leading to this node
					curConf[i] = expNode.confSoFar[i];
				}
				curLevelNum = expNode.level;//get the corresponding level
				
				//if the current min node is at the last level, we have found a full conformation
				if (curLevelNum==numTreeLevels-1){
					curExpansion.delete(expNode);//delete this node to set-up for the next min conformation (next run of A*)
					done = true;
				}
				else {//we are not at the last level, so we can expand the current node to the next level
					
					curLevelNum++;	//since the new nodes are at the next level, increment the current level number
					
					
					int numChanges = 0;
					int curPruningInd = 0;
					if (!singleSeq) { //multiple sequences
						//Check the number of differences from the default node sequence: should not exceed numMaxChanges
						for (int i=0; i<curLevelNum; i++){ //first, get the actual conf for the assigned curConf
							int curNodeInd = 0;
							int actualNode = -1;
							for (int curRot=0; curRot<numForRes[i]; curRot++){
								if (!prunedNodes[curPruningInd]){								
									if (curNodeInd==curConf[i])
										actualNode = curRot;
									curNodeInd++;
								}
								curPruningInd++;
							}
							int index = getAAIndex(actualNode,i,resDefault,sysLR,residueMap,rl);
							if (index!=nodesDefault[i])
								numChanges++;
						}
					}
					
					for (int curNode=0;curNode<numNodesForLevel[curLevelNum];curNode++){
						
						int numChanges2 = numChanges;
						if (!singleSeq) //multiple sequences, so compute the difference from WT
							numChanges2 = getNewNumChanges(curPruningInd,curLevelNum,numForRes,prunedNodes,curNode,resDefault,sysLR,residueMap,nodesDefault,numChanges,rl);
						
						
						if ((singleSeq)||(numChanges2<=numMaxChanges)){ //continue only if (singleSeq) or (num differences does not exceed the max one)
							
							if ((stericF==null)||(stericF.checkAllowedSteric(curLevelNum,curConf,curNode))){//do not do a steric check if backbone minimization
		
								curConf[curLevelNum] = curNode;//add curNode to the conformation so far
								
								//compute f for the current node
								gScore = gCompute (curLevelNum, curConf);
								hScore = hCompute (curLevelNum, curConf);
								fScore = gScore + hScore;
								
								//create a queueNode with the corresponding information
								newNode = new QueueNode (curNode, curLevelNum, curConf, fScore);
								
								//insert in the expansion list
								curExpansion.insert(newNode);
							}
						}
					}

					//delete the expanded node from the queue, since it has already contributed
					curExpansion.delete(expNode);
				}	
			}
		}
		
		return curConf;
	}
	
	private int getAAIndex(int rotIndex, int curRes, String resDefault[], StrandRotamers sysLR, int residueMap[], RotamerLibrary rl){
		
		int rotSum = 0;
		for (int i=0; i<sysLR.getNumAllowable(residueMap[curRes]); i++){
			int curRot = rl.getNumRotamers(rl.getAAName(sysLR.getIndexOfNthAllowable(residueMap[curRes],i)));
			if (curRot==0) //GLY or ALA
				curRot = 1;
			rotSum += curRot;			
			if (rotSum>rotIndex)
				return (sysLR.getIndexOfNthAllowable(residueMap[curRes],i));
		}
		return -1; //the AA in the last position
	}
	
	private int getNewNumChanges(int curPruningInd, int curLevelNum, int numForRes[], boolean prunedNodes[], int curNode, String resDefault[],
			StrandRotamers sysLR, int residueMap[], int nodesDefault[], int numCurChanges, RotamerLibrary rl){
		
		//Get the actual node number for curNode
		int actualNode = -1;
		int index = -1;
		int curPruningInd2 = curPruningInd;
		if (curLevelNum!=(numTreeLevels-1)) {//not the ligand level	
			int curNodeInd = 0;
			for (int curRot=0; curRot<numForRes[curLevelNum]; curRot++){
				if (!prunedNodes[curPruningInd2]){								
					if (curNodeInd==curNode)
						actualNode = curRot;
					curNodeInd++;
				}
				curPruningInd2++;
			}							
			
			index = getAAIndex(actualNode,curLevelNum,resDefault,sysLR,residueMap,rl);
			if (index!=nodesDefault[curLevelNum])
				numCurChanges++;
		}
		
		return numCurChanges;
	}
	
	//Updates and prints the state of the queue
	private void printState(QueueNode expNode){
		
		numExpanded++;
		if (expNode.level+1>topL){
			topL = expNode.level+1;
			numTopL = 1;
		}
		else if (expNode.level+1==topL)
			numTopL++;
		
		if((numExpanded%1000)==0){
			System.out.print(curExpansion.numNodes+" "+expNode.fScore+" "+expNode.level+" ");
			for (int i=0;i<expNode.level;i++)System.out.print(expNode.confSoFar[i]+" ");
			System.out.println(topL+" "+numTopL);
		}
	}
//////////////////////////////////////////////////////////////////////////
	
//////////////////////////////////////////////////////////////////////////
	
	//Compute the h(n) score for the new node expanded by expNode
	//		called by doAStar(.)
	private double hCompute (int dLevel, int conf[]){
		
		double hn = 0.0;

		for (int curLevel=dLevel+1;curLevel<numTreeLevels;curLevel++){
			hn += EnergyAtLevel(dLevel, curLevel, conf);
		}
			
		return hn;
	}
	
	//Called by hCompute(.)
	private double EnergyAtLevel(int topLevel, int curLevel, int conf[]){
		
		double minE = Math.pow(10,30);
		double curE;		
		int index1;
		
		double minShellResE;
		double minIndVoxE;			//formula term 1
		double sumMinPairE;			//formula term 2
		double sumMinMinPairE;		//formula term 3
		
		for (int i1=0; i1<numNodesForLevel[curLevel];i1++){		//the rotamers at j
			
			index1 = nodeIndexOffset[curLevel]+i1;	//the index of s at j
				
			minShellResE = pairwiseMinEnergyMatrix[numTotalNodes][index1];//the shell-residue E is in the last row
			minIndVoxE = pairwiseMinEnergyMatrix[index1][numTotalNodes];//the intra-energy is in the last column
			sumMinPairE = hSumMinPVE (topLevel, index1, conf);
			sumMinMinPairE = sumMinMinPVE(topLevel+1, curLevel, index1);
			
			curE = minShellResE + minIndVoxE + sumMinPairE + sumMinMinPairE;
			if (curE<minE)		//compare to the min energy found so far
				minE = curE;
		}
		
		return minE;
		
	}
	
	//Called by EnergyAtLevel(.)
	private double hSumMinPVE (int topLevel, int index1, int conf[]){
		
		double sum = 0.0;
		int index2;
		
		for (int level=0; level<=topLevel; level++){
			
			index2 = nodeIndexOffset[level] + conf[level]; //the index of r at i
			
			sum += pairwiseMinEnergyMatrix[index2][index1]; //the pairwise energy between the two nodes
		}
		
		return sum;
	}
	
	//Called by EnergyAtLevel(.)
	private double sumMinMinPVE(int startLevel, int jLevel, int firstIndex){
		
		double sum = 0.0;
		for (int level=jLevel+1; level<numTreeLevels; level++){
			sum += indMinMinPVE(level, firstIndex);
		}
		
		return sum;
	}
	
	//Called by sumMinMinPVE(.)
	private double indMinMinPVE (int kLevel, int firstIndex){
		
		double minEn = Math.pow(10,30);
		double curEn;
		int secondIndex;
		
		for (int i2=0; i2<numNodesForLevel[kLevel]; i2++){ //u at k
			
			secondIndex = nodeIndexOffset[kLevel]+i2;
				
			curEn = pairwiseMinEnergyMatrix[firstIndex][secondIndex];
			if (curEn<minEn)
				minEn = curEn;
		}
		
		return minEn;
	}
//////////////////////////////////////////////////////////////////////////
	
//////////////////////////////////////////////////////////////////////////
	
	//Compute the g(n) score for the new node expanded by expNode;
	//		called by doAStar(.)
	private double gCompute (int dLevel,int conf[]){
		
		double gn = 0.0;
		int index1;
		
		double minShellResE;
		double minIndVoxE;		//formula term 1
		double sumMinPairE;		//formula term 2
		
		for (int curLevel=0; curLevel<=dLevel; curLevel++){//copute using the formula
			
			index1 = nodeIndexOffset[curLevel] + conf[curLevel];//index of r at i
			
			minShellResE = pairwiseMinEnergyMatrix[numTotalNodes][index1];//the shell-residue E is in the last row
			minIndVoxE = pairwiseMinEnergyMatrix[index1][numTotalNodes];//the intra-energy is in the last column
			
			sumMinPairE = gSumMinPVE(dLevel, curLevel+1, index1, conf);
			
			gn += (minShellResE + minIndVoxE + sumMinPairE);
		}
		
		return gn;
	}
	
	//Called by gCompute(.)
	private double gSumMinPVE(int topLevel, int startLevel, int index1, int conf[]){
		
		int index2;
		double sum = 0.0;
		
		for (int level=startLevel; level<=topLevel; level++){
			
			index2 = nodeIndexOffset[level] + conf[level];	//s at j
			sum += pairwiseMinEnergyMatrix[index1][index2];
		}
		
		return sum;
	}
}
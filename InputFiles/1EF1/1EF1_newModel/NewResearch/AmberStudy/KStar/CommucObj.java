////////////////////////////////////////////////////////////////////////////////////////////
// CommucObj.java
//
//  Version:           0.1
//
//
// authors:
//    initials    name            organization                email
//   ---------   --------------  ------------------------    ------------------------------
//     RHL        Ryan Lilien     Dartmouth College           ryan.lilien@dartmouth.edu
//
////////////////////////////////////////////////////////////////////////////////////////////

/*
 *
 * Written by Ryan Lilien (2002-2004)
 * 
 * The CommucObj class is a data structure used in communication between the master
 *  and slave nodes.
 * It is a 'no frills' object. It's basically just a data container.
 * It allows the master to specify what type of search the slave should perform
 *  and it allows the slave to return the result of the computation to the master.
 *
 */

/*
  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

 Contact Info:
   Bruce Donald
   HB 6211
   Computer Science Department
   Dartmouth College
   Hanover, NH 03755
   brd@cs.dartmouth.edu

 If you use or publish any results derived from the use of this program please cite:
   Ryan H. Lilien, Brian W. Stevens, Amy C. Anderson,
   Bruce R. Donald. "A Novel Ensemble-Based Scoring and
   Search Algorithm for Protein Redesign, and its
   Application to Modify the Substrate Specificity of
   the Gramicidin Synthetase A Phenylalanine Adenylation
   Enzyme." Proc. The Eighth Annual International Conference
   on Research in Computational Molecular Biology (RECOMB),
   San Diego, pp 46-57 (2004). 

 Copyright (C) 2004  Ryan H. Lilien and Bruce R. Donald

 <signature of Bruce Donald>, 16 May, 2004
 Bruce Donald, Professor of Computer Science
*/

import java.io.Serializable;
import java.util.Hashtable;
import java.math.BigDecimal;
import java.math.BigInteger;

public class CommucObj implements Serializable
{
	public class ConfInfo implements Serializable {
		String AA[] = null;
		int rot[] = null;
		float minBound = 0.0f;
		float unMinE = 0.0f;
		float minE = 0.0f;
		
		ConfInfo (int treeLevels){
			rot = new int[treeLevels];
			AA = new String[treeLevels];
		}
	}
	
	ConfInfo confSeq[] = null;
	int E_searchNumConfsTotal = 0;
	int E_searchNumConfsPrunedByE = 0;
	int E_searchNumConfsPrunedByS = 0;
	int E_searchNumConfsEvaluated = 0;
	int E_searchNumConfsLeft = 0;
	int E_searchNumPrunedMinDEE = 0;
	float E_searchBestEnergyFound = 99999.0f;
	int EL_searchNumConfsTotal = 0;
	int EL_searchNumConfsPrunedByE = 0;
	int EL_searchNumConfsPrunedByS = 0;
	int EL_searchNumConfsEvaluated = 0;
	int EL_searchNumConfsLeft = 0;
	int EL_searchNumPrunedMinDEE = 0;
	float EL_searchBestEnergyFound = 99999.0f;
	BigDecimal q_E = new BigDecimal(0.0);
	BigDecimal q_EL = new BigDecimal(0.0);
	BigDecimal bestScore = new BigDecimal(0.0);
	double bestBoundEMin = 0.0;		// The best minimized bound energy
	double bestUnBoundEMin = 0.0;	// The best minimized unbound energy
	double bestBoundE = 0.0;		// The best bound energy (no minimization)
	double bestUnBoundE = 0.0;		// The best unbound energy (no minimization)
		// passed from the master to the slave
	
	String currentMutation[] = null;
	int residueMap[] = null;
	String resDefault[] = null;
	int rotamerIndexOffset[] = null;
	int ligNum = -1;
	String ligType = null;
	boolean ligPresent = false;
	int numMutations = 0;
	String arpFilenameMin = null;
	String arpFilenameMax = null;
	boolean minDEEtypeMS = false;
	int algOption = 0;
	int numSplits = 0;
	String AAallowed[] = null;
	boolean resMutatable[] = null;
	String minDEEfile = null;
	float initEw = 0.0f;
	float pruningE = (float)Math.pow(10,38);
	boolean E_repeatEw = false;
	boolean EL_repeatEw = false;
	boolean EL_allPruned = false;
	boolean E_allPruned = false;
	int mutationNumber = -1;
	Hashtable params = null;
	boolean PEMcomp = false;
	boolean entropyComp = false; //this *must* be false for the pairwise matrix energy computation
	boolean compASdist = false;
	boolean asDist[] = null;
	float dist = 0.0f;
	boolean compEref = false;
	
	double gamma = 0.01;
	float epsilon = 0.03f;
	int numResidues = 0;
	float stericThresh = 1.5f;
	int numInAS = 0;
	int numTotalRotamers = 152;
	int numResAllowed = 0;
	boolean isLigAA = true;
	boolean computeEVEnergy = true;
	boolean doMinimization = true;
	boolean minimizeBB = false;
	boolean repeatSearch = true;
	boolean calculateVolumes = true;
	boolean approxMinGMEC = false;
	float lambda = (float)Math.pow(10, 38);
	boolean distDepDielect = true;
	double dielectConst = 1.0;
	boolean doDihedE = false;
	boolean doSolvationE = false;
	double solvScale = 1.0;
	double vdwMult = 1.0;
	boolean scaleInt = false;
	float maxIntScale = 1.0f;
	boolean useEref = false;
	float eRef[] = null;
	
	// Timing info (in seconds)
	int q_E_Time = 0;
	int q_EL_Time = 0;

	// Identification information
	int slaveNum = -1;
		// the number of the slave with which this packet communicates 1..n
	int portNum = -1;
		// the port number at which the slave can be reached perhaps 10000...
	String machineName = null;
		// the machine name at which the slave can be reached

	// Stopping information
	boolean workToDo = true;
		// when workToDo is false the slave exits
	
	
	//Variables specific to PEM computation
	int resMut[] = null;
	String flagMutType = null;	
	SamplingEEntries compEE[] = null;//initialized by the slave node, not by the master
	int numLigRotamers = 0;	
	int elapsedTime = 0; // timing info (in seconds)
	
	//Variables specific to distributed DACS and distributed DEE computations
	boolean prunedRot[] = null;
	boolean useSF = false;
	boolean distrDACS = false;
	boolean distrDEE = false;
	boolean splitFlags[][] = null;
	String rotFileIn = null;
	String sfFileIn = null;
	String sfFileOut = null;
	int numSpPos = -1;
	int msp[] = null;
	int typeDEE = -1;
	int maxDepth = -1;
	int diffFact = -1;
	double minRatioDiff = 0.0;
	BigInteger numInitUnprunedConf = null;
	String outputPruneInfo = null;
	String outputConfInfo = null;
	int partIndex = -1;

	CommucObj() {
	}

}

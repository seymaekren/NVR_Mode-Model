////////////////////////////////////////////////////////////////////////////////////////////
// RotamerLibrary.java
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
 * This class implements a rotamer library reader
 *
 * It reads from an input file created by Ryan Lilien from the data in the 
 *   Richardsons' rotamer library articles. I created my own data file rather 
 *   than using the supplied ones as they had a very complicated and 
 *   non-standard format which would have been very difficult to parse.
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
    Duke University
	Department of Computer Science
	Levine Science Research Center (LSRC)
	Durham
	NC 27708-0129 
	USA
	brd@cs.duke.edu

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

import java.io.*;
import java.util.*;


public class RotamerLibrary implements Serializable {
	
	// If the debug flag is set to true then additional debug statements are
	//  printed to standard out.
	public static final boolean debug = false;
	
	String aaNames[];       // Array of amino-acid names
	private int numAAallowed;              // Number of AA types read
	private int numDihedrals[];     // Number of dihedrals per AA
	private int numRotamers[];      // Number of rotamers per AA
	
	private String dihedralAtomNames[][][];  // Names of atoms involved in the dihedrals for each amino acid
	private int rotamerValues[][][];  // Actual angle values for each rotamer for each amino acid
	private float rotamerVolumes[][];	// Volumes of each rotamer for each amino acid
	
	private int totalNumRotamers; //AAs with 0 rotamers are counted as 1 rotamer	
	private int rotamerIndexOffset[] = null; //the rotamer index offset for each amino acid (AAs with 0 rotamers are counted as 1 rotamer)
	
	
	// Generic constructor
	RotamerLibrary(String rotFilename, String volFilename) {		
		
		try {
			readRotLibrary(rotFilename);
			readRotVol(volFilename);
		}
		catch (Exception e){
			System.out.println("ERROR reading rotamer library/volumes: "+e);
			System.exit(1);
		}
	}
	
	//Read in all of the rotamers for all amino acids from the rotFilename file
	private void readRotLibrary(String rotFilename) throws Exception {
		
		// HANDLE THE NORMAL AAs	
		FileInputStream is = new FileInputStream( rotFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;
		int curAA = 0;
		numAAallowed = 0;
		totalNumRotamers = 0;
		// Skip over comments (lines starting with !)
		curLine = bufread.readLine();
		while( curLine.charAt(0) == '!' )
			curLine = bufread.readLine();
		
		numAAallowed = (new Integer(getToken(curLine,1))).intValue(); //the first non-comment line is the number of AA types
		curLine = bufread.readLine();
		
		aaNames = new String[numAAallowed];
		numDihedrals = new int[numAAallowed];
		numRotamers = new int[numAAallowed];
		rotamerIndexOffset = new int[numAAallowed];
		dihedralAtomNames = new String[numAAallowed][][];
		rotamerValues = new int[numAAallowed][][];
		
	  	while( curLine != null ) {
			
	  		aaNames[curAA] = getToken(curLine,1);
			numDihedrals[curAA] = (new Integer(getToken(curLine,2))).intValue();
			numRotamers[curAA] = (new Integer(getToken(curLine,3))).intValue();
			
			dihedralAtomNames[curAA] = new String[numDihedrals[curAA]][4];
			rotamerValues[curAA] = new int[numRotamers[curAA]][numDihedrals[curAA]];

			// Read in the actual dihedrals
			for(int q=0;q<numDihedrals[curAA];q++) {
				curLine = bufread.readLine();
				dihedralAtomNames[curAA][q][0] = getToken(curLine,1);
				dihedralAtomNames[curAA][q][1] = getToken(curLine,2);
				dihedralAtomNames[curAA][q][2] = getToken(curLine,3);
				dihedralAtomNames[curAA][q][3] = getToken(curLine,4);
			}
			// Read in the actual rotamers
			for(int q=0;q<numRotamers[curAA];q++) {
				curLine = bufread.readLine();
				for(int w=0;w<numDihedrals[curAA];w++) {
					rotamerValues[curAA][q][w] = (new Integer(getToken(curLine,(w+1)))).intValue();
				}
			}
			
			rotamerIndexOffset[curAA] = totalNumRotamers;
			
			totalNumRotamers += numRotamers[curAA];
			if (numRotamers[curAA]<=0) //ALA or GLY
				totalNumRotamers += 1;
			
			curAA++;
			curLine = bufread.readLine();
		}
	  	
	  	changeHIP(); //change HIS/HID/HIE to HIP in aaNames[]
	  	
		bufread.close();
		
		if (curAA!=numAAallowed){
			System.out.println("ERROR: not all amino acid types read from rotamer library");
			System.exit(1);
		}
	}
	
	//Read in all of the rotamer volumes for all amino acids from the volFilename file
	private void readRotVol(String volFilename) throws Exception {
		
		FileInputStream is = new FileInputStream( volFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null;
		
		rotamerVolumes = new float[numAAallowed][];
		
		is = new FileInputStream( volFilename );
		bufread = new BufferedReader(new InputStreamReader(is));
		
		// Skip over comments (lines starting with !)
		curLine = bufread.readLine();
		while( curLine.charAt(0) == '!' )
			curLine = bufread.readLine();
  		
		int curResult = 0;
		while( curLine != null ) {
							
			int aaIndex = getAARotamerIndex(getToken(curLine,1));
			
			int numRotVol = numRotamers[aaIndex];
			if (numRotamers[aaIndex]==0)
				numRotVol++;
			
			rotamerVolumes[aaIndex] = new float[numRotVol];
			for (int j=0; j<numRotVol; j++)
				rotamerVolumes[aaIndex][j] = new Float(getToken(curLine,j+2)).floatValue();
		
			curLine = bufread.readLine();
			curResult++;
		}
		
		bufread.close();
		
		if (curResult!=numAAallowed){
			System.out.println("ERROR: not all amino acid types read from rotamer volumes file");
			System.exit(1);
		}
	}	

	// This function returns the rotamer index for rotamers of
	//  amino acid aaName; returns -1 if name not found
	public int getAARotamerIndex(String aaName) {

		if (aaName.equalsIgnoreCase("HID") || aaName.equalsIgnoreCase("HIE") ||
			aaName.equalsIgnoreCase("HIS")) {
			aaName = "HIP";
			if(debug)
				System.out.println("ASSUMING HID/E/S is " + aaName + " for rotamer purposes.");
		}
		if (aaName.equalsIgnoreCase("CYX")) {
			aaName = "CYS";
			if(debug)
				System.out.println("ASSUMING CYX is " + aaName + " for rotamer purposes.");
		}

		for(int q=0;q<numAAallowed;q++) {
			if (aaNames[q].equalsIgnoreCase(aaName))
				return q;
		}
		return -1;
	}

	// Returns the name of the amino acid at index
	public String getAAName(int aaIndex){
		if (aaNames[aaIndex].equalsIgnoreCase("his") ||
			aaNames[aaIndex].equalsIgnoreCase("hie") ||
			aaNames[aaIndex].equalsIgnoreCase("hid") ||
			aaNames[aaIndex].equalsIgnoreCase("hip"))
			return(new String("HIP"));

		return(aaNames[aaIndex]);
	}

	// Returns the number of rotamers for the specified amino acid type
	public int getNumRotForAAtype(int aaTypeInd){
		return(numRotamers[aaTypeInd]);
	}

	// This function returns the number of rotamers for a given
	//  amino acid type (by name)
	public int getNumRotamers(String aaName) {
		int aaNum = getAARotamerIndex(aaName);
		return(numRotamers[aaNum]);
	}

	// This function returns the number of rotamers for a given
	//  amino acid type (by single char AA type)
	public int getNumRotamers(char aaName) {
		int aaNum = getAARotamerIndex(getThreeLetAACode(aaName));
		return(numRotamers[aaNum]);
	}

	// Returns the number of dihedrals for the amino acid
	//  number aaTypeInd. NOTE: aaTypeInd is not an index
	//  into a molecule, it's 0..19
	public int getNumDihedrals(int aaTypeInd){
		if (aaTypeInd != -1)
			return(numDihedrals[aaTypeInd]);
		else
			return(0);
	}


	// Returns the residue local atom numbers for dihedral dihedNum of
	//  residue resNum of strand strNum
	public int[] getDihedralInfo(Molecule m, int strNum, int resNum, int dihedNum){

		Residue localResidue = m.strand[strNum].residue[resNum];
		String tmpName = null;
		if(localResidue.name.equalsIgnoreCase("CYX"))
			tmpName = "CYS";
		else
			tmpName = localResidue.name;
				
		int aaNum = getAARotamerIndex(tmpName);
		int atNum[] = new int[4];

		// Find atoms involved in the dihedral
		for(int q=0;q<localResidue.numberOfAtoms;q++) {
			for(int w=0;w<4;w++) {
				if (localResidue.atom[q].name.equalsIgnoreCase(dihedralAtomNames[aaNum][dihedNum][w])) {
					atNum[w] = q;
				}
			}
		}
		return(atNum);
	}
	
	// This function simply converts from a 3-letter AA name
	//  to the one character AA symbol
	public char getOneCharAACode(String aaName){
	
		char oneLet = '?';
		
		if (aaName.equalsIgnoreCase("ala"))
			oneLet = 'A';
		else if (aaName.equalsIgnoreCase("cys"))
			oneLet = 'C';
		else if (aaName.equalsIgnoreCase("asp"))
			oneLet = 'D';
		else if (aaName.equalsIgnoreCase("glu"))
			oneLet = 'E';
		else if (aaName.equalsIgnoreCase("phe"))
			oneLet = 'F';
		else if (aaName.equalsIgnoreCase("gly"))
			oneLet = 'G';
		else if (aaName.equalsIgnoreCase("his"))
			oneLet = 'H';
		else if (aaName.equalsIgnoreCase("hie"))
			oneLet = 'H';
		else if (aaName.equalsIgnoreCase("hid"))
			oneLet = 'H';
		else if (aaName.equalsIgnoreCase("hip"))
			oneLet = 'H';
		else if (aaName.equalsIgnoreCase("ile"))
			oneLet = 'I';
		else if (aaName.equalsIgnoreCase("lys"))
			oneLet = 'K';
		else if (aaName.equalsIgnoreCase("leu"))
			oneLet = 'L';
		else if (aaName.equalsIgnoreCase("met"))
			oneLet = 'M';
		else if (aaName.equalsIgnoreCase("asn"))
			oneLet = 'N';
		else if (aaName.equalsIgnoreCase("pro"))
			oneLet = 'P';
		else if (aaName.equalsIgnoreCase("gln"))
			oneLet = 'Q';
		else if (aaName.equalsIgnoreCase("arg"))
			oneLet = 'R';
		else if (aaName.equalsIgnoreCase("ser"))
			oneLet = 'S';
		else if (aaName.equalsIgnoreCase("thr"))
			oneLet = 'T';
		else if (aaName.equalsIgnoreCase("val"))
			oneLet = 'V';
		else if (aaName.equalsIgnoreCase("trp"))
			oneLet = 'W';
		else if (aaName.equalsIgnoreCase("tyr"))
			oneLet = 'Y';
		else
			oneLet = '?';
		
		return oneLet;
	}

	// This function simply converts from a 3-letter AA name
	//  to the one character AA symbol
	public String getThreeLetAACode(char oneLet){
	
		String aaName = "???";
		
		if (oneLet == 'A')
			aaName = "ALA";
		else if (oneLet == 'C')
			aaName = "CYS";
		else if (oneLet == 'D')
			aaName = "ASP";
		else if (oneLet == 'E')
			aaName = "GLU";
		else if (oneLet == 'F')
			aaName = "PHE";
		else if (oneLet == 'G')
			aaName = "GLY";
		else if (oneLet == 'H')
			aaName = "HIS";
		else if (oneLet == 'I')
			aaName = "ILE";
		else if (oneLet == 'K')
			aaName = "LYS";
		else if (oneLet == 'L')
			aaName = "LEU";
		else if (oneLet == 'M')
			aaName = "MET";
		else if (oneLet == 'N')
			aaName = "ASN";
		else if (oneLet == 'P')
			aaName = "PRO";
		else if (oneLet == 'Q')
			aaName = "GLN";
		else if (oneLet == 'R')
			aaName = "ARG";
		else if (oneLet == 'S')
			aaName = "SER";
		else if (oneLet == 'T')
			aaName = "THR";
		else if (oneLet == 'V')
			aaName = "VAL";
		else if (oneLet == 'W')
			aaName = "TRP";
		else if (oneLet == 'Y')
			aaName = "TYR";

		return aaName;
	}

	// This function returns the xth token in string s
	private String getToken(String s, int x) {
	
		int curNum = 1;	
		StringTokenizer st = new StringTokenizer(s," ,;\t\n\r\f");
		
		while (curNum < x) {
			curNum++;
			if (st.hasMoreTokens())
			  st.nextToken();
			else {
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken
	
	//Change HIS/HID/HIE to HIP in aaNames[]
	private void changeHIP(){
		for (int i=0; i<aaNames.length; i++){
			if (aaNames[i].equalsIgnoreCase("his") || aaNames[i].equalsIgnoreCase("hie") ||
					aaNames[i].equalsIgnoreCase("hid") || aaNames[i].equalsIgnoreCase("hip"))
						aaNames[i] = "HIP";
		}
	}
	
	public int getNumAAallowed(){
		return numAAallowed;
	}
	
	public float [][] getRotVol(){
		return rotamerVolumes;
	}
	
	public String getDihedralAtomNames(int i, int j, int k){
		return dihedralAtomNames[i][j][k];
	}
	
	public int getRotamerValues(int i, int j, int k){
		return rotamerValues[i][j][k];
	}
	
	public int [] getRotamerIndexOffset(){
		return rotamerIndexOffset;
	}
	
	public int getTotalNumRotamers(){
		return totalNumRotamers;
	}
	
	public String [] getAAtypesAllowed(){
		return aaNames;
	}
}

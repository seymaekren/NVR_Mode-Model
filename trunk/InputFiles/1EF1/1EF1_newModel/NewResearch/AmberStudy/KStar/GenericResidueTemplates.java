////////////////////////////////////////////////////////////////////////////////////////////
// GenericResidueTemplates.java
//
//  Version:           0.1
//
//
// authors:
//    initials     name                   organization                email
//   ---------    --------------     ------------------------    ------------------------------
//     RHL         Ryan Lilien         Dartmouth College           ryan.lilien@dartmouth.edu
//	   ISG		   Ivelin Georgiev	   Duke University			   ivelin.georgiev@duke.edu     
//
////////////////////////////////////////////////////////////////////////////////////////////

/** 
 * Written by Ryan Lilien (2001-2004);
 * 	minor modifications by Ivelin Georgiev (2004-2007)
 * 
 * This class reads from a generic residue file that includes element type,
 *  AMBER atom type, limited connectivity, and partial charge
 * The format of the residue file is similar to the PARM AMBER datafiles
 * Identical to the all_amino94.in
 * These can then be used to apply these atom types and partial charges to
 *  residues of another molecule
 * This file now includes the standard AMBER nucleic acid force field parameterization (as given by all_nuc94.in)
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
 
import java.io.InputStream;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.util.*;
import java.io.Serializable;

public class GenericResidueTemplates implements Serializable{

	public static final int MAX_NUM_RES = 50;

	String grFilename = "all_nuc94_and_gr.in";

	Residue grResidues[];   // arrays of residues
	int numGRs;

	GenericResidueTemplates() throws Exception {

		FileInputStream is = new FileInputStream( grFilename );
		BufferedReader bufread = new BufferedReader(new InputStreamReader(is));
		String curLine = null, tmpName = null;
		int dumPresent = 0;
		numGRs = 0;
		grResidues = new Residue[MAX_NUM_RES];
		// Skip over first 2 lines of header info
		curLine = bufread.readLine();
		curLine = bufread.readLine();

  	while( curLine != null ) {
			// Skip over first line which is the long residue name
			curLine = bufread.readLine();
			if (curLine.length() >= 4)
				if (curLine.substring(0,4).equalsIgnoreCase("stop")) {
					curLine = bufread.readLine();
					continue;
				}
			// Skip blank line
			curLine = bufread.readLine();
			// The next line contains the 3 letter amino acid name
			curLine = bufread.readLine();
			tmpName = getToken(curLine,1);
			Residue newRes = new Residue();
			newRes.name = tmpName;
			// Skip next 2 lines
			curLine = bufread.readLine();
			curLine = bufread.readLine();
			// Now we're into the section with atoms
			curLine = bufread.readLine();
			// Skip the dummy atoms
			dumPresent = 0;
			while (getToken(curLine,2).equalsIgnoreCase("DUMM")) {
				dumPresent++;
				curLine = bufread.readLine();
			}
			dumPresent++; // to adjust for 0-based
			while (!getToken(curLine,2).equals("")) {
				Atom at = new Atom();
				tmpName = getToken(curLine,2);
				at.name = tmpName;
				tmpName = tmpName.substring(0,1);		
				at.elementType = tmpName;
				at.forceFieldType = getToken(curLine,3);
				at.charge = (float) (new Float(getToken(curLine,11)).floatValue());
				at.isBBatom = at.setIsBBatom();
				at.addBond(((new Integer(getToken(curLine,5))).intValue())-dumPresent);
				newRes.addAtom(at);  // add atom
				curLine = bufread.readLine();  // read next line
			}
			// Eventually we might want to be able to handle the improper
			//  torsions listed here

			// Add the residue to the array
			grResidues[numGRs++] = newRes;
			// Read until the end of the residue
			boolean atDone = false;
			if (curLine.length() >= 4)
				atDone = curLine.substring(0,4).equalsIgnoreCase("done");
			else
				atDone = false;
			while (!atDone) {
				curLine = bufread.readLine();
				if (curLine.length() >= 4)
					atDone = curLine.substring(0,4).equalsIgnoreCase("done");
				}
		}
		bufread.close();
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
				return(new String(""));
			}
		}

		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken

}

////////////////////////////////////////////////////////////////////////////////////////////
// SciptParser.java
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
 * This software written by Ryan Lilien (2001-2004)
 * 
 * The ScriptParser class handles reading in a script, handles loops, etc...
 * 
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

import java.io.*;
import java.util.Vector;
import java.util.Stack;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.io.PrintStream;
import java.io.BufferedInputStream;

public class ScriptParser
{

	int curLine;			// Current line being executed
	int numLines;			// Number of lines in the current program
	Vector lines;			// The actual lines of the script



	// Constructor, script filename is sname
	ScriptParser(String sname) {

	BufferedReader bufread = null;
	String s = new String("");  // line being parsed
	boolean done = false;

	// Setup and open the file
	try{
		File file = new File(sname);
		FileReader fr = new FileReader(file);
		bufread = new BufferedReader(fr);
	}
	catch (FileNotFoundException e) {
		System.out.println("ERROR: File " + sname + " Not Found");
	}

	// Read all the lines
	lines = new Vector();
	while (!done) {
		try {
			s = bufread.readLine();
			if (s == null) { // stop if we've reached EOF
				done = true;
				bufread.close();
				s = new String("");
			}
			lines.addElement(s);
		}
		catch ( Exception e ){
		 System.out.println("ERROR: An error occured while reading input");
		 System.exit(0);
		}
	}

	// Prepare the loop information
	// Each end is replaced by a line of the following form:
	//   #*# process X
	// #*# Is a marker that this is internal code
	// X is the line number of the corresponding FOR
	numLines = lines.size();
	Stack lStack = new Stack();
	Integer tmpInt = null;
	for(int i=0;i<numLines;i++) {
		s=(String)lines.elementAt(i);
		StringTokenizer st = new StringTokenizer(s," ;\t\n\r\f");
		String firstToken = new String("");
		if (st.hasMoreTokens())
			firstToken = st.nextToken();  // snag a copy of the first token
		if (firstToken.equalsIgnoreCase("for"))
			lStack.push(new Integer(i));
		if (firstToken.equalsIgnoreCase("end")) {
			if (lStack.isEmpty()) {
				System.out.println("Syntax Error: Unmatched END");
				System.exit(0);
			}
			tmpInt = (Integer)lStack.pop();
			String tmpstg = new String("#*# process " + tmpInt.intValue());
			lines.setElementAt(tmpstg,i);
		}
	}
	if (!lStack.isEmpty()) {
		System.out.println("Syntax Error: Unmatched FOR");
		System.exit(0);
	}

	} // end constructor


	// *****************************************************
	// Returns the next executable line
	// This function takes into account loops and modifies loop 
	//   variables appropriately
	public String getNextLine(Hashtable curScope) {

	if (curLine >= numLines)
	 return(new String("EOF"));
	 
	String curStr = (String)lines.elementAt(curLine++);
	ScopeObject so;
	String retStr = curStr;

	StringTokenizer st = new StringTokenizer(curStr," ;\t\n\r\f");
	String firstToken = new String("");
	if (st.hasMoreTokens())
		firstToken = st.nextToken();  // snag a copy of the first token
	if (firstToken.equalsIgnoreCase("for")) {
		if (numTokens(curStr) > 5) {
			// Make new loop variable (smashes old declaration)
			if (getToken(curStr,2).equalsIgnoreCase("int")) {
				// Add int object to scope
				so = new ScopeObject(new Integer(parseToInt(getToken(curStr,4),curScope)),getToken(curStr,3),11);
				curScope.put(getToken(curStr,3),so);
			}
			else if (getToken(curStr,2).equalsIgnoreCase("float") || getToken(curStr,2).equalsIgnoreCase("double")) {
				// Add double object to scope
				so = new ScopeObject(new Double(parseToDouble(getToken(curStr,4),curScope)),getToken(curStr,3),12);
				curScope.put(getToken(curStr,3),so);
			}
		retStr = getNextLine(curScope);		
		}
		else {
			System.out.println("Syntax Error: Not enough parameters for FOR");
			System.exit(0);		
		}
	}
	else if (firstToken.equalsIgnoreCase("#*#")) {
		if (getToken(curStr,2).equalsIgnoreCase("process")) {
			int targetLine = parseToInt(getToken(curStr,3),curScope);
			String tmpStr = (String)lines.elementAt(targetLine);
			if (getToken(tmpStr,2).equalsIgnoreCase("int")) {
				// Get int variable object from scope
				so = (ScopeObject) curScope.get(getToken(tmpStr,3));
				Integer tmpInt = (Integer)so.data;
				int lcvValue = tmpInt.intValue();
				// Get increment value from line
				int incValue = parseToInt(getToken(tmpStr,6),curScope);
				// Get max value from line
				int maxValue = parseToInt(getToken(tmpStr,5),curScope);
				lcvValue += incValue;
				if (lcvValue<=maxValue) {
					// Put int object back in scope
					so = new ScopeObject(new Integer(lcvValue),getToken(tmpStr,3),11);
					curScope.put(getToken(tmpStr,3),so);
					curLine = targetLine+1;
				}
				retStr = getNextLine(curScope);
			}
			else if (getToken(tmpStr,2).equalsIgnoreCase("float") || getToken(tmpStr,2).equalsIgnoreCase("double")) {
				// Get double variable object from scope
				so = (ScopeObject) curScope.get(getToken(tmpStr,3));
				Double tmpDoub = (Double)so.data;
				double lcvValue = tmpDoub.doubleValue();
				// Get increment value from line
				double incValue = parseToDouble(getToken(tmpStr,6),curScope);
				// Get max value from line
				double maxValue = parseToDouble(getToken(tmpStr,5),curScope);
				lcvValue += incValue;
				if (lcvValue<=maxValue) {
					// Put double object back in scope
					so = new ScopeObject(new Double(lcvValue),getToken(tmpStr,3),12);
					curScope.put(getToken(tmpStr,3),so);
					curLine = targetLine+1;
				}
				retStr = getNextLine(curScope);
			}
		}
	}

	return(retStr);

	} // end getNextLine()


	/******************************/
	// This function returns true if the variable named
	//  s is in the current scope
	private boolean inScope(String s, Hashtable curScope) {
		
		if (curScope.get(s) != null)
			return true;
		return false;
	} // end inScope()

	/******************************/
	// These functions convert between various number types
	private float double2float(double d) {
		return( (new Float(d)).floatValue() );
	}
		
	private int float2int (float f) {
		return( (new Float(f)).intValue() );
	}

	private int double2int (double d) {
		return( (new Double(d)).intValue() );
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
		
		Boolean tmpBool;
		boolean d = false;
		ScopeObject so = null;
			
		if (inScope(s,curScope)) {
			so = (ScopeObject) curScope.get(s);
			if (so == null) {
				System.out.println("ERROR: undefined variable " + getToken(s,1));
				return false;
			}
			if (so.objType != 13) {
				System.out.println("ERROR: variable " + s + " is not of type boolean");
				return false;
			}
			tmpBool = (Boolean) so.data;	
			d = tmpBool.booleanValue();
		}
		else
			d = new Boolean(s).booleanValue();
		
		return d;
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
				System.out.println("ERROR: Unable to access parameter " + x);
				return(new String(""));
			}
		}
			
		if (st.hasMoreTokens())		
			return(st.nextToken());
		return(new String(""));

	} // end getToken


} // end class ScriptParser
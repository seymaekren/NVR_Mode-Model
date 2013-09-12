///////////////////////////////////////////////////////////////////////////////////////////////
//	EEF1.java
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
 * Manages the computation of exp(x) for large values of x, using BigDecimal;
 * 		For large values of x, the standard Math.exp(x) cannot be used, since it only has double precision;
 * Implements pow() for integer powers of a BigDecimal number and an approximation to the natural logarithm of a BigDecimal number
 * 
*/

import java.math.*;

public class ExpFunction {
	
	BigDecimal exp = new BigDecimal("2.71828182845904523536"); //Euler's number to 20 decimal digits
	
	final int maxPrecision = 8; //the number of decimal digits to which the BigDecimal numbers must be accurate

	//constructor
	ExpFunction(){
	}
	
	//Computes exp(x) using BigDecimal arithmetic for large x or the standard Math.exp() function for small x;
	//		If x is large, it is divided into its integer and fractional part; BigDecimal is used to compute
	//			exp to the power of the integer-part of x; the standard Math.exp() function is used to compute
	//			exp(fractional power), since the fractional part < 1; the two results are then multiplied to
	//			obtain the BigDecimal exp(x)
	public BigDecimal exp(double x){
		
		final double minX = 0.0; //only use BigDecimal arithmetic for (x>=minX)
		
		BigDecimal expX = null;
		
		if (x<minX){ //x is small enough, so use the standard Math.exp() function
			expX = new BigDecimal(Math.exp(x));
		}
		else { //x is large, so use the BigDecimal arithmetic approach
			int integerPart = (int)Math.floor(x);
			double fractionPart = x - integerPart;
			
			BigDecimal intPartExp = pow(exp,integerPart);
			BigDecimal fractPartExp = new BigDecimal(Math.exp(fractionPart));
			
			expX = intPartExp.multiply(fractPartExp);
		}
		
		expX = expX.setScale(maxPrecision,4); //rounding is ROUND_HALF_UP (standard rounding: up for next digit >=5, down otherwise)
		
		return expX;
	}
	
	//Returns the BigDecimal number num to the power a
	public BigDecimal pow(BigDecimal num, int a){
		
		BigDecimal expPow = new BigDecimal("1.0");
		for (int i=0; i<a; i++)
			expPow = expPow.multiply(num);
		
		return expPow;
	}
	
	//Returns an approximation to the natural logarithm of the BigDecimal number num
	public BigDecimal log(BigDecimal num){
		
		final BigDecimal ONE = new BigDecimal("1.0");
		
		BigDecimal sum = new BigDecimal("0.0");
		BigDecimal x = num;
		
		if (num.compareTo(new BigDecimal(Math.pow(10, 38)))<0){ //num is small, so use the standard Math.log() function
			if (num.compareTo(new BigDecimal("0.00001"))<0)
				sum = new BigDecimal("0.0");
			else
				sum = new BigDecimal(Math.log(num.doubleValue()));
		}
		else { //num is large, so compute an approximation to the natural logarithm
			
			double t = 0.0;
			
			boolean done = false;
			while (!done){
				if (x.compareTo(exp)>0){
					t += 1.0;
				}
				else {
					sum = sum.add(new BigDecimal(t+Math.log(x.doubleValue())));
					done = true;
				}
				x = x.divide(exp,4);
			}
		}
		
		return sum;
	}
}

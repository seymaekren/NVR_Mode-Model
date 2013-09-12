////////////////////////////////////////////////////////////////////////////////////////////
// RotMatrix.java
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

/** 
 * This class implements rotation matricies.
 * Written by: Ryan Lilien  (2001-2004)
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

public class RotMatrix
{


	RotMatrix()
	{
	
	}

	// Rotates an array of points around the +x axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void xAxisRotate(float thetaDeg, float theCoords[],int numCoords) {

		axisRotate(1.0f,0.0f,0.0f,thetaDeg,theCoords,numCoords);
	}

	// Rotates an array of points around the +x axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void xAxisRotate(double thetaDeg, float theCoords[],int numCoords) {

		axisRotate(1.0f,0.0f,0.0f,new Double(thetaDeg).floatValue(),theCoords,numCoords);
	}

	// Rotates an array of points around the +y axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void yAxisRotate(float thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,1.0f,0.0f,thetaDeg,theCoords,numCoords);
	}

	// Rotates an array of points around the +y axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void yAxisRotate(double thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,1.0f,0.0f,new Double(thetaDeg).floatValue(),theCoords,numCoords);
	}

	// Rotates an array of points around the +z axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void zAxisRotate(float thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,0.0f,1.0f,thetaDeg,theCoords,numCoords);
	}

	// Rotates an array of points around the +z axis
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void zAxisRotate(double thetaDeg, float theCoords[],int numCoords) {

		axisRotate(0.0f,0.0f,1.0f,new Double(thetaDeg).floatValue(),theCoords,numCoords);
	}


	// Rotates an array of points around the axis ax, ay, az
	//  using a right handed rotation of thetaDeg degrees
	// theCoords are the coordinates where the ith coodinates
	//  are in position 3*i .. (3*(i+1))-1
	public void axisRotate(float ax, float ay, float az, float thetaDeg,
	 float theCoords[],int numCoords) {
		
			float tx,ty,tz;

			float[][] rot_mtx = new float[3][3];
			getRotMatrix(ax,ay,az,(float) thetaDeg,rot_mtx);

			int ix3 = 0;
			for(int i=0;i<numCoords;i++){
				tx=theCoords[ix3];
				ty=theCoords[ix3+1];
				tz=theCoords[ix3+2];
			
				theCoords[ix3] = tx * rot_mtx[0][0] + ty * rot_mtx[0][1] + tz * rot_mtx[0][2];
				theCoords[ix3+1] = tx * rot_mtx[1][0] + ty * rot_mtx[1][1] + tz * rot_mtx[1][2];
				theCoords[ix3+2] = tx * rot_mtx[2][0] + ty * rot_mtx[2][1] + tz * rot_mtx[2][2];
			
				ix3 += 3;
			}
	}	
	
	// Translates an array of points by the specified amount
	public void translate(float tx, float ty, float tz, float theCoords[],
	 int numCoords) {

		int ix3 = 0;
		for(int i=0;i<numCoords;i++){
			theCoords[ix3] += tx;
			theCoords[ix3+1] += ty;
			theCoords[ix3+2] += tz;
			ix3 += 3;
		}
	}

	// This function constructs a rotation matrix from a rotation in
	//  axis-angle notation
	public void getRotMatrix(float fx, float fy, float fz, float angle,
		float[][] rot_mtx) {

		// First convert the axisangle to a quaternion
		float sin_a = (float) Math.sin(angle*3.14159265/180/2);
		float cos_a = (float) Math.cos(angle*3.14159265/180/2);
		float tmp = (float) Math.sqrt(fx*fx + fy*fy + fz*fz);
		float qx = fx / tmp * sin_a;
		float qy = fy / tmp * sin_a;
		float qz = fz / tmp * sin_a;	
		float qw = cos_a;
		tmp = (float) Math.sqrt(qx*qx + qy*qy + qz*qz + qw*qw);
		qx /= tmp;
		qy /= tmp;
		qz /= tmp;
		qw /= tmp;
		float xx = qx * qx;
		float xy = qx * qy;
		float xz = qx * qz;
		float xw = qx * qw;

		float yy = qy * qy;
		float yz = qy * qz;
		float yw = qy * qw;

		float zz = qz * qz;
		float zw = qz * qw;
			
		rot_mtx[0][0] = 1 - 2 * (yy + zz);
		rot_mtx[0][1] = 2 * (xy - zw);
		rot_mtx[0][2] = 2 * (xz + yw);

		rot_mtx[1][0] = 2 * (xy + zw);
		rot_mtx[1][1] = 1 - 2 * (xx + zz);
		rot_mtx[1][2] = 2 * (yz - xw);

		rot_mtx[2][0] = 2 * (xz - yw);
		rot_mtx[2][1] = 2 * (yz + xw);
		rot_mtx[2][2] = 1 - 2 * (xx + yy);
	}

}

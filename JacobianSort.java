package Jacobian;

import Jama.*;

import java.util.ArrayList;

/**
 * @author Anthony Stange
 */

public class JacobianSort {
    
    //Returns offB values of inputted matrix
    public static ArrayList<Double> getOffBs(Matrix m) {
        ArrayList<Double> values = new ArrayList<Double>();
        while (getOffB(m) > .000000001 ) {
            values.add(getOffB(m));
            m = gTAG(m);
        }       
        return values;
    }
    
    
    //returns offA value
    public static double offA(Matrix orig) {
        return getOffB(orig);
    }

    
    //finds the max by checking the upper right diagonal
    private static int[] getMax(Matrix matrix) {
    	
		int[] max = {0,1};
		
		for (int i = 0; i < matrix.getRowDimension(); i++) {
			
		    for (int j = i; j < matrix.getColumnDimension(); j++) {
		    	
	                if (Math.abs(matrix.get(i,j)) > Math.abs(matrix.get(max[0], max[1]))) {
	                	
			            if (i != j) {
				            max[0] = i;
				            max[1] = j;
			            }
	                }
	            }
	        }
		return max;
    }



    //makes the 2 by 2 matrix by using the upper right diagonal method
    private static Matrix get2X2(Matrix matrix) {
		Matrix returnMatrix = new Matrix(2,2);
		int[] max = getMax(matrix);
		returnMatrix.set(0,1,matrix.get(max[0],max[1]));
		returnMatrix.set(1,0,matrix.get(max[1],max[0]));
		returnMatrix.set(0,0,matrix.get(max[0],max[0]));
		returnMatrix.set(1,1,matrix.get(max[1],max[1]));
	
		return returnMatrix;
    }

    //returns Matrix V
    private static Matrix getEigenVectors2x2(Matrix matrix) {
    	//calculates M+ using (a+d)/2 + (b^2 + (a-d)/2)^(1/2)
    	double mPlus = (((matrix.get(0, 0) + matrix.get(1, 1))/2 + 
    			Math.sqrt(Math.pow(matrix.get(1, 0), 2) + Math.pow(((matrix.get(0, 0) - matrix.get(1, 1))/2), 2))));
    	//gets the top left hand corner and does a-eigenvalue(I)
    	double topLeft = matrix.get(0, 0);
    	topLeft = topLeft - mPlus;
    	
    	//produces u1
    	double[][] u1 = new double[1][2];
    	u1[0][0] = matrix.get(0, 1);
    	u1[0][1] = topLeft;
    	
    	//  u1/|u1|
    	double u1Size = Math.sqrt(Math.pow(u1[0][0], 2) + Math.pow(u1[0][1], 2));
    	u1[0][0] = u1[0][0]/u1Size;
    	u1[0][1] = u1[0][1]/u1Size;
    	
    	// u2 = u1 orthognal
    	double[][] u2  = new double[1][2];
    	u2[0][0] = -u1[0][1];
    	u2[0][1] = u1[0][0];
    	
    	// forms V
    	double[][] V = new double[2][2];
    	V[0][0] = u1[0][0];
    	V[1][0] = u2[0][0];
    	V[0][1] = u1[0][1];
    	V[1][1] = u2[0][1];
    	
    	//creates matrix and returns it
    	Matrix eigVectMatrix = new Matrix(V,2,2);
		return eigVectMatrix;
    }


    //Computes Givens matrix and returns it
    private static Matrix getGivens(Matrix matrix) {
    	
		Matrix identity = Matrix.identity(matrix.getRowDimension(),matrix.getColumnDimension());
		
		Matrix twoByTwo = get2X2(matrix);
		Matrix subValues = getEigenVectors2x2(twoByTwo);
		int[] max = getMax(matrix);
		identity.set(max[0],max[0], subValues.get(0,0));
		identity.set(max[0],max[1], subValues.get(0,1));
		identity.set(max[1],max[0], subValues.get(1,0));
		identity.set(max[1],max[1], subValues.get(1,1));
		return identity;
    }
 

    //Computes G^t * A * G and returns it
    private static Matrix gTAG(Matrix matrix) {
		Matrix g = getGivens(matrix);
		Matrix gTranspose = g.transpose();
		Matrix result = gTranspose.times(matrix);
		result = result.times(g);
		return result;
    }


    //Calculates the offB value of the matrix and returns it
    private static double getOffB(Matrix matrix) {
		double offVals = 0;
		for (int i = 0; i < matrix.getRowDimension(); i++) {
		    for (int j = 0; j < matrix.getColumnDimension(); j++) {
				if (i != j) {
				    offVals += Math.pow(matrix.get(i,j),2);
				}
		    }
		}
		return offVals;
    }
}

package Jacobian;

import Jama.*;

import java.util.ArrayList;

/**
 * @author Anthony Stange
 */
public class JacobianNoSort {
    
	//locations of upper right values so that the no sorting alg. can sweep through them
	//when it is getting the off b values.
    private static int[][] coordinates = {{0,1},{0,2},{0,3},{0,4},{1,2},{1,3},{1,4},{2,3},{2,4},{3,4}};
    private static int index;
        
    //gets the OffB values for the inputted matrix
    public static ArrayList<Double> getOffBs(Matrix m) {
    	index = 0;
        ArrayList<Double> values = new ArrayList<Double>();
        
        //will keep running until the value is less than 10^-9.
        while (getOffB(m) > .000000001 ) {
            values.add(getOffB(m));
            m = gTAG(m);
            index++;
            if (index >= coordinates.length){
                index = 0;
            }
        }
        
        return values;
    }
    
    
    //returns offA of the matrix
    public static double offA(Matrix orig) {
        return getOffB(orig);
    }

    
    //finds max off diagonal value for the matrix so it can get a givens.
    private static int[] getMaxDiagonal(Matrix matrix) {
		int[] max = {0,1};
		
		//for loops iterate through the matrix
		for (int i = 0; i < matrix.getRowDimension(); i++) {
		    for (int j = i; j < matrix.getColumnDimension(); j++) {
		    	
		    		//if the current i, j absolute value is greater than the previous i, j then it is maximum
	                if (Math.abs(matrix.get(i,j)) > Math.abs(matrix.get(max[0], max[1]))) {
	                	
	                	//only if the i and j are not equal
			            if (i != j) {
				            max[0] = i;
				            max[1] = j;
			            }
			        }
	            }
	        }
		return max;
    }

    //using the coordinates array we create a matrix from the max diagonal specified in the array
    private static Matrix get2X2(Matrix matrix) {
		Matrix returnMatrix = new Matrix(2,2);
		
        int[] max = coordinates[index];
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

    //Returns the Givens matrix for the inputted matrix 
    private static Matrix getGivens(Matrix matrix) {
    	
    	//Makes Identity matrix that is the same size as the inputted matrix
		Matrix identity = Matrix.identity(matrix.getRowDimension(),matrix.getColumnDimension());
		
		Matrix twoByTwo = get2X2(matrix);
		
		//Retrieves matrix U from the created 2 by 2 matrix
		Matrix subValues = getEigenVectors2x2(twoByTwo);
		int[] max = getMaxDiagonal(matrix);
		identity.set(max[0],max[0], subValues.get(0,0));
		identity.set(max[0],max[1], subValues.get(0,1));
		identity.set(max[1],max[0], subValues.get(1,0));
		identity.set(max[1],max[1], subValues.get(1,1));
		return identity;
    }
 

    //Returns G^t * A * G
    private static Matrix gTAG(Matrix matrix) {
		Matrix g = getGivens(matrix);
		Matrix gTranspose = g.transpose();
		Matrix result = gTranspose.times(matrix);
		result = result.times(g);
		return result;
    }


    //retrieves the offB values for the inputted matrix
    private static double getOffB(Matrix matrix) {
		double offVals = 0;
		for (int i = 0; i < matrix.getRowDimension(); i++) {
			    for (int j = 0; j < matrix.getColumnDimension(); j++) {
					if (i != j) {
						//Off(B) = |Bi,j|^2
					    offVals += Math.pow(matrix.get(i,j),2);
					}
			    }
		}
		return offVals;
    }
}

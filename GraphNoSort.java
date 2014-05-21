package Jacobian;

import java.awt.*;
import javax.swing.*;

import Jama.Matrix;

import java.util.ArrayList;

/**
 *
 * @author Anthony Stange
 */
public class GraphNoSort extends JPanel{
    
	private static final long serialVersionUID = -7365703762651574779L;
	private ArrayList<Double> bK;
    private int run;
    private String output;
    private Matrix m, orig;
    
    /**
     * Creates the Panel to store the graph
     * and sets run to 0 and the current matrix to 1
     */
    public GraphNoSort() {
        setPreferredSize(new Dimension(700,300));
        setBackground(Color.white);
        run = 0;
    }
    
    /**
     * When my driver tells the program to run the method generates
     * 10 random matrices and runs the jacobi no sort algorithm
     * on them
     */
    public void run(){
    	//resets the output string
    	output = "";
    	
    	//resets the bK variable
    	bK = new ArrayList<Double>();
    	
        output += "OffB values for No Sort Method: " + "\n";
        
    	//Creates a random matrix
    	m = Matrix.random(5, 5);
    	
    	//makes the matrix symmetric
    	m = m.transpose().times(m);
    	
    	//puts the untouched matrix in the original array that contains
    	//the originals of all the matrices
    	orig = m;
    	
    	//Performs the offB diagonalization on the matrix
        ArrayList<Double> list = JacobianNoSort.getOffBs(m);
        
        //adds the  ln(offBs) to the bK array list
        for (double item : list) {
            bK.add(Math.log(item));
        }
        
        //Also puts all the offBs in a string format so we can see them
        //printed on the screen
        output += "\n";
        for(double item: list){
        	output += item + "\n";	
        }
        output += "\n";
        
        //increments run so the graph will display the bk values
        run++;
    }
    
    //returns the offBs for all ten matrices
    public String getOffB(){
		return output;
    }
    
    //If run() has been run the graph will display the boundary
    //in blue and the matrix in black
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        if(run == 1){
            g.setColor(Color.blue);
            g.drawLine(0,150,700,150);
            g.drawLine(0,(int)(-Math.log(JacobianNoSort.offA(orig)) + 150),700, 75 - (int)Math.log(JacobianNoSort.offA(orig))+150); 
            
            g.setColor(Color.BLACK);
            for (int i = 0; i < bK.size(); i++) {
                g.drawOval(i * 5, (int) (-bK.get(i) + 150),2,2);
            }   
        } 
        //resets variables
        run = 0;
        output = "";
    }
}

package Jacobian;

import java.awt.*;

import javax.swing.*;

import Jama.Matrix;

import java.util.ArrayList;

/**
 * @author Anthony Stange
 */
public class GraphSort extends JPanel{

	private static final long serialVersionUID = -5523188430116811298L;
	private ArrayList<Double> bK;
    private int run;
    private Matrix orig, m;
    private String output;
    
    /**
     * Creates JPanel for graph and initializes variables
     */
    public GraphSort() {
        setPreferredSize(new Dimension(300,300));
        setBackground(Color.white);
        run = 0;
        output = "";
    }
    
    public void run(){
    	//resets possibly previously used variables
    	output = "";
    	bK = new ArrayList<Double>();        
        
        output += "OffB values for Sort Method: " + "\n";
        	
    	//Creates a random matrix and makes it symmetric
    	m = Matrix.random(5, 5);
    	m = m.transpose().times(m);
    	
    	//puts untouched matrix into an original matrix so that it
    	//can be used to get the off a values
    	orig = m;
    	
    	//Performs sorted offB diagonalization for matrix
        ArrayList<Double> list = JacobianSort.getOffBs(m);
        for (double item : list) {
            bK.add(Math.log(item));
        }

        output += "\n";
        for(double item: list){
        	output += item + "\n";	
        }
        output += "\n";
        
        run++;
    }
    
    //returns offB values in an array
    public String getOffB(){
    	return output;
    }
    
    //If run() has been run the graph will display the boundary
    //in blue and the matrix offB values in black
    public void paintComponent(Graphics g) {
        super.paintComponent(g);
        if(run == 1){
            g.setColor(Color.blue);
            g.drawLine(0,150,300,150);
            g.drawLine(0,(int)(-Math.log(JacobianSort.offA(orig)) + 150),300, 32 - (int)Math.log(JacobianSort.offA(orig))+150); 

            g.setColor(Color.BLACK);
            for (int i = 0; i < bK.size(); i++) {
                g.drawOval((i * 5), (int) (-bK.get(i) + 150),2,2);
            }   
            run = 0;
            output = "";
        } 
    }

}

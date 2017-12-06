package drawing;

import java.awt.Color;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import util.ArrayOps;
import util.NumFormat;
import util.Printer;
import util.geom.Distance;
import util.regression.ACM_CustomFunctions;
import util.regression.ACM_NonLinearFitter;
import util.regression.ACM_NonLinearFitter.FittingResult;

public class LayerCutDrawer {
	
	
	GraphDrawerCart g;
//	Layer t = null;
	LayerViewer lv = null;
	int npts;
	JFileChooser fc;
	
	double[][] points; //[n][2]
	double[] s; //distance along the line
	double[] spec; //[nspec]
	boolean drawThin = true;
	
	double[] sVert = new double [2], yVert = new double [2]; //This is for drawing a vertical line. sVert will be the position of
	//the cursor in a LayerViewer along the line
	double min, max;
	boolean drawingVert = true;
	Color vertC = Color.green;
	
	boolean drawingLines = true; //else, render the line cut as an image
	int fx, fy; //image size ratios
	
	private double[] tempx, tempy;
	
	Color def = Color.RED;

	
	boolean evaluateBiCubic = false;
	
	//refresh the number of points
	boolean refreshNPTS = false;
	double spacing = 1;
	
	//do fitting
	boolean doingFitting = false;
	FittingResult fit = null;
	String function = "FermiLine";
	Color fitC = Color.blue;
	
	public LayerCutDrawer(LayerViewer lv, int nspec, double[][] line)
	{
		this.lv = lv;
		this.npts = nspec;
		points = new double [nspec][2];
		tempx = new double [nspec];
		tempy = new double [nspec];
		s = new double [npts];
		putPoints(line);
//		nspec = points.length;
		spec = new double [nspec];
		
		evaluateSpec();
		g = new GraphDrawerCart("Line Cut Viewer", true);
		setUpDrawer();
		g.showWindow();
	}
	
	public void setNpts(int npts)
	{
		double[][] line = new double [2][2];
		line[0] = points[0];
		line[1] = points[this.npts-1];
		
		this.npts = npts;
		tempx = new double [npts];
		tempy = new double [npts];
		points = new double [npts][2];
		s = new double [npts];
		putPoints(line);
		spec = new double [npts];
		evaluateSpec();
		setUpDrawer();
	}
	public void evaluateSpec()
	{
		if (!evaluateBiCubic)
			for (int i = 0; i < npts; i++)
					spec[i] = lv.t.evaluateAt(points[i]);
		else
			for (int i = 0; i < npts; i++)
				spec[i] = lv.t.evaluateBiCubic(points[i]);
			
		
		if (g != null){
			//The x and y data of the line cut
			g.setDomain(0, s[s.length-1]);
			g.setRange(spec);
		}

		if (doingFitting)
		{
			fit = ACM_NonLinearFitter.fitToFunctionFull(s, spec, function);
			g.setXY2(fit.x, fit.yCalc, true);
			g.setColor(fitC, 2);
		}
		if (drawingVert)
		{
			yVert[0] = ArrayOps.min(spec);
			yVert[1] = ArrayOps.max(spec);
		}
	}
	public void setUpDrawer()
	{
		max = ArrayOps.max(spec); min = ArrayOps.min(spec);
		yVert[0] = min;
		yVert[1] = max;
		g.setNumPlots(npts);
		//Draws the line cut spectra
		g.setXY(new GraphDrawerCart.GraphObject(s, spec), false, false, 0);
		g.setColor(def, 0);
		g.setThin(drawThin, 0);
		//
		if (drawingVert)
		{
			g.setXY(new GraphDrawerCart.GraphObject(sVert, yVert), false, false, 1);
			g.setColor(vertC, 1);
			g.setThin(true, 1);
		}
		g.setDomain(0, s[s.length-1]);
		g.setRange(min, max);
		g.repaint();
	}
	
	public void refresh()
	{
		evaluateSpec();
		String message = "Line from " + Printer.vectorP(points[0]) + " to " + Printer.vectorP(points[points.length-1]) + 
				"; " + npts + " points; ";
		if (doingFitting)
		{
			String[] paramNames =  ACM_CustomFunctions.getNew(function).plist;
			
			for (int i = 0; i < paramNames.length; i++)
				message += paramNames[i] + ": " + NumFormat.scientific(fit.fitParams[i], 3) + "; ";
		
		}
		g.setTitle(message);
		
		g.repaint();
	}
	public void dispose()
	{
		g.dispose();
	}
	public void putPoints(double[][] line)
	{
		
		double dist = Distance.distance(lv.t.getMetricCoords(line[0]), lv.t.getMetricCoords(line[1]));
		
		if (refreshNPTS)
		{
			points[0][0] = line[0][0];
			points[0][1] = line[0][1];
			points[points.length-1][0] = line[1][0];
			points[points.length-1][1] = line[1][1];
			int n = (int)(dist/spacing);
			
			if (n != npts){
				setNpts(n);
				return;
			}
			
		}
		//System.out.println(npts);
	
		ArrayOps.putArrayInclBoth(line[0][0], line[1][0], npts, tempx);
		ArrayOps.putArrayInclBoth(line[0][1], line[1][1], npts, tempy);
		ArrayOps.putArrayInclBoth(0, dist, npts, s);
		for (int i = 0; i < npts; i++)
		{
			points[i][0] = tempx[i];
			points[i][1] = tempy[i];
		}
		
	}

	public void processKeyStroke(char ch) {
		// TODO Auto-generated method stub
		if (ch == 'n') {
			this.setNpts(Integer.parseInt(JOptionPane.showInputDialog("Enter the desired number of points")));
		}
		if (ch == 'T'){
			drawingVert = !drawingVert;
			g.plot[1].isDrawing = drawingVert;
		}
		if (ch == 'i'){
			drawThin = !drawThin;
			for (int i = 0; i < npts; i++)
			{
				g.setThin(drawThin, i);
			}
			g.repaint();
		}
		if (ch == 'v') //print x and y to the console
		{
//			double[] temp;
			for (int i = 0; i < npts; i++)
				System.out.println(s[i] + "\t" + spec[i]);
		}
		if (ch == 'V') //save everything
		{
			double[] x = new double [npts], y = new double [npts];
			double[] r = new double [2];
		}
		if (ch == 'U')
		{
			evaluateBiCubic = !evaluateBiCubic;
		}
		if (ch == 'D' && refreshNPTS)
		{
			spacing = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired constant spacing between points (pixel units)"));
		}
		if (ch == 'R')
		{
			refreshNPTS = !refreshNPTS;
			evaluateSpec();
		}
		if (ch == 'F')
		{
			doingFitting = !doingFitting;
			if (doingFitting)
			{
				g.setColor(fitC, 2);
				evaluateSpec();
			}
			else
			{
				g.plot[2].isDrawing = false;
				evaluateSpec();
			}
		}
	}
	public void setFC(JFileChooser fc)
	{
		this.fc = fc;
	}
	public void setSVert(double s)
	{
		sVert[0] = s;
		sVert[1] = s;
		
		if (drawingVert) g.repaint();

	}
}

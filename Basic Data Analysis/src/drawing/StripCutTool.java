package drawing;

import java.awt.Color;

import javax.swing.JOptionPane;

import util.ArrayOps;
import util.calc.StripCut;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.Topomap;

public class StripCutTool {
	
	
	GraphDrawerCart g;
	Topomap t = null;
	int nspec;
	
	double dh; //height offset.
	double[] x; //energy numbers.
	double[][] points; //[n][2]
	StripCut[] cut;
	//	double[][] spec; //[t.nlayers][nspec]. This is because it is a stripcut (I vs k plot);
	double[][] offsetSpec;
	boolean drawThin = true;
	
	private double[] tempx, tempy;
	
	Color def = Color.RED;
	
	public StripCutTool(Topomap t, int nspec, double r0, double r1, double width)
	{
		this.t = t;
		x = ArrayOps.generateArrayInclBoth(r0, r1, nspec);
		this.nspec = nspec;
		points = new double [nspec][2];
		tempx = new double [nspec];
		tempy = new double [nspec];
		cut = new StripCut[t.nlayers];
		for (int i = 0; i < cut.length; i++)
			cut[i] = new StripCut(t.data[i], 0, r0, r1, width, nspec);
		offsetSpec = new double [t.nlayers][nspec];
		evaluateSpec();
		setDefaultDH();
		setOffsetSpec();
		g = new GraphDrawerCart("Line Cut Viewer", true);
		setUpDrawer();
		g.showWindow();
	}
	public StripCutTool(Layer l, int nspec, double r0, double r1, double width)
	{
		this.t = Topomap.newTopomap(new Layer[] {l});
		x = ArrayOps.generateArrayInclBoth(r0, r1, nspec);
		this.nspec = nspec;
		points = new double [nspec][2];
		tempx = new double [nspec];
		tempy = new double [nspec];
		cut = new StripCut[t.nlayers];
		for (int i = 0; i < cut.length; i++)
			cut[i] = new StripCut(t.data[i], 0, r0, r1, width, nspec);
		offsetSpec = new double [t.nlayers][nspec];
		
		evaluateSpec();
		setDefaultDH();
		setOffsetSpec();
		g = new GraphDrawerCart("Line Cut Viewer", true);
		setUpDrawer();
		g.showWindow();
	}
	
	private void setOffsetSpec() {
		// TODO Auto-generated method stub
		for (int i = 0; i < nspec; i++)
			for (int j = 0; j < t.nlayers; j++)
				offsetSpec[j][i] = cut[j].result[i] + dh*j;
	}
	public void changeNspec(int nspec)
	{
		double[][] line = new double [2][2];
		line[0] = points[0];
		line[1] = points[this.nspec-1];
		this.nspec = nspec;
		tempx = new double [nspec];
		tempy = new double [nspec];
		points = new double [nspec][2];
		x = ArrayOps.generateArrayInclBoth(cut[0].r0, cut[0].r1, nspec);
		for (int i = 0; i < cut.length; i++)
			cut[i].setNPTS(nspec);
			offsetSpec = new double [t.nlayers][nspec];
		evaluateSpec();
		this.setDefaultDH();
		setOffsetSpec();
		setUpDrawer();
	}
	public void setDefaultDH()
	{
		dh = (ArrayOps.max(cut[0].result) - ArrayOps.min(cut[0].result))/nspec;
	}
	public void evaluateSpec()
	{
			for (int j = 0; j < t.nlayers; j++)
				cut[j].makeCut();
	}
	public void setUpDrawer()
	{
		double max = ArrayOps.max(offsetSpec), min = ArrayOps.min(offsetSpec);
		
		g.setNumPlots(t.nlayers);
		for (int i = 0; i < t.nlayers; i++)
		{
			g.setXY(new GraphDrawerCart.GraphObject(x, offsetSpec[i]), false, false, i);
			g.setColor(def, i);
			g.setThin(drawThin, i);
		}
		g.setDomain(x[0], x[x.length-1]);
		g.setRange(min, max);
		g.repaint();
	}
	
	public void refresh()
	{
		evaluateSpec();
		setOffsetSpec();
		g.repaint();
	}
	public void dispose()
	{
		g.dispose();
	}
	public void setR0(double r0)
	{
		for (int i = 0; i < t.nlayers; i++)
			cut[i].setR0(r0);
		x = ArrayOps.generateArrayInclBoth(cut[0].r0, cut[0].r1, nspec);
	}
	public void setR1(double r1)
	{
		for (int i = 0; i < t.nlayers; i++)
			cut[i].setR1(r1);
		x = ArrayOps.generateArrayInclBoth(cut[0].r0, cut[0].r1, nspec);
	}
	public void setPhi(double phi) {
		// TODO Auto-generated method stub
		for (int i = 0; i < t.nlayers; i++)
			cut[i].setPhi(phi);
	}
	public void setWidth(double width)
	{
		for (int i = 0; i < t.nlayers; i++)
			cut[i].setWidth(width);
	}
	public void processKeyStroke(char ch) {
		// TODO Auto-generated method stub
		if (ch == 'n') {
			this.changeNspec(Integer.parseInt(JOptionPane.showInputDialog("Enter the desired number of spectra")));
		}
		if (ch == 'w')
			setWidth(Double.parseDouble(JOptionPane.showInputDialog("Enter the width of the cut.")));
		if (ch == 't'){
			this.dh = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired offset in 'arbitrary units'."));
			setOffsetSpec();
			double max = ArrayOps.max(offsetSpec), min = ArrayOps.min(offsetSpec);
			g.setRange(min, max);
			g.repaint();
		}
		if (ch == 'i'){
			drawThin = !drawThin;
			for (int i = 0; i < nspec; i++)
			{
				g.setThin(drawThin, i);
				g.repaint();
			}
		}
		if (ch == 'v') //save everything
		{
			double[] temp;
			double[][] results = new double [nspec+1][t.nlayers+2];
			for (int i = 2; i < t.nlayers+2; i++){
				results[0][i] = t.v[i-2];
				for (int j = 0; j < nspec; j++)
					results[j+1][i] = cut[i-2].result[j];
			}
			for (int j = 0; j < nspec; j++)
			{
				temp = t.getMetricCoords(points[j]);
				results[j+1][0] = temp[0];
				results[j+1][1] = temp[1];
			}
			FileOps.writeTableASCII(results);
		}
		
	}
	public void setData(double[][][] data) {
		for (int i = 0; i < t.nlayers; i++)
			cut[i].changeMap(data[i]);
	}
	public void putLine(double[][] line) {

		line[0][0] = cut[0].r0*Math.cos(cut[0].phi) + t.nx/2;
		line[0][1] = cut[0].r0*Math.sin(cut[0].phi) + t.ny/2;		
		line[1][0] = cut[0].r1*Math.cos(cut[0].phi) + t.nx/2;
		line[1][1] = cut[0].r1*Math.sin(cut[0].phi) + t.ny/2;		
	}
}

package drawing;


import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;

import util.ArrayOps;
import util.FieldOps;
import util.Printer;
import util.SpectraUtil;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.PointSpectra;
import util.fileops.Topomap;

public class LineCutDrawer extends JFrame implements KeyStrokeProccessor {
	
	
	GraphDrawerCart g;
	Topomap t = null;
	PointSpectra tAlt = null;
	int tnlayers = 0;
	double[][][] dataset;
	double[] mean;
	int nspec;
	public JFileChooser fc = new JFileChooser(Topomap.stddir);;
	
	double dh; //height offset.
	double oldDH; //in case of drawing an image, dh will go to zero. 
	double[] x; //energy numbers.
	double[][] points; //[n][2]
	double[][] spec; //[nspec][t.nlayers];
	double[][] offsetSpec;
	boolean drawThin = true;
	
	boolean subtractFirstSpectrum = false;
	boolean usingMap = true;
	boolean drawingLines = true; //else, render the line cut as an image
	int fx, fy; //image size ratios
	
	int scaleIndex = 0;
	
	private double[] tempx, tempy;
	
	Color def = Color.RED;
	
	JMenuBar menuBar;
	
	public LineCutDrawer(Topomap t, int nspec, double[][] line)
	{
		this.t = t;
		tnlayers = t.nlayers;
		x = t.v;
		setDataset(t.data);
		this.nspec = nspec;
		points = new double [nspec][2];
		tempx = new double [nspec];
		tempy = new double [nspec];
		putPoints(line);
//		nspec = points.length;
		spec = new double [nspec][t.nlayers];
		offsetSpec = new double [nspec][t.nlayers];
		
		evaluateSpec();
		setDefaultDH();
		setOffsetSpec();
		g = new GraphDrawerCart("Line Cut Viewer", true);
		setUpDrawer();
		g.showWindow();
		
		//Setting up menu bar
		menuBar = new JMenuBar();
		
	    // File Menu, F - Mnemonic
	    JMenu fileMenu = new JMenu("File");
	    fileMenu.setMnemonic(KeyEvent.VK_F);
	    menuBar.add(fileMenu);
		
	    //File->Save data as csv
	    JMenuItem saveAsCsvMI = new JMenuItem("Save as comma-separated text file");
	    fileMenu.add(saveAsCsvMI);
	    class saveAsCsvAL implements ActionListener{
	    	double[][] data;
	    	saveAsCsvAL(double[][] dataPrime){
	    		data=dataPrime;
	    	}
	    	public void actionPerformed(ActionEvent e){
	    		ColumnIO.writeTableCSV(data, FileOps.selectSave(fc).toString());
	    	}
	    }
	    saveAsCsvMI.addActionListener(new saveAsCsvAL(spec));
	    
	    setJMenuBar(menuBar);
	}
	public LineCutDrawer(PointSpectra t)
	{
		this.t = null;
		tnlayers = t.nlayers;
		x = t.v;
		this.tAlt = t;
		this.dataset = null;
		this.mean = null;
		this.nspec = t.nspec;
		usingMap = false;
		spec = new double [nspec][tnlayers];
		offsetSpec = new double [nspec][tnlayers];
		evaluateSpec();
		setDefaultDH();
		setOffsetSpec();
		g = new GraphDrawerCart("Line Cut Viewer", true);
		setUpDrawer();
		g.extra = this;
		g.showWindow();
	}
	public void setDataset(double[][][] data) {
		this.dataset = data;
		this.mean = new double[data.length];
		for (int i = 0; i < data.length; i++)
			mean[i] = FieldOps.mean(dataset[i]);
	}
	public LineCutDrawer(Layer l, int nspec, double[][] line)
	{
		this.t = Topomap.newTopomap(new Layer[] {l});
		x = t.v;
		this.nspec = nspec;
		points = new double [nspec][2];
		tempx = new double [nspec];
		tempy = new double [nspec];
		putPoints(line);
//		nspec = points.length;
		spec = new double [nspec][tnlayers];
		offsetSpec = new double [nspec][tnlayers];
		
		evaluateSpec();
		setDefaultDH();
		setOffsetSpec();
		g = new GraphDrawerCart("Line Cut Viewer", true);
		setUpDrawer();
		g.showWindow();
	}
	
	private void setOffsetSpec() {
		// TODO Auto-generated method stub
		if (!subtractFirstSpectrum)
		for (int i = 0; i < nspec; i++)
			for (int j = 0; j < tnlayers; j++)
				offsetSpec[i][j] = spec[i][j] + dh*i;
		else
			for (int i = 0; i < nspec; i++)
				for (int j = 0; j < tnlayers; j++)
					offsetSpec[i][j] = spec[i][j] + dh*i - spec[0][j];
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
		putPoints(line);
		spec = new double [nspec][tnlayers];
		offsetSpec = new double [nspec][tnlayers];
		evaluateSpec();
		if (drawingLines) this.setDefaultDH();
		setOffsetSpec();
		setUpDrawer();
	}
	public void setDefaultDH()
	{
		dh = (ArrayOps.max(spec) - ArrayOps.min(spec))/nspec;
	}
	public void evaluateSpec()
	{
		if (usingMap){
		for (int i = 0; i < nspec; i++)
			for (int j = 0; j < tnlayers; j++)
				spec[i][j] = FieldOps.getValueAt(dataset[j], points[i][0], points[i][1], mean[j]);
		}
		else
		{
			for (int i = 0; i < nspec; i++)
				for (int j = 0; j < tnlayers; j++)
					spec[i][j] = tAlt.data[i][j];
		}
		if (g != null && !drawingLines) g.imageScale = ColorScales.getNew(spec, scaleIndex);
	}
	public void setUpDrawer()
	{
		double max = ArrayOps.max(offsetSpec), min = ArrayOps.min(offsetSpec);
		
		g.setNumPlots(nspec);
		for (int i = 0; i < nspec; i++)
		{
			g.setXY(new GraphDrawerCart.GraphObject(x, offsetSpec[i]), false, false, i);
			g.setColor(def, i);
			g.setThin(drawThin, i);
		}
		int nl = (this.t == null ? tAlt.nlayers : this.t.nlayers);
		if (x[0] < x[nl-1])
			g.setDomain(x[0], x[nl-1]);
		else
			g.setDomain(x[nl-1], x[0]);
//		g.setDomain(x[0], x[tnlayers-1]);
		g.setRange(min, max);
		g.repaint();
	}
	
	public void refresh()
	{
		evaluateSpec();
		setOffsetSpec();
		g.setTitle("Line from " + Printer.vectorP(points[0]) + " to " + Printer.vectorP(points[points.length-1]));
		g.repaint();
	}
	public void dispose()
	{
		g.dispose();
	}
	public void putPoints(double[][] line)
	{
		ArrayOps.putArrayInclBoth(line[0][0], line[1][0], nspec, tempx);
		ArrayOps.putArrayInclBoth(line[0][1], line[1][1], nspec, tempy);
		for (int i = 0; i < nspec; i++)
		{
			points[i][0] = tempx[i];
			points[i][1] = tempy[i];
		}
	}

	public void processKeyStroke(char ch) {
		// TODO Auto-generated method stub
		if (ch == 'n') {
			this.changeNspec(Integer.parseInt(JOptionPane.showInputDialog("Enter the desired number of spectra")));
		}
		if (ch == 't'){
			this.dh = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired offset in 'arbitrary units'. \r\n Current value = " + dh + "."));
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
			}
			g.repaint();
		}
		if (ch == 'I')
		{
			switchLineDrawing();
		}
		if (ch == 'v') //save everything
		{
//			if (t == null) FileOps.selectSave(fc);
			double[] temp;
			double[][] results = new double [nspec+1][tnlayers+2];
			double dh = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired offset in 'arbitrary units'. \r\n Current value = " + this.dh + "."));
			double[] x = new double [nspec], y = new double [nspec];
			for (int i = 2; i < tnlayers+2; i++){
				results[0][i] = t == null ? tAlt.v[i-2] : t.v[i-2];
				for (int j = 0; j < nspec; j++)
					results[j+1][i] = spec[j][i-2] + j*dh - (subtractFirstSpectrum ? spec[0][i-2] : 0);
			}
			for (int j = 0; j < nspec; j++)
			{
				temp = t != null ? t.getMetricCoords(points[j]) : new double [] {0, tAlt.y[j]};
				results[j+1][0] = temp[0];
				results[j+1][1] = temp[1];
				x[j] = temp[0];
				y[j] = temp[1];
			}
			FileOps.writeTableASCII(fc, results);
			PointSpectra ps = new PointSpectra(spec, t == null ? tAlt.v : t.v, x, y);
			PointSpectra.writeBIN(ps, fc.getSelectedFile().toString() + ".bin");
			for (int i = 2; i < tnlayers+2; i++){
				results[0][i] = t == null ? tAlt.v[i-2] : t.v[i-2];
				for (int j = 0; j < nspec; j++)
					results[j+1][i] = spec[j][i-2];
			}
			ColumnIO.writeTable(results, fc.getSelectedFile().toString() + "_unshifted.txt");
			
		}
		if (ch == 'V') //save everything
		{
			double[] x = new double [nspec], y = new double [nspec];
			double[] r = new double [2];
			for (int i = 0; i < nspec; i++)
			{
				r = t.getMetricCoords(points[i]);
				x[i] = r[0];
				y[i] = r[1];
			}
			PointSpectra ps = new PointSpectra(spec, t.v, x, y);
			PointSpectra.writeBIN(ps, FileOps.selectSave(fc).toString());
		}
		if (ch == 'R')
		{
			subtractFirstSpectrum = !subtractFirstSpectrum;
			refresh();			
		}
		
	}
	public void setFC(JFileChooser fc)
	{
		this.fc = fc;
	}
	private void switchLineDrawing() {
		// TODO Auto-generated method stub
		drawingLines = !drawingLines;
		g.paintingArrays = drawingLines;
		if (!drawingLines) {oldDH = dh; dh = 0; g.imageScale = ColorScales.getNew(spec, scaleIndex);}
		else setDefaultDH();//dh = oldDH;
		refresh();
	}
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		PointSpectra t = PointSpectra.open(fc);
		LineCutDrawer lcd = new LineCutDrawer(SpectraUtil.averageEachPoint(t.splitByPosition()));
		lcd.fc = fc;
	}


}

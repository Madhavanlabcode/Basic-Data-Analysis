package util.fourier;

import image.ImageEditing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import drawing.GraphDrawerCart;

import main.SRAW;

import util.ArrayOps;
import util.FieldOps;
import util.color.ColorScale;
import util.color.ColorScale2d;
import util.color.ColorScales;
import util.fileops.FileOps;
import util.fileops.Topomap;


public class ImpurityEmergenceEditor extends JFrame implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	
	
	Topomap t;
	int N;
	double fmax, fmin, fdelta;
	double[][] drawField;
	double[][][] drawFieldC;
	boolean real = true;
	
	GraphDrawerCart spec;
	double[] spectrum;
	double[] twoV;
	
	BufferedImage image;
	boolean refresh = true;
	
	ColorScale scale;
	ColorScale2d cscale;
	double scalemin, scalemax;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	
	int[] writepoint = {ox + sx + 50, oy + 10};
	int linesize = 15;
	
	int zoomLevel = 0, zoomFactor = 1;
	int sizeratio = 1;
	
	int WIDTH = 1600, HEIGHT = 1080;
	
	//para is the current index
	int para = 0;
	
	int currentx, currenty;
	int currenti, currentj;
	int calcx = 0, calcy = 0;
	String dir;
	double[] pbounds;
	SliderPanel s;
	int snpts;
	
	ArrayList<Impurity> imps;
	int selectedImp = 0;
	
	String imageOutDir;
	
	//these are for image output purposes only.
	BufferedImage imagePink;
	double[][] pinkRing;
	double ringR = 10, ringT = 0.5;
	
	public ImpurityEmergenceEditor(Topomap t, String dir)
	{
		snpts = t.nlayers-1;
		this.t = t;
		this.dir = dir;
		
		imps = getImpurities(FileOps.selectOpen(null), 1);
		
		spectrum = new double [t.nlayers];
//		t.putSpectrum(0, 0, spectrum);
		putImpuritySpectrum(4);
		
		spec = new GraphDrawerCart("Spectrum", t.v, spectrum);
		twoV = new double[] {t.v[0], t.v[0]};
		spec.setXY(new GraphDrawerCart.GraphObject(twoV, new double[] {ArrayOps.min(t.data), ArrayOps.max(t.data)}), false, false, 1);
		spec.showWindow();
		
		JOptionPane.showMessageDialog(null, "Select the image output folder.");
		imageOutDir = FileOps.selectDir(null);
		
		N = t.nx;
		sx = N;
		drawField = t.data[0];
		while (drawField.length > 1024)
		{
			this.drawField = FieldOps.reduce(2, this.drawField);
			sx = this.drawField.length;
			zoomLevel++;
			zoomFactor *= 2;
		}
		while (sizeratio*sx < 1024)
		{
			sizeratio*=2;
		}
		
		if (sx*sizeratio == 1024) writepoint[0] += 512;
		
		sy = this.drawField[0].length;
		setFieldInfo();
		image = new BufferedImage(sx*sizeratio, sy*sizeratio, BufferedImage.TYPE_INT_RGB);
		
		imagePink = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		pinkRing = new double [t.nx][t.ny];
		
		resetColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data), true);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		showWindow();
		setTitle("Title");

		s = new SliderPanel(this, new JFrame());
		s.show();
		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		addMouseListener(this);
	}
	public void setFieldInfo()
	{
//		double[] bounds = calc.getMinMax(real);
		double[] bounds = null;
		bounds = new double[] {ArrayOps.min(t.data[para]), ArrayOps.max(t.data[para])};
			
		fmax = bounds[1]; fmin = bounds[0];
		fdelta = fmax - fmin;
		setTitle("Range: [" + fmin + ", " + fmax + "]" + "     " + "p = " + para);
	}
	public void resetColorScale(double downnum, double upnum, boolean hardLimits)
	{
		scalemax = fmin + fdelta*upnum;
		scalemin = fmin + fdelta*downnum;
		
		if (hardLimits) scale = new ColorScales.LinearBRYW(upnum, downnum);
		else scale = new ColorScales.LinearBRYW(scalemax, scalemin);
	}
	public void resetColorScale()
	{
		if (real) scale = new ColorScales.LinearBRYW(scalemax, scalemin);
		else cscale = new ColorScales.MYC2d(scalemax, scalemin, 2*Math.PI);
	}

	public void formImage()
	{
		if (real) SRAW.writeImage(image, drawField, scale, sizeratio);
		else SRAW.writeImage(image, drawFieldC, cscale);
	}
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		repaint();
		show();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}	
	public void paint(Graphics g)
	{
		if (refresh){
			g.clearRect(0, 0, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			refresh = false;
		}
		drawRings(g);
		drawText(g);
	}
	public void drawText(Graphics g)
	{
		g.clearRect(writepoint[0], writepoint[1]-linesize, sx, sy);
		String it = "";
		String[] lines = it.split("\r\n");
		for (int i = 0; i < lines.length; i++)
			g.drawString(lines[i], writepoint[0], writepoint[1] + i*linesize);
	}
	public void drawRings(Graphics g)
	{
		double rx, ry;
		int r = (int)(ringR*sizeratio);
		g.setColor(Color.MAGENTA);
		for (int i = 0; i < imps.size(); i++)
		{
			if (i == selectedImp)
			{
				rx = windowX(imps.get(i).position[0]);
				ry = windowY(imps.get(i).position[1]);
				g.drawOval((int)(rx+0.5)-r, (int)(ry+0.5)-r, 2*r, 2*r);
			}
		}
	}
	public double windowX(double datax)
	{
		return sizeratio*datax + ox;
	}
	public double windowY(double datay)
	{
		return sizeratio*datay + oy;
	}


	public void resetGraphics()
	{
		drawField = t.data[para];
		while (drawField.length > 1024)
		{
			this.drawField = FieldOps.reduce(2, this.drawField);
			sx = this.drawField.length;
			zoomLevel++;
			zoomFactor *= 2;
		}
		setFieldInfo();
//	    resetColorScale();
	    this.formImage();
	    refresh = true;
	    repaint();
	}
	//if 
	public void keyPressed(KeyEvent arg0) {
	}
	public void keyReleased(KeyEvent arg0) {
	}
	public void keyTyped(KeyEvent arg0) {
		System.out.println(arg0.getKeyChar());
		if (arg0.getKeyChar() == ' '){
			resetColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data), true);
		}
		if (arg0.getKeyChar() == 'n')
		{
			imps.get(selectedImp).setEmergence(para, t.v[para]);
			//prepare pink ring image:
			SRAW.writeImage(imagePink, t.data[t.nlayers-1], scale);
			addPinkRing(selectedImp);
			ImageEditing.moveEachPixelTowardAColor(imagePink, pinkRing, Color.MAGENTA, imagePink);
			removePinkRing(selectedImp);
			SRAW.writeImage(imageOutDir + "imp" + selectedImp, imagePink);
			//print graph:
			SRAW.writeImage(imageOutDir + "spec" + selectedImp, (BufferedImage)spec.dbimage);

			selectedImp++;
			if (selectedImp >= imps.size())
			{
				System.out.println(getImpList(""));
				System.exit(0);
			}
			putImpuritySpectrum(4);
			spec.repaint();
			refresh = true;
			repaint();
			
		}
	}
	public void mouseWheelMoved(MouseWheelEvent arg0) {
	}
	public void mouseDragged(MouseEvent arg0) {
	}
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		calcx = currentx - ox;
		calcy = currenty - oy;
		currenti = calcx/sizeratio;
		currentj = calcy/sizeratio;
		
//		if (withinBounds(currenti, currentj, t.data[0]))
//			t.putSpectrum(currenti, currentj, spectrum);
//		spec.repaint();
		
//		refresh = false;
//		repaint();
	}
	public void mouseClicked(MouseEvent arg0) {
	}
	public void mouseEntered(MouseEvent arg0) {
	}
	public void mouseExited(MouseEvent arg0) {
	}
	public void mousePressed(MouseEvent arg0) {
	}
	public void mouseReleased(MouseEvent arg0) {
	}
	
	public void putImpuritySpectrum(double radius)
	{
		for (int j = 0; j < t.nlayers; j++)
			spectrum[j] = FieldOps.getLocalAvgCircle(t.data[j], imps.get(selectedImp).position[0], imps.get(selectedImp).position[1], radius);
	}
	public String getImpList(String name)
	{
		String implist = "Impurity " + name + " position list.\r\n";
		implist += "()\r\n";
		implist += "X \t Y \t Index of emergence\t Energy of Emergence \r\n";
		for (int i = 0; i < imps.size(); i++)
		{
			implist += imps.get(i).toString() + "\r\n";
		}
		return implist;
	}
	private void makeSmallPinkRing()
	{
		double rangemax = ringR + 4*ringT;
		double g, r;
		for (int i = 0; i < sx; i++){
			for (int j = 0; j < sy; j++){ pinkRing[i][j] = 0;
				for (int k = 0; k < imps.size(); k++)
				{
					if (Math.abs(i - imps.get(k).possmall[0]) < rangemax && Math.abs(j - imps.get(k).possmall[1]) < rangemax && k == selectedImp)
					{
						r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
						g = Math.exp(-(r-ringR)*(r-ringR)/(ringT*ringT));
						pinkRing[i][j] += g;
					}
				}
				pinkRing[i][j] = Math.min(pinkRing[i][j], 1);
			}
		}
	}
	private void removePinkRing(int k)
	{
		double rangemax = ringR + 4*ringT;
		double g, r;
		int imin, imax, jmin, jmax;
		imin = (int)Math.max(imps.get(k).possmall[0]-rangemax, 0);
		jmin = (int)Math.max(imps.get(k).possmall[1]-rangemax, 0);
		imax = (int)Math.min(imps.get(k).possmall[0]+rangemax, sx-1);
		jmax = (int)Math.min(imps.get(k).possmall[1]+rangemax, sy-1);
		for (int i = imin; i < imax; i++){
			for (int j = jmin; j < jmax; j++){
				r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
				g = Math.exp(-(r-ringR)*(r-ringR)/(ringT*ringT));
				pinkRing[i][j] -= g;
				pinkRing[i][j] = Math.max(pinkRing[i][j], 0);
			}
		}
	}
	private void addPinkRing(int k)
	{
		double rangemax = ringR + 4*ringT;
		double g, r;
		int imin, imax, jmin, jmax;
		imin = (int)Math.max(imps.get(k).possmall[0]-rangemax, 0);
		jmin = (int)Math.max(imps.get(k).possmall[1]-rangemax, 0);
		imax = (int)Math.min(imps.get(k).possmall[0]+rangemax, sx-1);
		jmax = (int)Math.min(imps.get(k).possmall[1]+rangemax, sy-1);
		for (int i = imin; i < imax; i++){
			for (int j = jmin; j < jmax; j++){
				r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
				g = Math.exp(-(r-ringR)*(r-ringR)/(ringT*ringT));
				pinkRing[i][j] += g;
				pinkRing[i][j] = Math.min(pinkRing[i][j], 1);
			}
		}
	}

	public static class SliderPanel extends JPanel implements ChangeListener
	{
		ImpurityEmergenceEditor parent;
		JFrame frame;
		public JSlider s;
		public JSlider min, max;
		int parentmember;
		int oldvalue = 0;
		int oldminv = 0, oldmaxv = 999;
		static int npts = 1001;
		
		public SliderPanel(ImpurityEmergenceEditor parent, JFrame frame)
		{
			npts = parent.snpts;
			this.frame = frame;
			this.parent = parent;
			this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			this.setBorder(new LineBorder(Color.GRAY));
			s = new JSlider(0, npts, 0);
			s.setValue(0);
			
			s.setSnapToTicks(true);
			s.setMinorTickSpacing(1);
			s.addChangeListener(this);
			
			min = new JSlider(0, 1000, oldminv);
			max = new JSlider(0, 1000, oldmaxv);
			min.addChangeListener(this);
			max.addChangeListener(this);
			add(s);
			add(min);
			add(max);
			frame.setSize(5*npts+40, 80);
			frame.add(this);
		}
		public void show(){
			frame.show();
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
//			if (s.getValue() == oldvalue && min.getValue() == oldminv && max.getValue() == oldmaxv) return;
			if (min.getValue() != oldminv || max.getValue() != oldmaxv)
			{
				parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000), false);
				parent.formImage();
				parent.refresh = true;
				parent.repaint();
				oldminv = min.getValue();
				oldmaxv = max.getValue();
				return;
			}
//			else if (s.getValue() == oldvalue) return;
//			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
			parent.para = s.getValue();
			parent.twoV[0] = parent.t.v[parent.para];
			parent.twoV[1] = parent.t.v[parent.para];
			parent.spec.repaint();
			
			parent.resetGraphics();
			frame.setTitle(parent.t.names[parent.para]);
			System.out.println(s.getValue());
			oldvalue = s.getValue();
		}
	}
	 public static class Impurity
	 {
		 public double[] position;
		 double[] possmall;
		 boolean isinbox = false;
		 
		 int emergenceIndex;
		 double emergenceEnergy;
		 double[] intensity;
		 
		 public Impurity(double x, double y, int zoomfactor, int emergenceIndex, double intensity, int ilength)
		 {
			 position = new double[] {x, y};
			 possmall = new double[] {x/zoomfactor, y/zoomfactor};
			 this.emergenceIndex = emergenceIndex;
			 this.intensity = new double[ilength];
			 this.intensity[0] = intensity;
		 }
		 public void setIntensity(double i, int index)
		 {
			 intensity[index] = i;
		 }
		 public void setEmergence(int i, double d)
		 {
			 emergenceIndex = i;
			 emergenceEnergy = d;
		 }
		 public String toString()
		 {
			 String x = "";
			 
			 x += position[0] + "\t" + position[1] + "\t" + emergenceIndex + "\t" + emergenceEnergy;
			 for (int i = 0; i < intensity.length; i++)
				 x += "\t" + intensity[i];
			 return x;
		 }
	 }
	 
	public static ArrayList<Impurity> getImpurities(File f, int zoomFactor) {
		ArrayList<Impurity> imps = new ArrayList<Impurity>();
		Scanner in = null;
		String line;
		String[] words;
		try {
			in = new Scanner(f);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (int i = 0; i < 3; i++)
			in.nextLine();
		while(in.hasNextLine())
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			imps.add(new Impurity(Double.parseDouble(words[0]), Double.parseDouble(words[1]), 1, 1, 1, 1));
		}
		return imps;
	}
	public static ArrayList<Impurity> getImpuritiesEmergence(File f) {
		ArrayList<Impurity> imps = new ArrayList<Impurity>();
		Scanner in = null;
		String line;
		String[] words;
		try {
			in = new Scanner(f);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		for (int i = 0; i < 3; i++)
			in.nextLine();
		while(in.hasNextLine())
		{
			line = in.nextLine();
			line = line.trim();
			words = line.split("\t");
			imps.add(new Impurity(Double.parseDouble(words[0]), Double.parseDouble(words[1]), 1, 1, 1, 1));
			imps.get(imps.size()-1).setEmergence(Integer.parseInt(words[2]), Double.parseDouble(words[3]));
		}
		return imps;
	}
	
	public static void makeEmergenceHist()
	{
		ArrayList<Impurity> emer = getImpuritiesEmergence(FileOps.selectOpen(null));
		double max = 0.15, min = 0;
		int npts = 100;
		double dx = (max-min)/(npts-1);
		double e;
		int[] nimps = new int [npts];
		for (int i = 0; i < npts; i++)
		{
			e = min + dx*i;
			for (int j = 0; j < emer.size(); j++)
				if (emer.get(j).emergenceEnergy < e)
					nimps[i]++;
		}
		System.out.println("Energy \t Number of visible impurities");
		for (int i = 0; i < npts; i++)
		{
			System.out.println((min+dx*i) + "\t" + nimps[i]);
		}
	}
	public static boolean withinBounds(int i, int j, double[][] array)
	{
		return i > 0 && j > 0 && i < array.length && j < array[i].length;
	}
	public static void main(String[] args)
	{
//		new ImpurityEmergenceEditor(Topomap.open(), Topomap.stddir);
		makeEmergenceHist();
	}

}

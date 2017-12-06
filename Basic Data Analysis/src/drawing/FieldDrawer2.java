package drawing;

import java.awt.Graphics;
import java.awt.event.*;
import java.awt.image.BufferedImage;

import javax.swing.JFileChooser;
import javax.swing.JFrame;

import main.SRAW;
import util.color.ColorScale;
import util.color.ColorScale2d;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.RHKFileOps;
import util.fourier.FFTOps;
import util.*;

//This class must draw a sx*sy section of a field, zoomed in to arbitrary amounts.

public class FieldDrawer2 extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	private static final long serialVersionUID = 7335768275347252739L;
	//The dimensions of field will be data.length/(2**zoomedIn)
	double field[][]; double fmax, fmin, fdelta;
	double phase[][] = null;
	
	double[][] data;
	int sizeFactor;
	//zoomedIn is the number of times the field has been zoomed in
	//from the data. each time is a power of 2.
	int zoomedIn;
	
	int[] mouseDelta = null;
	
	boolean twoComp = false;
	
	BufferedImage image, drawImage;
	
	ColorScale scale;
	ColorScale2d scale2;
	
	GradCalculator gradcalc;
	CurlCalculator curlcalc = null;
	
	int[] stringr = {600, 100};
	int liney = 20;
	int[] disppt = new int [2];
	//we assume sx = sy.
	int sx = 512, sy = 512;
	int Nx, Ny;
	int fx, fy;
	int ox = 20, oy = 40;
	
	//physical dimensions:
	double Lx, Ly;
	double ax, ay; //per pixel
	
	int[] corner = {0, 0}; //The field coordinates of the corner. The pixel coordinates of the
	//corner are ox, oy.
	
	//local statistics;
	int imin = 287, jmin = 192, imax = 307, jmax = 212;
//	int imin = 240, jmin = 240, imax = 260, jmax = 260;
	double[] centroid = new double [2], sigma = new double [2];
	double fmean, fstdev, fminloc, fmaxloc;
	int[] minR = new int [2], maxR = new int [2];
	int[] statpt = {600, 300};
	
	int WIDTH = 800, HEIGHT = 600;
	
	public DataZoomPanel panel = null;
	
	public FieldDrawer2(double[][] field)
	{
		this.field = field;
		sx = field.length;
		sy = field[0].length;
		setFieldInfo();
		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		gradcalc = new GradCalculator(this, 0);
		resetColorScale(0, 1);
		formImage();
		//gradcalc.activate();
		showWindow();
		
		addKeyListener(this);
		addMouseWheelListener(this);
	}
	public FieldDrawer2(double[][] field, double[][] phase)
	{
		this.field = field;
		this.phase = phase;
		twoComp = true;
		sx = field.length;
		sy = field[0].length;
		setFieldInfo();
		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		curlcalc = new CurlCalculator(this);
		resetColorScale(0, 1);
		formImage();
		//gradcalc.activate();
		showWindow();
		
		addKeyListener(this);
		addMouseWheelListener(this);
	}
	
	//Here we assume that the data is a power of 2 times sx x sy
	public FieldDrawer2(String binFile)
	{
		data = ColumnIO.readSquareTable(binFile);
		Nx = data.length; Ny = data[0].length;
		sizeFactor = Nx/sx;
		zoomedIn = log2(sizeFactor);
		drawImage = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		formField();
//		setFieldInfo();
//		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
//		gradcalc = new GradCalculator(this, 0);
//		resetColorScale(0, 1);
//		formImage();
		//gradcalc.activate();
		showWindow();
		
		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseListener(this);
		addMouseMotionListener(this);
	}
	public FieldDrawer2(String binFile, double Lx, double Ly)
	{
		this.Lx = Lx; this.Ly = Ly;
		
		data = ColumnIO.readSquareTable(binFile);
		Nx = data.length; Ny = data[0].length;
		sizeFactor = Nx/sx;
		zoomedIn = log2(sizeFactor);
		drawImage = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		formField();
//		setFieldInfo();
//		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
//		gradcalc = new GradCalculator(this, 0);
//		resetColorScale(0, 1);
//		formImage();
		//gradcalc.activate();
		showWindow();
		
		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseListener(this);
		addMouseMotionListener(this);
	}
	public FieldDrawer2(double[][] data, double Lx, double Ly)
	{
		this.Lx = Lx; this.Ly = Ly;
		
		this.data = data;
		Nx = data.length; Ny = data[0].length;
		sizeFactor = Nx/sx;
		zoomedIn = log2(sizeFactor);
		drawImage = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		formField();
		getLocalStatistics();
//		setFieldInfo();
//		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
//		gradcalc = new GradCalculator(this, 0);
//		resetColorScale(0, 1);
//		formImage();
		//gradcalc.activate();
		showWindow();
		
		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseListener(this);
		addMouseMotionListener(this);
		@SuppressWarnings("unused")
		FieldPanel panel = new FieldPanel(this);
	}
	
	void formField()
	{
		int factor = twoTT(zoomedIn);
		fx = Nx/factor; fy = Ny/factor;
		ax = Lx/fx; ay = Ly/fy;
		field = new double[fx][fy];
		//factor is always >= 1.
		if (factor > 1)
		{
			FieldOps.reduce(factor, data, field);
		}
		else
			field = data;
		
		setFieldInfo();
		resetColorScale(0, 1);
		//Now, we have the field which may be smaller than data by some power of 2.
		formImage();
		formDrawImage();
		setTitle("Lx = " + Lx + "\t   Ly = " + Ly);
	}
	
	public static int log2(int mult2)
	{
		int x = 0;
		int i = mult2;
		while(i > 1)
		{
			x++;
			i /= 2;
		}
		return x;
	}
	public static int twoTT(int x)
	{
		int i = x;
		int a = 1;
		while (i > 0)
		{
			a *= 2;
			i--;
		}
		return a;
	}
	
	//Assumes field size does not change
	public void resetField(double[][] field, boolean changeColors)
	{
		this.field = field;
		setFieldInfo();
		if (changeColors)
			resetColorScale(0, 1);
		formImage();
		repaint();
	}
	
	public void changeField(double[][] field, boolean changeColors)
	{
		this.field = field;
		setFieldInfo();
		if (changeColors)
			resetColorScale(0, 1);
		gradcalc.reset();
		formImage();
		repaint();
		gradcalc.active = false;
	}
	public void setFieldInfo()
	{
		fmax = max(field); fmin = min(field);
		fdelta = fmax - fmin;
		System.out.println("[" + fmin + ", " + fmax + "]");
	}
	public void redrawColorScale(double downnum, double upnum)
	{
		resetColorScale(downnum, upnum);
		formImage();
		formDrawImage();
		repaint();
	}
	public void resetColorScale(double downnum, double upnum)
	{
		if(!twoComp)
			scale = new ColorScales.LinearBRYW(fmin + fdelta*upnum, fmin + fdelta*downnum);
		else
			scale2 = new ColorScales.MYC2d(fmin + fdelta*upnum, fmin + fdelta*downnum, 2*Math.PI);
	}
	
	//This includes resizing the image.
	//It is necessary to hold the subimage
	public void formImage()
	{
//		image.
		image = new BufferedImage(fx, fy, BufferedImage.TYPE_INT_RGB);
		if(!twoComp)
			for (int i = 0; i < fx; i++)
				for (int j = 0; j < fy; j++)
					image.setRGB(i, j, scale.of(field[i][j]).getRGB());
		else
			for (int i = 0; i < fx; i++)
				for (int j = 0; j < fy; j++)
					image.setRGB(i, j, scale2.of(field[i][j], phase[i][j]).getRGB());
			
	}
	public void formDrawImage()
	{
//		image.
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				drawImage.setRGB(i, j, image.getRGB(corner[0]+i, corner[1]+j));
			
	}

	public int[] getFieldPoint(int x, int y)
	{
		return new int [] {x - ox, y - oy};
	}
	
	public int[] getQuad(int x, int y)
	{
		int[] r = getFieldPoint(x, y);
		r[0] *= 2; r[0] /= sx;
		r[1] *= 2; r[1] /= sy;
		return r;
	}
	public void drawImage(Graphics g)
	{
		g.drawImage(drawImage, ox, oy, null);
//		g.drawImage(image, ox, oy, ox+sx, oy+sy, corner[0], corner[1], corner[0]+sx, corner[1]+sy, null);
	}
	
	public void paint(Graphics g)
	{
		g.clearRect(0, 0, WIDTH, HEIGHT);
		drawImage(g);
		drawFieldPt(g);
		drawStatistics(g);
	}
	public void drawFieldPt(Graphics g)
	{
		if (disppt[0] >= 0 && disppt[1] >= 0 && disppt[0] < fx && disppt[1] < fy)
		{
			g.drawString("" +field[disppt[0]][disppt[1]], stringr[0], stringr[1]);
			g.drawString("at [" + disppt[0] + ", " + disppt[1] + "]", stringr[0], stringr[1] + liney);
			g.drawString("at [" + (disppt[0] - fx/2) + ", " + (disppt[1] - fy/2) + "] from center", stringr[0], stringr[1] + 2*liney);
			g.drawString("x = " + (disppt[0] - fx/2)*ax, stringr[0], stringr[1] + 3*liney);
			g.drawString("y = " + (disppt[1] - fy/2)*ay, stringr[0], stringr[1] + 4*liney);
			g.drawString("r = " + Complex.mag((disppt[1] - fy/2)*ay, (disppt[0] - fx/2)*ax), stringr[0], stringr[1] + 5*liney);
			g.drawString("2pi/r = " + 2*Math.PI/Complex.mag((disppt[1] - fy/2)*ay, (disppt[0] - fx/2)*ax), stringr[0], stringr[1] + 6*liney);
		}
	}
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		repaint();
		setVisible(true);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	@Override
	public void keyPressed(KeyEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyReleased(KeyEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyTyped(KeyEvent arg0) {
		// TODO Auto-generated method stub
		if (arg0.getKeyChar() == 'g')
		{
			System.out.println('g');
			if (gradcalc.active) gradcalc.deactivate();
			else if (!gradcalc.active) gradcalc.activate();
		}
		if (arg0.getKeyChar() == 'c')
		{
			System.out.println('c');
			if (curlcalc.active) curlcalc.deactivate();
			else if (!curlcalc.active) curlcalc.activate();
		}
		if (arg0.getKeyChar() == 's')
		{
			System.out.println('s');
			JFileChooser s = new JFileChooser();
			if(s.showSaveDialog(this) == JFileChooser.APPROVE_OPTION);
				SRAW.writeImage(s.getSelectedFile().toString(), drawImage);
		}
//		if (arg0.getKeyChar() == 's' || arg0.getKeyChar() == 'S' && arg0.getModifiers() == KeyEvent.CTRL_DOWN_MASK)
//		{
//			JFileChooser chooser = new JFileChooser();
//			chooser.showSaveDialog(this);
//		}
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent arg0) {
		// TODO Auto-generated method stub
//		int x = arg0.getX();
//		int y = arg0.getY();
		int value = arg0.getWheelRotation();
		System.out.println(value);
		if (zoomedIn + value >= 0);
		{
			zoomedIn -= value;
			formField();
		}
		repaint();
		
		
			
	}

	@Override
	public void mouseDragged(MouseEvent arg0) {
		// TODO Auto-generated method stub
		System.out.print("[" + arg0.getX() + ", " + arg0.getY() + "]\t");
		if (mouseDelta == null){
			mouseDelta = new int [2];
			mouseDelta[0] = arg0.getX(); mouseDelta[1] = arg0.getY();
		}
		else if (corner[0] >= 0 && corner[1] >= 0 && corner[0] + sx < fx && corner[1] + sy < fy)
		{	
			mouseDelta[0] -= arg0.getX(); mouseDelta[1] -= arg0.getY();
			System.out.println("[" + mouseDelta[0] + ", " + mouseDelta[1] + "]\t");
			corner[0] += mouseDelta[0];
			corner[1] += mouseDelta[1];
			corner[0] = Math.max(corner[0], 0); corner[0] = Math.min(corner[0], fy - sy -1);
			corner[1] = Math.max(corner[1], 0); corner[1] = Math.min(corner[1], fy - sy -1);
			mouseDelta[0] = arg0.getX(); mouseDelta[1] = arg0.getY();
		}
		formDrawImage();
		repaint();
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		int x = arg0.getX();
		int y = arg0.getY();
		disppt[0] = corner[0] + x - ox; disppt[1] = corner[1] + y - oy;
		repaint();
		
	}

	@Override
	public void mouseClicked(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseEntered(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseExited(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mousePressed(MouseEvent arg0) {
		// TODO Auto-generated method stub
		if (mouseDelta == null)
			mouseDelta = new int [2];
		mouseDelta[0] = arg0.getX(); mouseDelta[1] = arg0.getY();
		
	}

	@Override
	public void mouseReleased(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	public static double max(double[][] x)
	{
		return ArrayOps.max(x);
	}
	public static double min(double[][] x)
	{
		return ArrayOps.min(x);
	}
	public static double[][][] grad(double[][] field)
	{
		return FieldOps.gradient(field);
	}
	public static double[][] curl(double[][][] field)
	{
		return FieldOps.curl(field);
	}
	public static double[][] copy(double[][] x)
	{
		double[][] copy = new double [x.length][x[0].length];
		for (int i = 0; i < x.length; i++)
			for (int j = 0; j < x[0].length; j++)
			{
				copy[i][j] = x[i][j];
			}
		return copy;
	}
	
	//replaces the field by its directional derivative.
	public static class GradCalculator
	{
		boolean active = false;
		public FieldDrawer2 parent;
		public double [][][] grad;
		public double [][] field;
		public double [][] dfds;
		public double theta;
		
		public GradCalculator(FieldDrawer2 parent, double theta)
		{
			this.parent = parent;
			this.field = copy(parent.field);
			grad = grad(field);
			dfds = new double[grad[0].length][grad[0][0].length];
			this.theta = theta;
			setDeriv();
		}
		
		//assumes the field is the same size.
		public void reset(){
			this.field = copy(parent.field);
			grad = grad(field);
			setDeriv();
		}
		
		public void setDeriv()
		{
			for (int i = 0; i < grad[0].length; i++)
				for (int j = 0; j < grad[0][0].length; j++)
				{
					dfds[i][j] = grad[0][i][j]*Math.cos(theta) + grad[1][i][j]*Math.sin(theta);
				}
		}
		public void resetTheta(double theta)
		{
			this.theta = theta;
			setDeriv();
			parent.resetField(dfds, false);
		}
		public void activate()
		{
			parent.resetField(dfds, true);
			active = true;
		}
		public void deactivate()
		{
			parent.resetField(field, true);
			active = false;
		}
	}
	public static class CurlCalculator
	{
		boolean active = false;
		public FieldDrawer2 parent;
		public double [][][] cart;
		public double [][] curl;
		public double theta;
		
		public CurlCalculator(FieldDrawer2 parent)
		{
			this.parent = parent;
			int N = parent.field.length, M = parent.field[0].length;
			cart = new double[N][M][2];
			setCart();
			curl = curl(cart);
			FieldOps.abs(curl);
		}
		
		void setCart()
		{
			for (int i = 0; i < parent.field.length; i++)
				for (int j = 0; j < parent.field[0].length; j++)
				{
					cart[i][j][0] = parent.field[i][j]*Math.cos(parent.phase[i][j]);
					cart[i][j][1] = parent.field[i][j]*Math.sin(parent.phase[i][j]);
				}
			
		}
		//assumes the field is the same size.
		public void reset(){
			setCart();
			curl = curl(cart);
		}
		
		public void activate()
		{
			parent.resetField(curl, true);
			parent.twoComp = false;
			active = true;
		}
		public void deactivate()
		{
			parent.resetField(null, true);
			active = false;
		}
	}
	 public double[][] getField() {
			return field;
		}
	 public BufferedImage getImage() {
			return image;
		}
	 
	 
	 //Run methods.
	 public static void main(String[] args)
	 {
//		 String dir = "C:\\data\\lintrans\\run146topo6_620\\cutoff_pt4\\";
////		 String file = "C:\\data\\lawlerhex\\8302010_0010\\topo.dat";
////		 double[] info = RHKFileOps.readTopoInfo("C:\\data\\lawlerhex\\8302010_0010\\topoinfo.txt");
////		 FieldDrawer2 drawer = new FieldDrawer2(file, info[0], info[1]);
//		 loadFFT(dir, "La-did", true);
	 }
	 public static void loadFFT(String dir, String name, boolean log)
	 {
		 double[] info = RHKFileOps.readTopoInfo(dir + name + ".txt");
////		 double[] info = RHKFileOps.topoInfo(dir + name + ".txt").toDarray();
		 double[][] fftmag = FFTOps.obtainFFTmagCent(ColumnIO.readSquareTable(dir + name + ".dat"));
		 if (log) FieldOps.log(fftmag);
		 int N = fftmag.length;
		 @SuppressWarnings("unused")
		FieldDrawer2 drawer = new FieldDrawer2(fftmag, 2*Math.PI*N/info[0], 2*Math.PI*N/info[1]);
	 }
	 
	 public void drawStatistics(Graphics g)
	 {
		 	g.drawString("For x = [" + imin + ", " + imax + "),   y = [" + jmin + ", " + jmax + ")", statpt[0], statpt[1]);
			g.drawString("<x> = " + centroid[0], statpt[0], statpt[1] + liney);
			g.drawString("<y> = " + centroid[1], statpt[0], statpt[1] + 2*liney);
			g.drawString("Sigma_X = " + sigma[0], statpt[0], statpt[1] + 3*liney);
			g.drawString("Sigma_Y = " + sigma[1], statpt[0], statpt[1] + 4*liney);
			g.drawString("Mean = " + fmean, statpt[0], statpt[1] + 5*liney);
			g.drawString("Sigma = " + fstdev, statpt[0], statpt[1] + 6*liney);
			g.drawString("Min = " + fminloc, statpt[0], statpt[1] + 7*liney);
			g.drawString(" at (" + minR[0] + ", " + minR[1] + ")" , statpt[0], statpt[1] + 8*liney);
			g.drawString("Max = " + fmaxloc, statpt[0], statpt[1] + 9*liney);
			g.drawString(" at (" + maxR[0] + ", " + maxR[1] + ")" , statpt[0], statpt[1] + 10*liney);
	 }
	 public void getLocalStatistics()
	 {
		 centroid[0] = FieldOps.centX(field, imin, imax, jmin, jmax);
		 centroid[1] = FieldOps.centY(field, imin, imax, jmin, jmax);
		 sigma[0] = FieldOps.sigmaX(field, imin, imax, jmin, jmax);
		 sigma[1] = FieldOps.sigmaY(field, imin, imax, jmin, jmax);
		 fmean = FieldOps.mean(field, imin, imax, jmin, jmax);
		 fstdev = FieldOps.sigma(field, imin, imax, jmin, jmax, fmean);
		 fminloc = FieldOps.min(field, imin, imax, jmin, jmax);
		 fmaxloc = FieldOps.max(field, imin, imax, jmin, jmax);
		 minR = FieldOps.minR(field, imin, imax, jmin, jmax);
		 maxR = FieldOps.maxR(field, imin, imax, jmin, jmax);
	 }
}

package drawing;

import image.ImageEditing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import main.SRAW;


import schrodinger.MovieMaker;
import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.FileOps;
import util.fileops.Topomap;
import util.geom.AltMask;
import util.geom.Mask;
import util.*;
//This class only draws real scalar fields.


public class GaussSquareImpurityAngleSuite extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	private static final long serialVersionUID = -7474746993486687545L;
	String dir;
	Topomap t;
	double[][] masterField; int N;
	double field[][]; double fmax, fmin, fdelta;
	boolean twoComp = false;
	BufferedImage image, imageft;
	
	int currentImp = 0;
	double[] imageCorner = new double [2];
	ArrayList<Impurity> imps;
	JFileChooser fc;
	
	int npts = 1000;
	double[][] refMask;
	double[][] tempMask;
	double theta = 0, tempTheta = 0;
	double blendp = 0;
	int width = 40; double[] squareCent = new double[] {64,64};
	double[][] matrix = new double [2][2];
	int[] rotO;
	double corr;
	double[] angles;
	double[][] correlations; //[i] is imp index, j is theta value.
	double[] maxCorrAngles;
	double sqThickness = 12;
	int maskType = 0;
	int nMasks = 4;
	double thmax = Math.PI/2;
	int[] bestLayers = new int []
			{
				395, 395, 395, 304, 258, 314, 395, 325, 271, 257,
				257, 328, 395, 395, 332, 332, 312, 312, 354, 305,
				283, 312, 304, 334, 292, 248, 318, 339, 291, 319,
				319, 376, 389, 368, 310, 269, 395, 395, 367, 391,
				340, 382, 356, 313, 340, 270, 391, 395, 374, 368,
				276, 361, 276, 380, 293, 332
			};
	
	ColorScale scale, ftscale;
	
	GraphDrawerCart spec;
	double[] spectrum;
	double[] twoV;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	int ftox = ox + sx + 10, ftoy = oy;
	
	int ftsizemin = 512;
	int ftsizeratio = 1;
	int sizeratio = 1;
	
	int[] writepoint = {ox, oy + sy + 10};
	
	int zoomLevel = 0, zoomFactor = 1;
	
	int WIDTH = 1080, HEIGHT = 570;
	
	double[][] ftmag; 
	int[] ftwritepoint = {ftox, ftoy};
	int boxsize;
	
	double ftmax, ftmin;
	double upfrac = 1, downfrac = 0;
	double ftscalemin, ftscalemax;
	
	int currentx, currenty;
	int[] inField = new int [2]; //the current position with respect to the image pixels, after scaling.

	//We intend to "up-sample" the impurity to at least 128x128.
	public GaussSquareImpurityAngleSuite(Topomap t, String dir)
	{
		this.dir = dir;
		System.out.println(dir);
		this.t = t;
		
		fc = new JFileChooser(dir);
		imps = getImpurities(FileOps.selectOpen(fc));	
		angles = ArrayOps.generateArrayNotInclUpper(0, thmax, npts);
		correlations = new double [imps.size()][npts];
		maxCorrAngles = new double [imps.size()];
		
		masterField = new double [128][128];
		ftmag = new double [128][128];
		refMask = new double [128][128];
		tempMask = new double [128][128];

		setUpForCurrentImp();
		putSquareMask();
		
		N = masterField.length;
		boxsize = N;
		sx = N;
		rotO = new int [] {N/2, N/2};
		this.field = masterField;
		while (sx > 512)
		{
			this.field = FieldOps.reduce(2, this.field);
			sx = this.field.length;
			zoomLevel++;
			zoomFactor *= 2;
		}
		
		sy = this.field[0].length;
		ftox = ox + ftsizemin + 10; ftoy = oy;
		writepoint = new int[] {ox, oy + sy + 10};
		setFieldInfo();
		
		sizeratio = ftsizemin/sx;
		ftsizeratio = ftsizemin/sx;
		if (sx < 512)
			image = new BufferedImage(ftsizemin, ftsizemin, BufferedImage.TYPE_INT_RGB);
		else image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		
		imageft = new BufferedImage(ftsizemin, ftsizemin, BufferedImage.TYPE_INT_RGB);
		System.out.println("image");
		
//		JOptionPane.showMessageDialog(null, "Select the mask (to be rotatable).");
//		refMask = FileOps.openBin(dir);
		putSquareMask();
		
		spectrum = FieldOps.getCorrelationvsTheta(masterField, refMask, npts, thmax);
		
		spec = new GraphDrawerCart("Spectrum", angles, spectrum);
		twoV = new double[] {0, 0};
		spec.setXY(new GraphDrawerCart.GraphObject(twoV, new double[] {ArrayOps.min(spectrum), ArrayOps.max(spectrum)}), false, false, 1);
		spec.showWindow();
		
		
		ftmax = ArrayOps.max(ftmag);
		ftmin = ArrayOps.min(ftmag);
		resetColorScale(0, 1);
		resetFTColorScale(0, 1);
		formFTImage();
//		resetFRColorScale(0, 1);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		showWindow();
		setTitle("FFT region: " + boxsize + "x" + boxsize);

		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		addMouseListener(this);
		JFrame f = new JFrame();
		BoundsPanel panel = new BoundsPanel(this, f);
		f.add(panel);
		f.setSize(200, 80);
		f.setVisible(true);
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
		formImage();
		repaint();
	}
	public void setFieldInfo()
	{
		fmax = max(field); fmin = min(field);
		fdelta = fmax - fmin;
		System.out.println("[" + fmin + ", " + fmax + "]");
	}
	public void redrawColorScale(double downnum, double upnum)
	{
		resetFTColorScale(downnum, upnum);
		formFTImage();
		repaint();
	}
	public void resetColorScale(double downnum, double upnum)
	{
		scale = new ColorScales.LinearBRYW(fmin + fdelta*upnum, fmin + fdelta*downnum);
	}
	public void resetFTColorScale(double downnum, double upnum)
	{
		upfrac = upnum;
		downfrac = downnum;
		ftscalemin = ftmin + (ftmax - ftmin)*downnum;
		ftscalemax = ftmin + (ftmax - ftmin)*upnum;
		ftscale = new ColorScales.LinearBRYW(ftscalemax, ftscalemin);
	}
	public void resetFTColorScale()
	{
		ftscale = new ColorScales.LinearBRYW(ftscalemax, ftscalemin);
	}
	public void formImage()
	{
		for (int i = 0; i < image.getWidth(); i++)
			for (int j = 0; j < image.getHeight(); j++)
				image.setRGB(i, j, scale.of(field[i/sizeratio][j/sizeratio]).getRGB());
	}
	public void formFTImage()
	{
		for (int i = 0; i < imageft.getWidth(); i++)
			for (int j = 0; j < imageft.getHeight(); j++)
				imageft.setRGB(i, j, ftscale.of(ftmag[i/ftsizeratio][j/ftsizeratio]).getRGB());
	}
	public void resetImages()
	{
		formImage();
		formFTImage();
		ImageEditing.fuse(image, imageft, blendp, image);
	}
	

	public int[] getFieldPoint(int x, int y)
	{
		return new int [] {x - ox, y - oy};
	}
	public int[] getImagePoint(int x, int y)
	{
		return new int [] {(x - ox)/sizeratio, ((y - oy)/sizeratio)};
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
		g.drawImage(image, ox, oy, null);
	}
	public void drawFTImage(Graphics g)
	{
		g.drawImage(imageft, ftox, ftoy, null);
	}
	
	public void paint(Graphics g)
	{
		g.clearRect(0, 0, 2000, 2000);
		drawImage(g);
		drawFTImage(g);
//		g.drawString(scale.toString(), writepoint[0], writepoint[1]);
//		g.drawString(ftscale.toString(), ftwritepoint[0], ftwritepoint[1]);
		g.drawString("Correlation is " + corr + "; Angle is " + Math.toDegrees(theta) + " degrees.", ftwritepoint[0], ftwritepoint[1]);
		
		g.setColor(java.awt.Color.BLUE);
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
		System.out.println(arg0.getKeyChar());
		
		boolean sqChanged = false;
		
		if (arg0.getKeyChar() == 'b') {width++; sqChanged = true;} 
		if (arg0.getKeyChar() == 'r') {width--; sqChanged = true;}
		if (arg0.getKeyChar() == 'w') {squareCent[1]--; sqChanged = true;}
		if (arg0.getKeyChar() == 'a') {squareCent[0]--; sqChanged = true;}
		if (arg0.getKeyChar() == 's') {squareCent[1]++; sqChanged = true;}
		if (arg0.getKeyChar() == 'd') {squareCent[0]++; sqChanged = true;}
		if (arg0.getKeyChar() == 't') {sqThickness++; sqChanged = true;}
		if (arg0.getKeyChar() == 'l') {sqThickness--; sqChanged = true;}
		if (arg0.getKeyChar() == 'm') {maskType++; maskType %= nMasks; sqChanged = true;}
		if (arg0.getKeyChar() == 'c') {
			resetSpectrum();
		}
		if (arg0.getKeyChar() == 'g'){
			squareCent[0] = inField[0];
			squareCent[1] = inField[1];
			sqChanged = true;
		}
		if (sqChanged)
		{
			changeSq();
			sqChanged = false;
		}
		
		if (arg0.getKeyChar() == ' ')
		{
			resetSpectrum();
			correlations[currentImp] = FieldOps.copy(spectrum);
			maxCorrAngles[currentImp] = angles[ArrayOps.maxIndex(correlations[currentImp])];
			imps.get(currentImp).position = this.getSquareCentInMap();
			imps.get(currentImp).angle = Math.toDegrees(maxCorrAngles[currentImp]);
			imps.get(currentImp).maxCorr = ArrayOps.max(correlations[currentImp]);
			imps.get(currentImp).sqSide = width;
			imps.get(currentImp).sqThick  = sqThickness;
			imps.get(currentImp).maskKey = maskType;
			
			System.out.println("SquareCent = (" + squareCent[0] + ", " + squareCent[1] + "); width = " + width + ", wall thickness = " + sqThickness);
			SRAW.writeImage(fc.getCurrentDirectory().toString() + "\\imppic" + MovieMaker.fromInt(currentImp), image);
			SRAW.writeImage(fc.getCurrentDirectory().toString() + "\\square" + MovieMaker.fromInt(currentImp), imageft);
			SRAW.writeImage(fc.getCurrentDirectory().toString() + "\\graph" + MovieMaker.fromInt(currentImp), (BufferedImage)spec.dbimage);
			currentImp++;
			
			if (currentImp == imps.size())
			{
//				Printer.printlnVertical(maxCorrAngles);
				
				JOptionPane.showMessageDialog(null, "Save the file with the correlation vs angle data.");
				FileOps.writeTableASCII(fc, correlations);
				
				JOptionPane.showMessageDialog(null, "Save the list of impurities. (includes updated positions, maximum angles, etc.");
				FileOps.writeString(fc, getImpText(imps));
				System.exit(0);
			}
			setUpForCurrentImp();
			setFieldInfo();
			resetColorScale(0, 1);
			resetImages();
			resetSpectrum();
			repaint();
		}
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent arg0) {
			
	}

	@Override
	public void mouseDragged(MouseEvent arg0) {
		// TODO Auto-generated method stub
		int x = arg0.getX();
		int y = arg0.getY();
//		System.out.println(x);
		currentx = x; currenty = y;
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		inField = this.getImagePoint(currentx, currenty);
		setTitle(titleString());
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
		
	}

	@Override
	public void mouseReleased(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	public void changeSq()
	{
		putSquareMask();
		resetImages();
		repaint();
//		FieldOps.putCorrelationvsTheta(masterField, refMask, thmax, spectrum);
//		setSpectrum();
//		spec.plot[1].y[0] = ArrayOps.min(spectrum);
//		spec.plot[1].y[1] = ArrayOps.max(spectrum);
//		spec.resetRange(0);
//		spec.repaint();
		setTitle(titleString());
	}
	public void resetSpectrum()
	{
		setSpectrum();
		spec.plot[1].y[0] = ArrayOps.min(spectrum);
		spec.plot[1].y[1] = ArrayOps.max(spectrum);
		spec.resetRange(0);
		spec.repaint();
	}
	public void setUpForCurrentImp()
	{
		int ix, iy;
		ix = (int)Math.round(imps.get(currentImp).position[0]);
		iy = (int)Math.round(imps.get(currentImp).position[1]);
		imageCorner[0] = ix-8;
		imageCorner[1] = iy-8;
		
//		int k = Integer.parseInt(JOptionPane.showInputDialog("Enter the index of the layer in which this impurity is the most visible."));
//		if (k >= t.nlayers) k = t.nlayers-1;
		int k = bestLayers[currentImp];
		FieldOps.expandBi(t.data[k], 8, ix-8, 16, iy-8, 16, masterField);
	}
	public String titleString()
	{
		return "Now at: (" + inField[0] + ", " + inField[1] + "); The center of the square is (" + squareCent[0] + ", " + squareCent[1] + "). The side length is " + width +". The thickness is " + sqThickness + ". The current impurity index is " + currentImp + ".";
	}
	public double[] getSquareCentInMap()
	{
		return new double[] {squareCent[0]/8 + imageCorner[0], squareCent[1]/8 + imageCorner[1]};
	}
	public void putSquareMask()
	{
		switch (maskType)
		{
		case 0:
			Mask.putGaussSidedSquareMaskRotated(N, N, width, squareCent, theta, sqThickness, ftmag, refMask, matrix);
			break;
		case 1:
			Mask.putSolidGaussSidedSquareMaskRotated(N, N, width, squareCent, theta, sqThickness, ftmag, refMask, matrix);
			break;
		case 2:
			Mask.putRectMaskRotated(N, N, width, width, squareCent, theta, refMask, ftmag, matrix);
			break;
		case 3:
			AltMask.putRectMaskRotated(N, N, width, squareCent, theta, ftmag, matrix);
			break;
		}
	}
	public static double max(double[][] x)
	{
		return ArrayOps.max(x);
	}
	public static double min(double[][] x)
	{
		return ArrayOps.min(x);
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
	 public double[][] getField() {
			return field;
		}
	 public BufferedImage getImage() {
			return image;
		}
	 public void setSpectrum()
	 {
		 tempMask = FieldOps.copy(ftmag);
		 tempTheta = theta;
		 for (int i = 0; i < npts; i++)
		 {
			 theta = angles[i];
			 putSquareMask();
			 spectrum[i] = FieldOps.correlation(field, ftmag);
		 }
		 ftmag = FieldOps.copy(tempMask);
		 theta = tempTheta;
	 }
	 
	 public static class BoundsPanel extends JPanel implements ChangeListener
		{
		private static final long serialVersionUID = 6507617978766359288L;
			GaussSquareImpurityAngleSuite parent;
			JFrame frame;
			public JSlider min, max, theta, blend;
			public BoundsPanel(GaussSquareImpurityAngleSuite parent, JFrame frame)
			{
				this.frame = frame;
				this.parent = parent;
				this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
				this.setBorder(new LineBorder(Color.GRAY));
				min = new JSlider(0, 1000, 0);
				min.addChangeListener(this);
				max = new JSlider(0, 1000, 1000);
				max.addChangeListener(this);
				theta = new JSlider(0, 1000, 0);
				theta.addChangeListener(this);
				blend = new JSlider(0, 1000, 0);
				blend.addChangeListener(this);
				
				add(blend);
				add(theta);
				add(min);
				add(max);
			}
			@Override
			public void stateChanged(ChangeEvent arg0) {
				if (theta.getValueIsAdjusting())
				{
					parent.theta = (parent.thmax/1000)*theta.getValue();
					parent.putSquareMask();
					parent.repaint();
					
					parent.twoV[0] = parent.theta;
					parent.twoV[1] = parent.theta;
					parent.resetImages();
					parent.spec.repaint();

				}
				if (blend.getValueIsAdjusting())
				{
					parent.blendp = blend.getValue()/1000.0;
					parent.resetImages();
					parent.repaint();
					parent.spec.repaint();
				}
				else{
					parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
					parent.resetImages();
					parent.repaint();
//					frame.setTitle("Min: " + parent.ftscalemin + "   Max: " + parent.ftscalemax);
				}
			}
		}
		public static ArrayList<Impurity> getImpurities(File f) {
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
				imps.add(new Impurity(Double.parseDouble(words[0]), Double.parseDouble(words[1])));
			}
			return imps;
		}
		public static ArrayList<Impurity> getImpuritiesFull(File f) {
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
				imps.add(new Impurity(line));
			}
			return imps;
		}
		public static String getImpText(ArrayList<Impurity> imps)
		{
			String s = "Impurity List (with angles)\r\n"
			 +"Mask ID code: 0 is hollow square with blurred Gaussian sides; 1 is solid square with blurred sides; 2 is sharp solid square; 3 is sharp solid square with no anti-aliasing.\r\n"
			 + "X \t Y \t Angle(deg.) \t Max. Corr. \t Square mask side length(*8) \t Gaussian side blur thickness(*8) \t Mask ID code \t Lattice X \t Lattice Y \r\n";
			for (int i = 0; i < imps.size(); i++)
				s += imps.get(i).toString() + "\r\n";
			return s;
		}

	 public static class Impurity
	 {
		 public double[] position;
		 public double[] latticePos;
		 public double angle;
		 public double maxCorr;
		 
		 //properties of the mask:
		 public int sqSide, maskKey;
		 public double sqThick;
		 
		 public Impurity(double x, double y, double angle, double maxCorr, int sqSide, double sqThick, int maskKey)
		 {
			 position = new double[] {x, y};
			 latticePos = new double[] {0, 0}; //blank for default.
			 this.angle = angle;
			 this.maxCorr = maxCorr;
			 this.sqSide = sqSide;
			 this.sqThick = sqThick;
			 this.maskKey = maskKey;
		 }
		 public Impurity(String line)
		 {
			 String[] words = line.split("\t");
			double x, y;
			 x = Double.parseDouble(words[0]);
			y =  Double.parseDouble(words[1]);
			position = new double[] {x, y};
			latticePos = new double[] {0, 0};
			angle =  Double.parseDouble(words[2]);
			maxCorr = Double.parseDouble(words[3]);
			sqSide = Integer.parseInt(words[4]);
			sqThick = Double.parseDouble(words[5]);
			maskKey = Integer.parseInt(words[6]);
			x =Double.parseDouble(words[7]);
			y = Double.parseDouble(words[8]);
			latticePos = new double[] {x, y};
		 }
		 public Impurity(double x, double y)
		 {
			 position = new double[] {x, y};
			 this.latticePos = new double [] {0, 0};
		 }
		 public String toString()
		 {
			 String x = "";
			 
			 x += position[0] + "\t" + position[1] + "\t" + angle + "\t" + maxCorr + "\t" + sqSide + "\t" + sqThick + "\t" + maskKey +"\t" + latticePos[0] + "\t" + latticePos[1];
			 return x;
		 }
	 }

	 public static void main(String[] args)
	 {
//	 	String dir = "C:\\data\\analysis\\";
//		JFileChooser fc = new JFileChooser(dir);
//		File f = null;
//		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
//			f = fc.getSelectedFile();
		new GaussSquareImpurityAngleSuite(Topomap.open(), Topomap.stddir);
	 }

}

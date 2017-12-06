package util.fourier;

import image.ImageEditing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Image;
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

import util.color.ColorScale;
import util.color.ColorScaleHolder;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.Topomap;
import util.geom.AtomicCoordinatesSet;
import util.*;
//This class only draws real scalar fields.

//Superseded: This class will soon be able to draw two-component fields in principle.
//It may draw a two-component field consisting of amplitude and phase data, under the assumption that
//the phase is periodic with period 2pi.

public class ImpurityListEditor extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{

	String dir;
	double[][] masterField; int N;
	double field[][]; double fmax, fmin, fdelta;
	boolean twoComp = false;
	BufferedImage image, imageft;
	BufferedImage imageConst;
	BufferedImage imagepink;

	ColorScale scale, ftscale;
	ColorScaleHolder csh;
	
	int sx = 1024, sy = 1024;
	int ox = 20, oy = 40;
	int ftox = ox + sx + 10, ftoy = oy;
	
	int ftsizemin = 1024;
	int ftsizeratio = 1;
	int ftsize = ftsizemin; //this is the size of the sub-window in pixels.
	
	int[] writepoint = {ox, oy + sy + 10};
	
	int zoomLevel = 0, zoomFactor = 1;
	
	int WIDTH = 1920, HEIGHT = 1080;
	
	Box box;
	Subset sub;
	ArrayList<Impurity> imps;
	double[][] ftmag;
	int[] ftwritepoint = {ftox, ftoy};
	int boxsize;
	
	//For ring, etc. generation
	double[][] pinkrings, pinkringssmall;
	double[][] greenring, greenringsmall;
	double[][] cyanline, cyanlinesmall;
	int nearestimpindex, oldnearestimpindex;
	
	double ringradius = 15, ringthickness = 1;
	double ringrsmall, ringtsmall;
	
	double ftmax, ftmin;
	double upfrac = 1, downfrac = 0;
	double ftscalemin, ftscalemax;
	boolean listChanged = false;
	
	double[][][] latticeSites = null; //2D array of lattice sites
	AtomicCoordinatesSet latt = null;
	boolean drawCircles = true;
	boolean drawLines = true;

	
	int currentx, currenty;
	int ix = 0, iy = 0; //Coordinates in the image
	
	int numdeletions = 0, numadditions = 0, numchanges = 0;
	int initialsize = 0;
	
	//coordinates of the mouse in the main image, based on its location in the box.
	double mxd, myd;
	int mxi, myi;
	
	//For double buffering
	Image dbimage;
	Graphics dbg;
	
	JFileChooser fc;
	//The field must be 2^n by 2^n;
	public ImpurityListEditor(String dir)
	{
		
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		File f = null;
		JOptionPane.showMessageDialog(null, "Select the reference topograph.");
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
			f = fc.getSelectedFile();
		double[][] field = Layer.openFree(f).data;


		this.field = field;
		this.dir = dir;
		System.out.println(dir);
		masterField = field;
		N = masterField.length;
		boxsize = N;
		sx = N;
		while (sx > 512)
		{
			this.field = FieldOps.reduce(2, this.field);
			sx = this.field.length;
			zoomLevel++;
			zoomFactor *= 2;
		}
		sy = this.field[0].length;
		JOptionPane.showMessageDialog(null, "Select the impurity list file.");
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
		{
			f = fc.getSelectedFile();
			imps = getImpurities(f, zoomFactor);
		}
		else
		{
			imps = new ArrayList<Impurity>();
			imps.add(new Impurity(1,1,1));
		}
		

		JOptionPane.showMessageDialog(null, "Open the file containing the lattice info. (If there is none, just click cancel.)");
		String w = FileOps.openText(fc);
		if (w != null){ 
			latt = new AtomicCoordinatesSet(w);
			populateLatticeSites();
		}

		
		initialsize = imps.size();
		ringrsmall = ringradius/zoomFactor;
		ringtsmall = ringthickness/zoomFactor;
		
		pinkrings = new double [N][N];
		pinkringssmall = new double [sx][sy];
		greenring = new double [N][N];
		greenringsmall = new double [sx][sy];
		
		ftox = ox + sx + 10; ftoy = oy;
		writepoint = new int[] {ox, oy + sy + 10};
		setFieldInfo();
		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		imageConst = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		imagepink = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		imageft = new BufferedImage(N, N, BufferedImage.TYPE_INT_RGB);
		box = new Box(0, 0, N);
		System.out.println("image");
		sub = new Subset(0, 0, N, N, masterField);
		csh = new ColorScaleHolder(FieldOps.min(field), FieldOps.max(field));

		ftmag = sub.values;
		ftmin = FieldOps.min(ftmag);
		ftmax = FieldOps.max(ftmag);
		fmin = FieldOps.min(masterField);
		fmax = FieldOps.max(masterField);
		fdelta = fmax-fmin;
		doSub();
		resetColorScale(0, 1);
		resetFTColorScale(0, 1);
		formFTImage();
//		resetFRColorScale(0, 1);
		formImageConst();
		makeSmallPinkRings();
		drawPinkRingsSmall();
//		formFTImage();
		//gradcalc.activate();
		showWindow();
		setTitle("FFT region: " + boxsize + "x" + boxsize);

		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		addMouseListener(this);
		JFrame frame = new JFrame();
		BoundsPanel panel = new BoundsPanel(this, frame);
		frame.add(panel);
		frame.setSize(200, 80);
		frame.show();
	}
	 public void populateLatticeSites()
	 {
		 double[] atomc1 = latt.getAtomicCoords(0, 0);
		 double[] atomc2 = latt.getAtomicCoords(0, N);
		 double[] atomc3 = latt.getAtomicCoords(N, 0);
		 double[] atomc4 = latt.getAtomicCoords(N, N);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 latticeSites = new double [xmx-xmn][ymx-ymn][2];
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j] = latt.getPixelCoords(i+xmn, j+ymn);
	 }

	private void loadNewImpFile()
	{
		File f;
		JOptionPane.showMessageDialog(null, "Select the impurity list file.");
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
		{
			removeGreenRing(nearestimpindex);
			f = fc.getSelectedFile();
			imps = getImpurities(fc, zoomFactor);
			initialsize = imps.size();
			numdeletions = 0;
			numadditions = 0;
			numchanges = 0;
			makeSmallPinkRings();
			drawPinkRingsSmall();
			repaint();
		}
		
	}
	public static ArrayList<Impurity> getImpurities(JFileChooser fc, int zoomFactor) {
		ArrayList<Impurity> imps = new ArrayList<Impurity>();
		Scanner in = null;
		String line;
		String[] words;
		File f = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f = fc.getSelectedFile();
		
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
			imps.add(new Impurity(Double.parseDouble(words[0]), Double.parseDouble(words[1]), zoomFactor));
		}
		return imps;
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
			imps.add(new Impurity(Double.parseDouble(words[0]), Double.parseDouble(words[1]), zoomFactor));
		}
		return imps;
	}
	public static ArrayList<Impurity> getImpurities(int zoomFactor) {
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		File f = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f = fc.getSelectedFile();
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
			imps.add(new Impurity(Double.parseDouble(words[0]), Double.parseDouble(words[1]), zoomFactor));
		}
		return imps;
	}
	private void makeSmallPinkRings()
	{
		double rangemax = ringrsmall + 4*ringtsmall;
		double g, r;
		for (int i = 0; i < sx; i++){
			for (int j = 0; j < sy; j++){ pinkringssmall[i][j] = 0;
				for (int k = 0; k < imps.size(); k++)
				{
					if (Math.abs(i - imps.get(k).possmall[0]) < rangemax && Math.abs(j - imps.get(k).possmall[1]) < rangemax)
					{
						r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
						g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
						pinkringssmall[i][j] += g;
					}
				}
				pinkringssmall[i][j] = Math.min(pinkringssmall[i][j], 1);
			}
		}
	}
	private void removePinkRing(int k)
	{
		double rangemax = ringrsmall + 4*ringtsmall;
		double g, r;
		int imin, imax, jmin, jmax;
		imin = (int)Math.max(imps.get(k).possmall[0]-rangemax, 0);
		jmin = (int)Math.max(imps.get(k).possmall[1]-rangemax, 0);
		imax = (int)Math.min(imps.get(k).possmall[0]+rangemax, sx-1);
		jmax = (int)Math.min(imps.get(k).possmall[1]+rangemax, sy-1);
		for (int i = imin; i < imax; i++){
			for (int j = jmin; j < jmax; j++){
				r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
				g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
				pinkringssmall[i][j] -= g;
				pinkringssmall[i][j] = Math.max(pinkringssmall[i][j], 0);
			}
		}
	}
	private void addPinkRing(int k)
	{
		double rangemax = ringrsmall + 4*ringtsmall;
		double g, r;
		int imin, imax, jmin, jmax;
		imin = (int)Math.max(imps.get(k).possmall[0]-rangemax, 0);
		jmin = (int)Math.max(imps.get(k).possmall[1]-rangemax, 0);
		imax = (int)Math.min(imps.get(k).possmall[0]+rangemax, sx-1);
		jmax = (int)Math.min(imps.get(k).possmall[1]+rangemax, sy-1);
		for (int i = imin; i < imax; i++){
			for (int j = jmin; j < jmax; j++){
				r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
				g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
				pinkringssmall[i][j] += g;
				pinkringssmall[i][j] = Math.min(pinkringssmall[i][j], 1);
			}
		}
	}
	private void addGreenRing(int k)
	{
		double rangemax = ringrsmall + 4*ringtsmall;
		double g, r;
		int imin, imax, jmin, jmax;
		imin = (int)Math.max(imps.get(k).possmall[0]-rangemax, 0);
		jmin = (int)Math.max(imps.get(k).possmall[1]-rangemax, 0);
		imax = (int)Math.min(imps.get(k).possmall[0]+rangemax, sx-1);
		jmax = (int)Math.min(imps.get(k).possmall[1]+rangemax, sy-1);
		for (int i = imin; i < imax; i++){
			for (int j = jmin; j < jmax; j++){
				r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
				g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
				greenringsmall[i][j] += g;
				greenringsmall[i][j] = Math.min(greenringsmall[i][j], 1);
			}
		}
	}
	private void removeGreenRing(int k)
	{
		double rangemax = ringrsmall + 4*ringtsmall;
		double g, r;
		int imin, imax, jmin, jmax;
		imin = (int)Math.max(imps.get(k).possmall[0]-rangemax, 0);
		jmin = (int)Math.max(imps.get(k).possmall[1]-rangemax, 0);
		imax = (int)Math.min(imps.get(k).possmall[0]+rangemax, sx-1);
		jmax = (int)Math.min(imps.get(k).possmall[1]+rangemax, sy-1);
		for (int i = imin; i < imax; i++){
			for (int j = jmin; j < jmax; j++){
				r = Math.sqrt(Math.pow(i - imps.get(k).possmall[0], 2) + Math.pow(j - imps.get(k).possmall[1], 2));
				g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
				greenringsmall[i][j] -= g;
				greenringsmall[i][j] = Math.max(greenringsmall[i][j], 0);
			}
		}
	}
	private void drawPinkRingsSmall()
	{
		ImageEditing.moveEachPixelTowardAColor(imageConst, pinkringssmall, Color.MAGENTA, imagepink);
	}
	private void drawGreenRingSmall()
	{
		ImageEditing.moveEachPixelTowardAColor(imagepink, greenringsmall, Color.GREEN, image);
	}
	public void selectNearestImp(int mx, int my)
	{
		oldnearestimpindex = nearestimpindex;
		double r = Double.MAX_VALUE, oldr = r;
		for (int i = 0; i < imps.size(); i++)
		{
			r = Math.sqrt(Math.pow(mx - imps.get(i).possmall[0], 2) + Math.pow(my - imps.get(i).possmall[1], 2));
			if (r < oldr)
			{
				oldr = r;
				nearestimpindex = i;
			}
		}
	}
	private void setImpsInBox(){
		for (int i = 0; i < imps.size(); i++)
			imps.get(i).isinbox = imps.get(i).position[0] > box.x && imps.get(i).position[1] > box.y
				&& imps.get(i).position[0] < box.x+box.size && imps.get(i).position[1] < box.y + box.size;
	}
	
	
	public void doSub()
	{
		sub.setCorner(box.x, box.y);
		sub.doSub();
		ftmag = sub.values;
//		ftmax = ArrayOps.max(ftmag);
//		ftmin = ArrayOps.min(ftmag);
//		resetFTColorScale(0, 1);
		formFTImage();
	}
	
	
	public void resizeFT(int size)
	{
		sub = new Subset(box.x, box.y, size, size, masterField);
		boxsize = size;
		box.setSize(size);
		if (size <= ftsizemin)
			ftsizeratio = ftsizemin/size;
		if (size*ftsizeratio == ftsizemin)
			imageft = new BufferedImage(Math.max(size, ftsizemin), Math.max(size, ftsizemin), BufferedImage.TYPE_INT_RGB);
		else
			imageft = new BufferedImage(size*ftsizeratio, size*ftsizeratio, BufferedImage.TYPE_INT_RGB);
		
		ftsize = imageft.getWidth();
		ftmag = new double[size][size];
		setTitle("FFT region: " + size + "x" + size + ";    " + "Expansion ratio " + ftsizeratio);
		setImpsInBox();
		doSub();
		repaint();
	}
	public void removeImpurity(int index)
	{
		listChanged = true;
		numdeletions++;
		removePinkRing(index);
		removeGreenRing(index);
		imps.remove(index);
		drawPinkRingsSmall();
		drawGreenRingSmall();
		repaint();
	}
	public void addImpurity()
	{
		listChanged = true;
		numadditions++;
		removeGreenRing(nearestimpindex);
		if (mouseInFTImage())
			imps.add(new Impurity(zoomFactor*(mxd-ox), zoomFactor*(myd-oy), zoomFactor));
		else
			imps.add(new Impurity(zoomFactor*ix, zoomFactor*iy, zoomFactor));
		addPinkRing(imps.size()-1);
		setImpsInBox();
//		addGreenRing(imps.size()-1);
		drawPinkRingsSmall();
		drawGreenRingSmall();
		repaint();
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
		csh.resetColorScaleChange(ftscalemin, ftscalemax, false);
		
	}
	public void moveBox(int dx, int dy)
	{
		box.move(box.x + dx, box.y + dy, this);
		setImpsInBox();
		doSub();
		repaint();
	}
	public void formImage()
	{
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				image.setRGB(i, j, scale.of(field[i][j]).getRGB());
	}
	public void formImageConst()
	{
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				imageConst.setRGB(i, j, scale.of(field[i][j]).getRGB());
	}
	public void formFTImage()
	{
		for (int i = 0; i < imageft.getWidth(); i++)
			for (int j = 0; j < imageft.getHeight(); j++)
				imageft.setRGB(i, j, csh.getScale().of(ftmag[i/ftsizeratio][j/ftsizeratio]).getRGB());
	}

	public void drawImage(Graphics g)
	{
		g.drawImage(imageConst, ox, oy, null);
	}
	public void drawFTImage(Graphics g)
	{
		g.drawImage(imageft, ftox, ftoy, null);
	}
	
	public void paint(Graphics g)
	{
		dbg.clearRect(0, 0, 2000, 2000);
		drawImage(dbg);
		drawFTImage(dbg);
		dbg.drawString(scale.toString(), writepoint[0], writepoint[1]);
		dbg.drawString(csh.getScale().toString(), ftwritepoint[0], ftwritepoint[1]);
		
		dbg.setColor(java.awt.Color.BLUE);
		box.draw(dbg, this);
		
		if (mouseInFTImage())
		{
			drawPlus(dbg, mxi, myi);
		}
		if (latt != null) drawLittleCirclesAndLines(dbg);
		drawRings(dbg);
		g.drawImage(dbimage, 0, 0, this);
		
	}
	public void update(Graphics g)
	{
		paint(g);
	}
	public void drawPlus(Graphics g, int x, int y)
	{
		g.drawLine(x-5, y, x+5, y);
		g.drawLine(x, y-5, x, y+5);
	}
	public void drawRings(Graphics g)
	{
		double rx, ry;
		int r = (int)(ringradius*ftsizeratio);
		g.setColor(Color.MAGENTA);
		for (int i = 0; i < imps.size(); i++)
		{
			if (imps.get(i).isinbox)
			{
				rx = FTX(imps.get(i).possmall[0]);
				ry = FTY(imps.get(i).possmall[1]);
				g.drawOval((int)(rx+0.5)-r, (int)(ry+0.5)-r, 2*r, 2*r);
			}
		}
		//also the green ring:
		if (imps.get(nearestimpindex).isinbox)
		{
			rx = FTX(imps.get(nearestimpindex).possmall[0]);
			ry = FTY(imps.get(nearestimpindex).possmall[1]);
			g.setColor(Color.green);
			g.drawOval((int)(rx+0.5)-r, (int)(ry+0.5)-r, 2*r, 2*r);
			
		}
	}
	
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		setVisible(true);
		dbimage = createImage(WIDTH, HEIGHT);
		dbg = dbimage.getGraphics();
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		repaint();
		show();
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
		if (arg0.getKeyChar() == 'b' || arg0.getKeyChar() == 'B' || arg0.getKeyCode() == KeyEvent.VK_PLUS)
			resizeFT(boxsize*2);
		if (arg0.getKeyChar() == 'r' || arg0.getKeyCode() == KeyEvent.VK_MINUS)
			resizeFT(boxsize/2);
		if (arg0.getKeyChar() == 'n' || arg0.getKeyChar() == 'N')
			addImpurity();
		if (arg0.getKeyChar() == 'a')
			moveBox(-1, 0);
		if (arg0.getKeyChar() == 'd')
			moveBox(+1, 0);
		if (arg0.getKeyChar() == 'w')
			moveBox(0, -1);
		if (arg0.getKeyChar() == 's')
			moveBox(0, +1);
		if (arg0.getKeyChar() == 'A')
			moveBox(-box.size/2, 0);
		if (arg0.getKeyChar() == 'D')
			moveBox(+box.size/2, 0);
		if (arg0.getKeyChar() == 'W')
			moveBox(0, -box.size/2);
		if (arg0.getKeyChar() == 'S')
			moveBox(0, +box.size/2);
		if (arg0.getKeyChar() == 'R')
			ringradius = Integer.parseInt(JOptionPane.showInputDialog("Enter the ring radius in pixels (integer)", "" + ringradius));
		if (arg0.getKeyChar() == 'x' || arg0.getKeyChar() == 'X')
			resizeFT(Integer.parseInt(JOptionPane.showInputDialog(null, "Enter the desired box size.")));
		if (arg0.getKeyChar() == KeyEvent.VK_DELETE)
			removeImpurity(nearestimpindex);
		if (arg0.getKeyChar() == 'O' || arg0.getKeyChar() == 'o')
			loadNewImpFile();
		if ((arg0.getKeyChar() == 'c' || arg0.getKeyChar() == 'C') && latt != null)
		{
			drawCircles = !drawCircles;
			repaint();
		}
		if ((arg0.getKeyChar() == 'l' || arg0.getKeyChar() == 'L') && latt != null)
		{
			drawLines = !drawLines;
			repaint();
		}
		if (arg0.getKeyChar() == ' ')
		{
			numchanges = numdeletions+numadditions;
			if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION){
				SRAW.writeImage(fc.getSelectedFile().toString()+"rings", imagepink);
				String name = fc.getSelectedFile().toString().substring(fc.getSelectedFile().toString().lastIndexOf("\\")+1);
//				SRAW.writeImage(fc.getSelectedFile().toString()+"rescaled", ftmag);
				String implist = getImpList(name);
				ColumnIO.writeString(implist, fc.getSelectedFile() + "list.txt");
				
				String message = "The original list had " + initialsize + " members.\r\n" +
						"There were " + numdeletions + " deletions and " + numadditions + " additions. (" + numchanges + " changes total)\r\n" +
						"The final list has " + imps.size() + " members.";
				JOptionPane.showMessageDialog(null, message);
			}
		}
		if (arg0.getKeyChar() == 'q')
		{
			csh.incrementScaleIndex();
			formImage();
			formFTImage();
			repaint();
		}
		if (arg0.getKeyChar() == 'e')
		{
			csh.decrementScaleIndex();
			formImage();
			formFTImage();
			repaint();
		}
		if (arg0.getKeyChar() == 'j')
			moveAllImps(-1.0/zoomFactor, 0);
		if (arg0.getKeyChar() == 'g')
			moveAllImps(1.0/zoomFactor, 0);
		if (arg0.getKeyChar() == 'y')
			moveAllImps(0, -1.0/zoomFactor);
		if (arg0.getKeyChar() == 'h')
			moveAllImps(0, 1.0/zoomFactor);
		
			

	}
	public void moveAllImps(double dx, double dy)
	{
		for (int i = 0; i < imps.size(); i++){
			imps.get(i).setPosition(imps.get(i).position[0]+dx, imps.get(i).position[1]+dy, zoomFactor);
//			if (imps.get(i).position[0] > N || imps.get(i).position[0] < 0 || imps.get(i).position[1] > N || imps.get(i).position[1] < 0)
//				imps.remove(i);
		}
		setImpsInBox();
//		addGreenRing(imps.size()-1);
		drawPinkRingsSmall();
		drawGreenRingSmall();
		repaint();

		
		
	}
	public String getImpList(String name)
	{
		String implist = "Impurity " + name + " position list.\r\n";
		implist += "(" + numchanges + " changes to the initial list.)\r\n";
		implist += "X \t Y \r\n";
		for (int i = 0; i < imps.size(); i++)
		{
			implist += imps.get(i).position[0] + "\t" + imps.get(i).position[1] + "\r\n";
		}
		return implist;
	}
	public boolean mouseInImage()
	{
		return ix >= 0 && iy >= 0 && ix < sx && iy < sy;
	}
	public boolean mouseInFTImage()
	{
		return currentx >= ftox && currenty >= ftoy && currentx < ftox+ftsize && currenty < ftoy + ftsize;
	}
	//returns the rounded x-coordinate of the point in the image corresponding to the mouse's location in ftimage.
	public double getMouseXfromFT()
	{
		return ( (((currentx - ftox)/(double)ftsize)*box.size/zoomFactor) + box.x/zoomFactor + ox);
	}
	public double getMouseYfromFT()
	{
		return ( (((currenty - ftoy)/(double)ftsize)*box.size/zoomFactor) + box.y/zoomFactor + oy);
	}
	//These give the location on screen of a point in the aux. window, given its coordinates in the main window and the location of the box.
	public double FTX(double imagex)
	{
		return (((imagex*zoomFactor) - box.x)/(double)box.size)*ftsize + ftox;
	}
	public double FTY(double imagey)
	{
		return (((imagey*zoomFactor) - box.y)/(double)box.size)*ftsize + ftoy;
	}
	public void drawLittleCirclesAndLines(Graphics g)
	{
		int cx, cy;
		int[][] fourNeighbors = new int [4][2];
		for (int i = 0; i < latticeSites.length; i++)
			for (int j = 0; j < latticeSites[0].length; j++)
			{
				if (box.contains(latticeSites[i][j]))
				{
					cx = xInFT(latticeSites[i][j][0]);
					cy = yInFT(latticeSites[i][j][1]);
					if (this.drawCircles) drawCircle(g, cx, cy, 4, Color.BLUE);
					//now draw the lines. Only draw once.
					if (this.drawLines && ((i%2==0 && j%2==0)||(i%2==1 && j%2==1)) && i >= 1 && i+1 < latticeSites.length && j >= 1 && j+1 < latticeSites[0].length)
					{
							fourNeighbors[0][0] = xInFT(latticeSites[i+1][j][0]);
							fourNeighbors[0][1] = yInFT(latticeSites[i+1][j][1]);
							fourNeighbors[1][0] = xInFT(latticeSites[i-1][j][0]);
							fourNeighbors[1][1] = yInFT(latticeSites[i-1][j][1]);
							fourNeighbors[2][0] = xInFT(latticeSites[i][j+1][0]);
							fourNeighbors[2][1] = yInFT(latticeSites[i][j+1][1]);
							fourNeighbors[3][0] = xInFT(latticeSites[i][j-1][0]);
							fourNeighbors[3][1] = yInFT(latticeSites[i][j-1][1]);
							g.setColor(Color.GREEN);
							for (int k = 0; k < 4; k++)
								g.drawLine(cx, cy, fourNeighbors[k][0], fourNeighbors[k][1]);
					}
				}
			}
	}
	public void drawCircle(Graphics g, int x, int y, int r, Color c)
	{
		g.setColor(c);
		g.drawOval(x-r, y-r, 2*r, 2*r);
	}
	
	public int xInFT(double xInMap)
	{
		return ftox + FieldOps.round(ftsizeratio*(xInMap-box.x));
	}
	public int yInFT(double yInMap)
	{
		return ftoy + FieldOps.round(ftsizeratio*(yInMap-box.y));
	}

	
	@Override
	public void mouseWheelMoved(MouseWheelEvent arg0) {}
	@Override
	public void mouseDragged(MouseEvent arg0) {
		// TODO Auto-generated method stub
		int x = arg0.getX();
		int y = arg0.getY();
		moveBox((x - currentx)*zoomFactor, (y - currenty)*zoomFactor);
//		System.out.println(x);
		currentx = x; currenty = y;
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		ix = currentx - ox;
		iy = currenty - oy;
		if (mouseInFTImage())
		{
			mxd = this.getMouseXfromFT();
			myd = this.getMouseYfromFT();
			mxi = (int)(mxd+0.5);
			myi = (int)(myd+0.5);
			selectNearestImp(mxi-ox,myi-oy);
			if (listChanged)
			{
				addGreenRing(nearestimpindex);
				drawGreenRingSmall();
			}
			else if (oldnearestimpindex != nearestimpindex)
			{
				removeGreenRing(oldnearestimpindex);
				addGreenRing(nearestimpindex);
				drawGreenRingSmall();
				setTitle("Current Imp: " + nearestimpindex);
			}
			repaint();
		}
		else
		{
			selectNearestImp(ix, iy);
			if (listChanged)
			{
				addGreenRing(nearestimpindex);
				drawGreenRingSmall();
				repaint();
			}
			else if (oldnearestimpindex != nearestimpindex)
			{
				removeGreenRing(oldnearestimpindex);
				addGreenRing(nearestimpindex);
				drawGreenRingSmall();
				repaint();
			}
		}
		listChanged = false;
	}

	@Override
	public void mouseClicked(MouseEvent arg0) {	}
	@Override
	public void mouseEntered(MouseEvent arg0) {	}
	@Override
	public void mouseExited(MouseEvent arg0) {	}
	@Override
	public void mousePressed(MouseEvent arg0) {	}
	@Override
	public void mouseReleased(MouseEvent arg0) {	}
	@Override
	public void actionPerformed(ActionEvent arg0) {	}

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
	 
	 private static class Box
	 {
		 int x, y, size;

		public Box(int x, int y, int size) {
			this.x = x;
			this.y = y;
			this.size = size;
		}
		
		public void draw(Graphics g, ImpurityListEditor parent)
		{
			g.drawRect(parent.ox + x/parent.zoomFactor, parent.oy + y/parent.zoomFactor, size/parent.zoomFactor, size/parent.zoomFactor);
		}
		
		public void move(int x, int y, ImpurityListEditor parent)
		{
			this.x = x > 0 ? x+size < parent.N ? x : parent.N - size : 0; this.y = y > 0 ? y+size < parent.N ? y : parent.N - size : 0;
		}
		public void setSize(int size)
		{
			this.size = size;
		}
		public boolean contains(double[] point)
		{
			return point[0] >= x && point[0] <= x+size && point[1] >= y && point[1] <= y+size;
		}

	 }
	 public static class Subset
	 {
		 int x, y, dy, dx; double[][] masterField;
		 double[][] values;
		 
		 public Subset(int x, int y, int dx, int dy, double[][] masterField)
		 {
			 this.x = x;
			 this.y = y;
			 this.dx = dx;
			 this.dy = dy;
			 this.masterField = masterField;
			 this.values = new double[dx][dy];
			 for (int i = 0; i < dx; i++)
				 for (int j = 0; j < dy; j++)
					 values[i][j] = masterField[x+i][y+j];
		 }
		 public void doSub()
		 {
			 for (int i = 0; i < dx; i++)
				 for (int j = 0; j < dy; j++)
					 values[i][j] = masterField[x+i][y+j];
		 }
		 public void setCorner(int x, int y)
		 {
			 this.x = x;
			 this.y = y;
			 
		 }
	 }
	 public static class Impurity
	 {
		 double[] position;
		 double[] possmall;
		 boolean isinbox = false;
		 
		 public Impurity(double x, double y, int zoomfactor)
		 {
			 position = new double[] {x, y};
			 possmall = new double[] {x/zoomfactor, y/zoomfactor};
		 }
		 public void setPosition(double x, double y, int zoomfactor)
		 {
			 position = new double[] {x, y};
			 possmall = new double[] {x/zoomfactor, y/zoomfactor};
		 }
		 public String toString()
		 {
			 return "Real position [" + position[0] + ", " + position[1] + "];  reduced [" + possmall[0] + ", " + possmall[1] + "]";  
		 }
		 public double[] getPosition()
		 {
			 return position;
		 }
	 }
	 public static class BoundsPanel extends JPanel implements ChangeListener
		{
			ImpurityListEditor parent;
			JFrame frame;
			public JSlider min, max;
			public BoundsPanel(ImpurityListEditor parent, JFrame frame)
			{
				this.frame = frame;
				this.parent = parent;
				this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
				this.setBorder(new LineBorder(Color.GRAY));
				min = new JSlider(0, 1000, 0);
				min.addChangeListener(this);
				max = new JSlider(0, 1000, 1000);
				max.addChangeListener(this);

				add(min);
				add(max);
			}
			@Override
			public void stateChanged(ChangeEvent arg0) {
				parent.redrawColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
				frame.setTitle("Min: " + parent.ftscalemin + "   Max: " + parent.ftscalemax);
			}
		}
	 public static class ListConverter
	 {
		 public static void convertPixelsToNM(double toposize, int npix, ArrayList<Impurity> imps)
		 {
			 //flip the y-coordinate also
			 double nmperpix = toposize/npix;
			 for (int i = 0; i < imps.size(); i++)
			 {
				 imps.get(i).position[0] *= nmperpix;
				 imps.get(i).position[1] *= nmperpix;
				 imps.get(i).position[1] = toposize - imps.get(i).position[1];
			 }
		 }
		 public static void convertLists(String dir)
		 {
			 	File f;
			 	ArrayList<Impurity> imps = new ArrayList<Impurity>();
			 	JFileChooser fc = new JFileChooser(dir);
			 	
			 	int npix = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of pixels."));
			 	double toposize = Double.parseDouble(JOptionPane.showInputDialog("Enter the size the topograph in nm."));
			 	
				JOptionPane.showMessageDialog(null, "Select the impurity list.");
				String name;
				while(fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
				{
					f = fc.getSelectedFile();
					name = f.toString().substring(f.toString().lastIndexOf("\\")+1, f.toString().length()-4);
					imps = getImpurities(f, 1);
					convertPixelsToNM(toposize, npix, imps);
					String implist = "Impurity " + name + " position list.\r\n";
					implist += "X \t Y \r\n";
					for (int i = 0; i < imps.size(); i++)
					{
						implist += imps.get(i).position[0] + "\t" + imps.get(i).position[1] + "\r\n";
					}
					if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
						ColumnIO.writeString(implist, fc.getSelectedFile().toString() + ".txt");
					JOptionPane.showMessageDialog(null, "Select the impurity list.");
				}

		 }
	 }
	 public static void main(String[] args)
	 {
	 	String dir = "C:\\Users\\Charlie\\Desktop\\spm data\\Practice\\processed data\\In-doped Bi2Se3\\1%\\";
		new ImpurityListEditor(dir);
//	 	ListConverter.convertLists(dir);
	 }

}

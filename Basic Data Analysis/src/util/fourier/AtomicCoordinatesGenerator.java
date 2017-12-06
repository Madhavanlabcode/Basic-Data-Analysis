package util.fourier;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.Topomap;
import util.geom.AtomicCoordinatesSet;
import util.*;
//This class only draws real scalar fields.

//Superseded: This class will soon be able to draw two-component fields in principle.
//It may draw a two-component field consisting of amplitude and phase data, under the assumption that
//the phase is periodic with period 2pi.

public class AtomicCoordinatesGenerator extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	
	static JFileChooser fc;
	String dir;
	double[][] masterField; int N;
	double field[][]; double fmax, fmin, fdelta;
	boolean twoComp = false;
	BufferedImage image, imageft;
	int linesize = 12;
	
	boolean justPaintText = false;
	ColorScale scale, ftscale;
	
	AtomicCoordinatesSet atom = null;
	String atomDesc;
	double[][] bragg = new double [2][2];
	double[][] latticeVectors = new double [2][2];
	double[][][] latticeSites = null; //2D array of lattice sites once the coordinate system shall have been chosen.

	int currentBragg = 0;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	int ftox = ox + sx + 10, ftoy = oy;
	
	int ftsx = 1024, ftsy = 1024;
	int ftsizeratio = 1;
	
	int[] writepoint = {ox, oy + sy + 10};
	
	int zoomLevel = 0, zoomFactor = 1;
	
	int WIDTH = 1920, HEIGHT = 1060;
	
	Box box;
	Subset sub;
	double[][] ftmag, ftmasterField;
	
	double[][] hiResFtmag = new double [ftsx][ftsy];
	boolean hiRes = false;
	int hiResRatio;
	int[] ftwritepoint = {ftox, ftoy};
	int boxsize;
	
	
	double ftmax, ftmin;
	double upfrac = 1, downfrac = 0;
	double ftscalemin, ftscalemax;
	
	int currentx, currenty;
	int currentift, currentjft;
	double currentxft, currentyft;
	double[] currentAtomft = new double [2];
	//The field must be 2^n by 2^n;
	public AtomicCoordinatesGenerator(double[][] field, String dir)
	{
		this.dir = dir;
		System.out.println(dir);
		ftsx = field.length; ftsy = field[0].length;
		masterField = field;
		N = masterField.length;
		boxsize = N;
		sx = N;
//		this.field = field;
		this.field = FFTOps.obtainFFTmagCent(masterField);
		FieldOps.log(this.field);
		ftmasterField = FieldOps.copy(this.field);
		while (sx > 512)
		{
			this.field = FieldOps.reduce(2, this.field);
			sx = this.field.length;
			zoomLevel++;
			zoomFactor *= 2;
		}
		int hasone = JOptionPane.showConfirmDialog(null, "Does a lattice file already exist?", "", JOptionPane.YES_NO_OPTION);
		if (hasone == JOptionPane.YES_OPTION){
			atom = new AtomicCoordinatesSet(FileOps.openText(fc));
			populateLatticeSites();
			currentBragg = 2;
			JOptionPane.showMessageDialog(null, "Press Spacebar.");
		}
		else
			JOptionPane.showMessageDialog(null, "Select the first bragg peak by zooming in, mousing over it, and then pressing the spacebar.");
		sy = this.field[0].length;
		ftox = ox + sx + 10; ftoy = oy;
		writepoint = new int[] {ftox + ftsx + 10, oy};
		setFieldInfo();
		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		imageft = new BufferedImage(N, N, BufferedImage.TYPE_INT_RGB);
		box = new Box(0, 0, field.length);
		System.out.println("image");
		sub = new Subset(0, 0, N, N, this.ftmasterField);
		
		ftmag = sub.values;
		doSub();
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
		f.show();
	}
	public void doSub()
	{	
		if (!hiRes)
		{	
			sub.setCorner(box.x, box.y);
			sub.doSub();
			ftmag = sub.values;
			ftmax = ArrayOps.max(ftmag);
			ftmin = ArrayOps.min(ftmag);
		}
		else
		{
			hiResRatio = ftsx/box.size;
			FieldOps.expandBi(masterField, hiResRatio, box.x, box.size, box.y, box.size, hiResFtmag);
			ftmax = ArrayOps.max(hiResFtmag);
			ftmin = ArrayOps.min(hiResFtmag);
		}
		justPaintText = false;
		resetFTColorScale();
		formFTImage();
	}
	
	
	public void resizeFT(int size)
	{
		if (currentBragg < 2)
			sub = new Subset(box.x, box.y, size, size, ftmasterField);
		else
			sub = new Subset(box.x, box.y, size, size, masterField);
		boxsize = size;
		box.setSize(size);
		if (size <= ftsx)
			ftsizeratio = ftsx/size;
		if (size*ftsizeratio == ftsx)
			imageft = new BufferedImage(Math.max(size, ftsx), Math.max(size, ftsx), BufferedImage.TYPE_INT_RGB);
		else
			imageft = new BufferedImage(size*ftsizeratio, size*ftsizeratio, BufferedImage.TYPE_INT_RGB);
			
		ftmag = new double[size][size];
		setTitle("FFT region: " + size + "x" + size + ";    " + "Expansion ratio " + ftsizeratio);
		doSub();
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
		ftscale = new ColorScales.LinearBRYW(ftscalemax, ftscalemin);
	}
	public void moveBox(int dx, int dy)
	{
		box.move(box.x + dx, box.y + dy, this);
		doSub();
		repaint();
	}
	public void moveAtomOrigin(int dx, int dy)
	{
		double dxd = dx/(double)ftsizeratio;
		double dyd = dy/(double)ftsizeratio;
		atom.moveOrigin(dxd, dyd);
		populateLatticeSites();
		justPaintText = false;
		repaint();
	}
	public void resetFTColorScale()
	{
		ftscale = new ColorScales.LinearBRYW(ftscalemax, ftscalemin);
	}
	public void formImage()
	{
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				image.setRGB(i, j, scale.of(field[i][j]).getRGB());
	}
	public void formFTImage()
	{
		if (!hiRes)
			for (int i = 0; i < imageft.getWidth(); i++)
				for (int j = 0; j < imageft.getHeight(); j++)
					imageft.setRGB(i, j, ftscale.of(ftmag[i/ftsizeratio][j/ftsizeratio]).getRGB());
		else
			for (int i = 0; i < imageft.getWidth(); i++)
				for (int j = 0; j < imageft.getHeight(); j++)
					imageft.setRGB(i, j, ftscale.of(hiResFtmag[i][j]).getRGB());
			
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
		g.drawImage(image, ox, oy, null);
	}
	public void drawFTImage(Graphics g)
	{
		g.drawImage(imageft, ftox, ftoy, null);
	}
	
	public void paint(Graphics g)
	{
		if (!justPaintText)
		{
			g.clearRect(0, 0, 2000, 2000);
			drawImage(g);
			drawFTImage(g);
			if (atom != null) drawLittleCirclesAndLines(g);

		}
		g.drawString(scale.toString(), writepoint[0], writepoint[1]);
		g.drawString(ftscale.toString(), ftwritepoint[0], ftwritepoint[1]);
		drawText(g);
		g.setColor(java.awt.Color.BLUE);
		box.draw(g, this);
	}
	public void drawText(Graphics g)
	{
		g.clearRect(writepoint[0], writepoint[1]-linesize, 500, sy);
		String it = "";
		it += "Current mouse position in map: (" + currentift + ", " + currentjft + ")" + "\r\n";
		it += "More precisely, (" + currentxft + ", " + currentyft + ") \r\n";
		it += "W.r.t. the center, (" + (currentift-ftsx/2) + ", " + (currentjft-ftsy/2) + ")\r\n";
		if (atom != null){
			it += "The position is " + Printer.vectorP(currentAtomft) + " in lattice coordinates.\r\n";
			it += atomDesc;
		}
		String[] lines = it.split("\r\n");
		for (int i = 0; i < lines.length; i++)
			g.drawString(lines[i], writepoint[0], writepoint[1] + (i+1)*linesize);
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
//					drawCircle(g, cx, cy, 4, Color.BLUE);
					//now draw the lines. Only draw once.
					if (((i%2==0 && j%2==0)||(i%2==1 && j%2==1)) && i >= 1 && i+1 < latticeSites.length && j >= 1 && j+1 < latticeSites[0].length)
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
	 public void  populateLatticeSites()
	 {
		 double[] atomc1 = atom.getAtomicCoords(0, 0);
		 double[] atomc2 = atom.getAtomicCoords(0, N);
		 double[] atomc3 = atom.getAtomicCoords(N, 0);
		 double[] atomc4 = atom.getAtomicCoords(N, N);
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
				 latticeSites[i][j] = atom.getPixelCoords(i+xmn, j+ymn);
	 }
	 
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
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
		if (arg0.getKeyChar() == 'r' || arg0.getKeyChar() == 'R' || arg0.getKeyCode() == KeyEvent.VK_MINUS)
			resizeFT(boxsize/2);
		if (arg0.getKeyChar() == 'n' || arg0.getKeyChar() == 'N')
			resizeFT(Integer.parseInt(JOptionPane.showInputDialog(null, "Enter the desired box size.")));
		if (arg0.getKeyChar() == 'a' || arg0.getKeyChar() == 'A' && atom != null)
			moveAtomOrigin(-1, 0);
		if (arg0.getKeyChar() == 'd' || arg0.getKeyChar() == 'D' && atom != null)
			moveAtomOrigin(+1, 0);
		if (arg0.getKeyChar() == 'w' || arg0.getKeyChar() == 'W' && atom != null)
			moveAtomOrigin(0, -1);
		if (arg0.getKeyChar() == 's' || arg0.getKeyChar() == 's' && atom != null)
			moveAtomOrigin(0, +1);
		if(arg0.getKeyChar() == 'u' || arg0.getKeyChar() == 'U')
		{
			hiRes = !hiRes;
			doSub();
			repaint();
		}
		if (arg0.getKeyChar() == ' ')
		{
			if (currentBragg < 2){
				bragg[currentBragg] = new double[] {currentift-N/2, currentjft-N/2};
				currentBragg++;
				if (currentBragg == 1)
					JOptionPane.showMessageDialog(null, "Now the other one.");
				else{
					JOptionPane.showMessageDialog(null, "Now pick the origin in the same way.");
					latticeVectors = new double [2][2];
					//so, actually, the lattice vector is perpindicular to the reciprocal lattice vector.
//					for (int i = 0; i < 2; i++)
//						for (int j = 0; j < 2; j++)
//						{
//							latticeVectors[i][j] = bragg[i][j];
//							latticeVectors[i][j] *= (N/(Complex.mag(bragg[i])*Complex.mag(bragg[i])));
//						}
					
					//the reciprocal lattice vector is expressed in pixels; we should therefore have Lattice dot Bragg = N.
					double dot;
					double m;
					for (int i = 0; i < 2; i++)
					{	//This is a clockwise rotation by 90 degrees with respect to the other bragg peak, given that the coordinate system is wrong-handed.
						latticeVectors[i][0] = -bragg[(i+1)%2][1];
						latticeVectors[i][1] = bragg[(i+1)%2][0];
						m = Complex.mag(latticeVectors[i]);
						for (int j = 0; j < 2; j++)
							latticeVectors[i][j] /= m; //which is the same as the magnitude of the Bragg peak of course.
						//Now we have unit vectors.
						dot = AtomicCoordinatesSet.dot(latticeVectors[i], bragg[i]);
						//the dot product should be N. Therefore, enlarge it by N/dot.
						for (int j = 0; j < 2; j++)
							latticeVectors[i][j] *= N/dot; 						
					}
					
					field = masterField;
					while (sx > 512)
					{
						this.field = FieldOps.reduce(2, this.field);
						sx = this.field.length;
						zoomLevel++;
						zoomFactor *= 2;
					}
					setFieldInfo();
//					sub.masterField = field;
					sub = new Subset(0, 0, N, N, masterField);
					moveBox(-box.x, -box.y);
					doSub();
					resetFTColorScale(0, 1);
					resetColorScale(0, 1);
					formImage();
					formFTImage();
					repaint();
				}
			}
			else if (currentBragg == 2)
			{
				currentBragg++;
				double[] atomO = new double[] {currentxft, currentyft};
				if (atom == null)
					atom = new AtomicCoordinatesSet(latticeVectors[0], latticeVectors[1], atomO);
				populateLatticeSites();
				JOptionPane.showMessageDialog(null, "To fine-tune the origin, use WASD. When satisfied, press Spacebar to save.");
			}
			else if (currentBragg == 3)
			{
				FileOps.writeString(fc, atom.toString());				
			}
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
		moveBox((x - currentx)*zoomFactor, (y - currenty)*zoomFactor);
//		System.out.println(x);
		currentx = x; currenty = y;
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		
		currentift = ((currentx - ftox)/ftsizeratio) + box.x;
		currentjft = ((currenty - ftoy)/ftsizeratio) + box.y;
		currentxft = ((currentx - ftox)/(double)ftsizeratio) + box.x;
		currentyft = ((currenty - ftoy)/(double)ftsizeratio) + box.y;
		if (atom != null) {
			atom.putAtomicCoords(currentxft, currentyft, currentAtomft);
			atomDesc = atom.forDisplay();
		}
		
		justPaintText = true;
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
		
		public void draw(Graphics g, AtomicCoordinatesGenerator parent)
		{
			g.drawRect(parent.ox + x/parent.zoomFactor, parent.oy + y/parent.zoomFactor, size/parent.zoomFactor, size/parent.zoomFactor);
		}
		
		public void move(int x, int y, AtomicCoordinatesGenerator parent)
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
	 public static class BoundsPanel extends JPanel implements ChangeListener
		{
			AtomicCoordinatesGenerator parent;
			JFrame frame;
			public JSlider min, max;
			public BoundsPanel(AtomicCoordinatesGenerator parent, JFrame frame)
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
				parent.justPaintText = false;
				parent.redrawColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
				frame.setTitle("Min: " + parent.ftscalemin + "   Max: " + parent.ftscalemax);
			}
		}
	 
	 //The offset is in lattice units
	 /**
	  * This returns a grid of lattice sites which includes the entire area of "real space" from 0 to N pixels.
	  * The lattice grid is a square table, so that in real space it would be slanted; the lattice area therefore includes points which are outside the real-space rectangle.
	  * @param a
	  * @param N
	  * @param offset
	  * @return
	  */
	 public static double[][][] getLatticeSites(AtomicCoordinatesSet a, int N, double[] offset)
	 {
		 if (offset == null) offset = new double [2];
		 double[] atomc1 = a.getAtomicCoords(0, 0);
		 double[] atomc2 = a.getAtomicCoords(0, N);
		 double[] atomc3 = a.getAtomicCoords(N, 0);
		 double[] atomc4 = a.getAtomicCoords(N, N);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 double[][][] latticeSites = new double [xmx-xmn][ymx-ymn][2];
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j] = a.getPixelCoords(i+xmn+offset[0], j+ymn+offset[1]);
		 return latticeSites;
	 }
	 public static void main(String[] args)
	 {
		 Topomap.setStdDir();
	 	String dir = Topomap.stddir;
		fc = new JFileChooser(dir);
		Layer l = Layer.openFree(fc);
		File f = fc.getSelectedFile();
		new AtomicCoordinatesGenerator(l.data, fc.getCurrentDirectory().toString());
	 }

}

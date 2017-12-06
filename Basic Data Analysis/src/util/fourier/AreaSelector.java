package util.fourier;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.*;
import java.awt.image.BufferedImage;

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



public class AreaSelector extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{

	String dir;
	double[][] masterField; int N;
	double field[][]; double fmax, fmin, fdelta;
	boolean twoComp = false;
	BufferedImage image, imageft;
	private int linesize = 12;
	
	boolean justPaintText = false;
	boolean drawCircles = true;
	boolean drawLines = true;
	
	AtomicCoordinatesSet latt = null;
	String atomDesc;
	ColorScale scale, ftscale;
	double[][][] latticeSites = null; //2D array of lattice sites
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	int ftox = ox + sx + 10, ftoy = oy;
	
	int ftsizemin = 1024;
	int ftsizeratio = 1;
	
	int[] writepoint = {ox + 256 + ftsizemin + 20, oy + sy + 10};
	
	int zoomLevel = 0, zoomFactor = 1;
	
	int WIDTH = 1920, HEIGHT = 1060;
	
	Box box;
	Subset sub;
	double[][] ftmag;

	
	double[][] hiResFtmag = new double [ftsizemin][ftsizemin];
	boolean hiRes = false;
	int hiResRatio;
	int[] ftwritepoint = {ftox, ftoy};
	int boxsize;
	
	double ftmax, ftmin;
	double upfrac = 1, downfrac = 0;
	double ftscalemin, ftscalemax;
	
	int currentx, currenty;
	private int currentift;
	private int currentjft;
	private double currentxft;
	private double currentyft;
	private double[] currentAtomft = null;
	public double[] currentPixft = null; //conversion back from atomic coordinates.
	//The field must be 2^n by 2^n;
	public AreaSelector(double[][] field, String dir)
	{
		this.dir = dir;
		System.out.println(dir);
		masterField = field;
		N = masterField.length;
		boxsize = N;
		sx = N;
		this.field = masterField;
		while (sx > 512)
		{
			this.field = FieldOps.reduce(2, this.field);
			sx = this.field.length;
			zoomLevel++;
			zoomFactor *= 2;
		}

		JOptionPane.showMessageDialog(null, "Open the file containing the lattice info. (If there is none, just click cancel.)");
		String w = FileOps.openText();
		if (w != null){ 
			latt = new AtomicCoordinatesSet(w);
			atomDesc = latt.forDisplay();
			currentAtomft = new double [2];
			currentPixft = new double [2];
			populateLatticeSites();
		}


		sy = this.field[0].length;
		ftox = ox + sx + 10; ftoy = oy;
		if (sx > 256)
			writepoint = new int[] {ox, oy + sy + 10};

		//		writepoint = new int[] {ox, oy + sy + 10};
		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		imageft = new BufferedImage(ftsizemin, ftsizemin, BufferedImage.TYPE_INT_RGB);
		box = new Box(0, 0, field.length);
		System.out.println("image");
		sub = new Subset(0, 0, N, N, masterField);
		ftsizeratio = ftsizemin/N;
		ftmag = sub.values;
		setFieldInfo();
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
		f.setVisible(true);
	}
	public AreaSelector(Layer l, String dir)
	{
		this(l.data, dir);
	}
	public void doSub()
	{	
		if (!hiRes)
		{	
			sub.setCorner(box.x, box.y);
			sub.doSub();
			ftmag = sub.values;
//			ftmax = ArrayOps.max(ftmag);
//			ftmin = ArrayOps.min(ftmag);
		}
		else
		{
			hiResRatio = ftsizemin/box.size;
			FieldOps.expandBi(masterField, hiResRatio, box.x, box.size, box.y, box.size, hiResFtmag);
//			ftmax = ArrayOps.max(hiResFtmag);
//			ftmin = ArrayOps.min(hiResFtmag);
		}
		resetFTColorScale();
		formFTImage();
	}
	
	
	public void resizeFT(int size)
	{
		sub = new Subset(box.x, box.y, size, size, masterField);
		justPaintText = false;
		boxsize = size;
		box.setSize(size);
		if (size <= ftsizemin)
			ftsizeratio = ftsizemin/size;
		if (size*ftsizeratio == ftsizemin)
			imageft = new BufferedImage(Math.max(size, ftsizemin), Math.max(size, ftsizemin), BufferedImage.TYPE_INT_RGB);
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
		ftmax = max(ftmag);
		ftmin = min(ftmag);
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
		justPaintText = false;
		doSub();
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
			box.draw(g, this);
			if (latt != null) drawLittleCirclesAndLines(g);
		}
//		g.drawString(scale.toString(), writepoint[0], writepoint[1]);
		drawText(g);
		
//		g.setColor(java.awt.Color.BLUE);
	}
	public void drawText(Graphics g)
	{
		g.setColor(Color.BLACK);
		g.clearRect(writepoint[0], writepoint[1]-linesize , 500, sy);
		String it = "";
		it += "Current mouse position in map: (" + currentift + ", " + currentjft + ")" + "\r\n";
		it += "More precisely, (" + currentxft + ", " + currentyft + ") \r\n";
		if (latt != null){
			it += "The position is " + Printer.vectorP(currentAtomft) + " in lattice coordinates.\r\n";
			it += "Pixel coordinates (inverse) gives " + Printer.vectorP(currentPixft) + " in pixels.\r\n";
			it += atomDesc;
//			it += latt.forDisplay2();
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
					if (drawCircles) drawCircle(g, cx, cy, 4, Color.BLUE);
					//now draw the lines. Only draw once.
					if (drawLines && ((i%2==0 && j%2==0)||(i%2==1 && j%2==1)) && i >= 1 && i+1 < latticeSites.length && j >= 1 && j+1 < latticeSites[0].length)
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
		if (arg0.getKeyChar() == 'b' || arg0.getKeyChar() == 'B' || arg0.getKeyCode() == KeyEvent.VK_PLUS)
			resizeFT(boxsize*2);
		if (arg0.getKeyChar() == 'r' || arg0.getKeyChar() == 'R' || arg0.getKeyCode() == KeyEvent.VK_MINUS)
			resizeFT(boxsize/2);
		if (arg0.getKeyChar() == 'n' || arg0.getKeyChar() == 'N')
			resizeFT(Integer.parseInt(JOptionPane.showInputDialog(null, "Enter the desired box size.")));
		if (arg0.getKeyChar() == 'a' || arg0.getKeyChar() == 'A')
			moveBox(-1, 0);
		if (arg0.getKeyChar() == 'd' || arg0.getKeyChar() == 'D')
			moveBox(+1, 0);
		if (arg0.getKeyChar() == 'w' || arg0.getKeyChar() == 'W')
			moveBox(0, -1);
		if (arg0.getKeyChar() == 's' || arg0.getKeyChar() == 's')
			moveBox(0, +1);
		if ((arg0.getKeyChar() == 'c' || arg0.getKeyChar() == 'C') && latt != null)
		{
			drawCircles = !drawCircles;
			justPaintText = false;
			repaint();
		}
		if ((arg0.getKeyChar() == 'l' || arg0.getKeyChar() == 'L') && latt != null)
		{
			drawLines = !drawLines;
			justPaintText = false;
			repaint();
		}
		if(arg0.getKeyChar() == 'u' || arg0.getKeyChar() == 'U')
		{
			hiRes = !hiRes;
			doSub();
			repaint();
		}
		if (arg0.getKeyChar() == ' ')
		{
			JFileChooser fc = new JFileChooser(dir);
			if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION){
				SRAW.writeImage(fc.getSelectedFile().toString(), ftmag, ftscale);
				SRAW.writeImage(fc.getSelectedFile().toString()+"rescaled", ftmag);
				SRAW.writeImage(fc.getSelectedFile().toString()+"highres", imageft);
				ColumnIO.writeBin(ftmag, fc.getSelectedFile().toString() + ".dat");
				JOptionPane.showMessageDialog(null, "Done");
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
		if (latt != null) {
			latt.putAtomicCoords(currentxft, currentyft, currentAtomft);
			latt.putPixelCoords(currentAtomft, currentPixft);
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
	 private static class Box
	 {
		 int x, y, size;

		public Box(int x, int y, int size) {
			this.x = x;
			this.y = y;
			this.size = size;
		}
		
		public void draw(Graphics g, AreaSelector parent)
		{
			g.setColor(Color.BLUE);
			g.drawRect(parent.ox + x/parent.zoomFactor, parent.oy + y/parent.zoomFactor, size/parent.zoomFactor, size/parent.zoomFactor);
		}
		
		public void move(int x, int y, AreaSelector parent)
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
			AreaSelector parent;
			JFrame frame;
			public JSlider min, max;
			public BoundsPanel(AreaSelector parent, JFrame frame)
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

	 public static void main(String[] args)
	 {
//		 new AreaSelector(FileOps.openBin(""), Topomap.stddir);
		 new AreaSelector(Layer.open(), Topomap.stddir);
	 }

}

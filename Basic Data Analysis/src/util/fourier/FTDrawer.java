package util.fourier;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.*;
//This class only draws real scalar fields.

//Superseded: This class will soon be able to draw two-component fields in principle.
//It may draw a two-component field consisting of amplitude and phase data, under the assumption that
//the phase is periodic with period 2pi.

public class FTDrawer extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{

	double[][] masterField; int N;
	double field[][]; double fmax, fmin, fdelta;
	boolean twoComp = false;
	BufferedImage image, imageft;
	
	ColorScale scale, ftscale;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	int ftox = ox + sx + 10, ftoy = oy;
	
	int ftsizemin = 1024;
	int ftsizeratio = 1;
	
	int[] writepoint = {ox, oy + sy + 10};
	
	int zoomLevel = 0, zoomFactor = 1;
	
	int WIDTH = 1920, HEIGHT = 1060;
	
	Box box;
	FFT2D_Subset fft;
	double[][] ftmag;
	boolean ftlog = true;
	int[] ftwritepoint = {ftox, ftoy};
	int ftsize;
	
	double ftmax, ftmin;
	double upfrac = 1, downfrac = 0;
	double ftscalemin, ftscalemax;
	
	int currentx, currenty;
	//The field must be 2^n by 2^n;
	public FTDrawer(double[][] field, boolean log)
	{
		ftlog = log;
		masterField = field;
		N = masterField.length;
		ftsize = N;
		sx = N;
		this.field = field;
		while (sx > 1024)
		{
			this.field = FieldOps.reduce(2, this.field);
			sx = this.field.length;
			zoomLevel++;
			zoomFactor *= 2;
		}
		sy = this.field[0].length;
		ftox = ox + sx + 10; ftoy = oy;
		writepoint = new int[] {ox, oy + sy + 10};
		setFieldInfo();
		image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		imageft = new BufferedImage(N, N, BufferedImage.TYPE_INT_RGB);
		box = new Box(0, 0, field.length);
		System.out.println("image");
		try {
			Thread.sleep(5000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		fft = new FFT2D_Subset(masterField, N);
		ftmag = new double[N][N];
		doFFT();
		resetColorScale(0, 1);
		resetFTColorScale(0, 1);
		formFTImage();
//		resetFRColorScale(0, 1);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		showWindow();
		setTitle("FFT region: " + ftsize + "x" + ftsize);

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
	public void doFFT()
	{
		fft.doFFT(box.x, box.y);
		fft.setFHat2();
		FieldOps.magnitude(fft.fHat2, ftmag);
		if (ftlog) FieldOps.log(ftmag);
		ftmax = ArrayOps.max(ftmag);
		ftmin = ArrayOps.min(ftmag);
		resetFTColorScale();
		formFTImage();
	}
	
	
	public void resizeFT(int size)
	{
		fft = new FFT2D_Subset(masterField, size);
		ftsize = size;
		box.setSize(size);
		if (size <= ftsizemin)
			ftsizeratio = ftsizemin/size;
		imageft = new BufferedImage(Math.max(size, ftsizemin), Math.max(size, ftsizemin), BufferedImage.TYPE_INT_RGB);
		ftmag = new double[size][size];
		setTitle("FFT region: " + size + "x" + size);
		doFFT();
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
	public void redrawColorScale(double downnum, double upnum, double downnump, double upnump)
	{
		resetFTColorScale(downnum, upnum);
		resetColorScale(downnump, upnump);
		formFTImage();
		formImage();
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
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				image.setRGB(i, j, scale.of(field[i][j]).getRGB());
	}
	public void formFTImage()
	{
		for (int i = 0; i < imageft.getWidth(); i++)
			for (int j = 0; j < imageft.getHeight(); j++)
				imageft.setRGB(i, j, ftscale.of(ftmag[i/ftsizeratio][j/ftsizeratio]).getRGB());
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
		g.clearRect(0, 0, 2000, 2000);
		drawImage(g);
		drawFTImage(g);
		g.drawString(scale.toString(), writepoint[0], writepoint[1]);
		g.drawString(ftscale.toString(), ftwritepoint[0], ftwritepoint[1]);
		
		g.setColor(java.awt.Color.BLUE);
		box.draw(g, this);
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
			resizeFT(ftsize*2);
		if (arg0.getKeyChar() == 's' || arg0.getKeyChar() == 'S' || arg0.getKeyCode() == KeyEvent.VK_MINUS)
			resizeFT(ftsize/2);
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent arg0) {
			
	}

	@Override
	public void mouseDragged(MouseEvent arg0) {
		// TODO Auto-generated method stub
		int x = arg0.getX();
		int y = arg0.getY();
		box.move(box.x + (x - currentx)*zoomFactor, box.y + (y - currenty)*zoomFactor, this);
		doFFT();
		repaint();
//		System.out.println(x);
		currentx = x; currenty = y;
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
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
		
		public void draw(Graphics g, FTDrawer parent)
		{
			g.drawRect(parent.ox + x/parent.zoomFactor, parent.oy + y/parent.zoomFactor, size/parent.zoomFactor, size/parent.zoomFactor);
		}
		
		public void move(int x, int y, FTDrawer parent)
		{
			this.x = x > 0 ? x+size < parent.N ? x : parent.N - size : 0; this.y = y > 0 ? y+size < parent.N ? y : parent.N - size : 0;
		}
		public void setSize(int size)
		{
			this.size = size;
		}
	 }
	 public static class BoundsPanel extends JPanel implements ChangeListener
		{
			FTDrawer parent;
			JFrame frame;
			public JSlider min, max;
			public JSlider minp, maxp;
			public BoundsPanel(FTDrawer parent, JFrame frame)
			{
				this.frame = frame;
				this.parent = parent;
				this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
				this.setBorder(new LineBorder(Color.GRAY));
				min = new JSlider(0, 1000, 0);
				min.addChangeListener(this);
				max = new JSlider(0, 1000, 1000);
				max.addChangeListener(this);
				minp = new JSlider(0, 1000, 0);
				minp.addChangeListener(this);
				maxp = new JSlider(0, 1000, 1000);
				maxp.addChangeListener(this);

				add(min);
				add(max);
				add(minp);
				add(maxp);
			}
			@Override
			public void stateChanged(ChangeEvent arg0) {
				parent.redrawColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000), ((double)minp.getValue()/1000), ((double)maxp.getValue()/1000));
				frame.setTitle("Min: " + parent.ftscalemin + "   Max: " + parent.ftscalemax + ";     " + parent.scale.getMin() + "  Max: " + parent.scale.getMax());
			}
		}

	 public static void main(String[] args)
	 {
	 	String dir = "C:\\data\\impuritycount\\test8302011\\topo_r\\fview\\";
		JFileChooser fc = new JFileChooser(dir);
		File f = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
			f = fc.getSelectedFile();
		new FTDrawer(ColumnIO.readSquareTable(f.toString()), true);
	 }

}

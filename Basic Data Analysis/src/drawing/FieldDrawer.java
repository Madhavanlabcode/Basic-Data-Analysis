package drawing;

import java.awt.Graphics;
import java.awt.event.*;
import java.awt.image.BufferedImage;

import javax.swing.JFrame;

import util.color.ColorScale;
import util.color.ColorScale2d;
import util.color.ColorScales;
import util.*;
//This class only draws real scalar fields.

//Superseded: This class will soon be able to draw two-component fields in principle.
//It may draw a two-component field consisting of amplitude and phase data, under the assumption that
//the phase is periodic with period 2pi.

public class FieldDrawer extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	
	private static final long serialVersionUID = 6892082910121387055L;
	double field[][]; double fmax, fmin, fdelta;
	double phase[][] = null;
	boolean twoComp = false;
	
	BufferedImage image;
	
	ColorScale scale;
	ColorScale2d scale2;
	
	GradCalculator gradcalc;
	CurlCalculator curlcalc = null;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	
	int WIDTH = 800, HEIGHT = 600;
	
	public DataZoomPanel panel = null;
	
	public FieldDrawer(double[][] field)
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
	public FieldDrawer(double[][] field, double[][] phase)
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
		repaint();
	}
	public void resetColorScale(double downnum, double upnum)
	{
		if(!twoComp)
			scale = new ColorScales.LinearBRYW(fmin + fdelta*upnum, fmin + fdelta*downnum);
		else
			scale2 = new ColorScales.MYC2d(fmin + fdelta*upnum, fmin + fdelta*downnum, 2*Math.PI);
	}
	public void formImage()
	{
		if(!twoComp)
			for (int i = 0; i < sx; i++)
				for (int j = 0; j < sy; j++)
					image.setRGB(i, j, scale.of(field[i][j]).getRGB());
		else
			for (int i = 0; i < sx; i++)
				for (int j = 0; j < sy; j++)
					image.setRGB(i, j, scale2.of(field[i][j], phase[i][j]).getRGB());
			
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
	
	public void paint(Graphics g)
	{
		drawImage(g);
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
//		if (arg0.getKeyChar() == 's' || arg0.getKeyChar() == 'S' && arg0.getModifiers() == KeyEvent.CTRL_DOWN_MASK)
//		{
//			JFileChooser chooser = new JFileChooser();
//			chooser.showSaveDialog(this);
//		}
	}

	@Override
	public void mouseWheelMoved(MouseWheelEvent arg0) {
		// TODO Auto-generated method stub
		int x = arg0.getX();
		int y = arg0.getY();
		int value = arg0.getWheelRotation();
		int[] r = getQuad(x, y);
		System.out.println(r[0] + ", " + r[1] + ", " + value);
		if (value < 0 && panel != null){
			panel.zoomInOnce(r[0], r[1]);
			changeField(panel.field, false);
		}
		
			
	}

	@Override
	public void mouseDragged(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
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
		public FieldDrawer parent;
		public double [][][] grad;
		public double [][] field;
		public double [][] dfds;
		public double theta;
		
		public GradCalculator(FieldDrawer parent, double theta)
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
		public FieldDrawer parent;
		public double [][][] cart;
		public double [][] curl;
		public double theta;
		
		public CurlCalculator(FieldDrawer parent)
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

		
}

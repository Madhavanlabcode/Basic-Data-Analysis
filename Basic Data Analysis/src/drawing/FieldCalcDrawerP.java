package drawing;

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

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import main.SRAW;

import util.FieldOps;
import util.calc.Calculation2D_1P;
import util.calc.GaussSmoothFT;
import util.color.ColorScale;
import util.color.ColorScale2d;
import util.color.ColorScales;

public class FieldCalcDrawerP extends JFrame implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	
	private static final long serialVersionUID = 7408964694534051181L;
	Calculation2D_1P calc;
	int N;
	double fmax, fmin, fdelta;
	double[][] drawField;
	double[][][] drawFieldC;
	boolean real = false;
	
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
	
	double para = 1;
	
	int currentx, currenty;
	int calcx = 0, calcy = 0;
	String dir;
	double[] pbounds;
	SliderPanel s;
	int snpts;
	public FieldCalcDrawerP(Calculation2D_1P calc, String dir, double[] pbounds, boolean realfirst, int slidernpts)
	{
		snpts = slidernpts;
		this.real = realfirst;
		this.calc = calc;
		this.dir = dir;
		this.pbounds = pbounds;
		if (real) N = calc.realField().length;
		else N = calc.complexField().length;
		calc.redoCalc(pbounds[1]);
		sx = N;
		this.drawField = calc.realField();
		this.drawFieldC = calc.complexField();
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
		resetColorScale(0, 1);
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
	public void doCalc()
	{
		calc.redoCalc(para);
	}
	public void setPara(int slider, int npts)
	{
		para = (pbounds[0] + slider*(pbounds[1]-pbounds[0]))/(npts);
	}
	public void setFieldInfo()
	{
		double[] bounds = calc.getMinMax(real);
		fmax = bounds[1]; fmin = bounds[0];
		fdelta = fmax - fmin;
		setTitle("Range: [" + fmin + ", " + fmax + "]" + "     " + "p = " + para);
	}
	public void resetColorScale(double downnum, double upnum)
	{
		scalemax = fmin + fdelta*upnum;
		scalemin = fmin + fdelta*downnum;
		
		if (real) scale = new ColorScales.LinearBRYW(scalemax, scalemin);
		else cscale = new ColorScales.MYC2d(scalemax, scalemin, 2*Math.PI);
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
		setVisible(true);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}	
	public void paint(Graphics g)
	{
		if (refresh){
			g.clearRect(0, 0, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			refresh = false;
		}
		drawText(g);
	}
	public void drawText(Graphics g)
	{
		g.clearRect(writepoint[0], writepoint[1]-linesize, sx, sy);
		String it = calc.getInfo(calcx, calcy);
		String[] lines = it.split("\r\n");
		for (int i = 0; i < lines.length; i++)
			g.drawString(lines[i], writepoint[0], writepoint[1] + i*linesize);
	}
	
	public void resetGraphics()
	{
		if (real) drawField = calc.realField();
		else drawFieldC = calc.complexField();
		if (real)
			while (drawField.length > 1024)
			{
				this.drawField = FieldOps.reduce(2, this.drawField);
				sx = this.drawField.length;
				zoomLevel++;
				zoomFactor *= 2;
			}
		else
			while (drawFieldC.length > 1024)
			{
				this.drawFieldC = FieldOps.reduce(2, this.drawFieldC);
				sx = this.drawFieldC.length;
				zoomLevel++;
				zoomFactor *= 2;
			}
		setFieldInfo();
//	    resetColorScale();
	    this.formImage();
	    refresh = true;
	    repaint();
	}
	public void getNewPBounds()
	{
		pbounds[0] = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the new minimum value."));
		pbounds[1] = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the new maximum value."));
	}
	public void getNewPara()
	{
		para = Double.parseDouble(JOptionPane.showInputDialog(null, "Enter the parameter value."));
		doCalc();
		resetGraphics();
		s.frame.setTitle("p = " + para);
	}
	//if 
	public void zoom(double factorX, double factorY)
	{
		double[] b = calc.getFieldBounds();
		if (b == null) return;
		double deltax = b[1] - b[0];
		double deltay = b[3] - b[2];
		double x = (((double)calcx)/sx)*deltax + b[0];
		double y = (((double)calcy)/sy)*deltay + b[2];
		double xminnew = x + (b[0]-x)/factorX;
		double xmaxnew = x + (b[1]-x)/factorX;
		double yminnew = y + (b[2]-y)/factorY;
		double ymaxnew = y + (b[3]-y)/factorY;
		calc.resize(xminnew, xmaxnew, yminnew, ymaxnew);
	}
	public void translate(int dx, int dy)
	{
		double[] b = calc.getFieldBounds();
		if (b == null) return;
		double deltax = b[1] - b[0];
		double deltay = b[3] - b[2];
		double xmin = (((double)dx)/sx)*deltax + b[0];
		double ymin = (((double)dy)/sy)*deltay + b[2];
		calc.resize(xmin, xmin + deltax, ymin, ymin+deltay);
	}
	
	public void keyPressed(KeyEvent arg0) {
	}
	public void keyReleased(KeyEvent arg0) {
	}
	public void keyTyped(KeyEvent arg0) {
		System.out.println(arg0.getKeyChar());
		if (arg0.getKeyChar() == ' '){
			calc.switchDisplayField();
			resetGraphics();
		}
		if (arg0.getKeyChar() == 's' || arg0.getKeyChar() == 'S')// && arg0.isControlDown())
			calc.save(dir);
		if (arg0.getKeyChar() == 'r' || arg0.getKeyChar() == 'R')// && arg0.isControlDown())
			real = !real;
		if (arg0.getKeyChar() == 'b' || arg0.getKeyChar() == 'B')// && arg0.isControlDown())
			getNewPBounds();
//		if (arg0.getKeyChar() == ' ')
//			{
//				real = !real;
//				resetGraphics();
//			}
		if (arg0.getKeyChar() == 'p' || arg0.getKeyChar() == 'P')// && arg0.isControlDown())
			getNewPara();
	}
	public void mouseWheelMoved(MouseWheelEvent arg0) {
		int notches = arg0.getWheelRotation();
	    boolean negative = notches < 0;
	    notches = Math.abs(notches);
	    double in = 1.1, out = 0.9;
	    if (true)
	    	negative = !negative;
	    if (arg0.isAltDown())
	    	for (int i = 0; i < notches; i++)
	    		zoom(1, negative ? out : in);
	    else if (arg0.isControlDown())
	    	for (int i = 0; i < notches; i++)
	    		zoom(negative ? out : in, 1);
	    else
	    	for (int i = 0; i < notches; i++)
	    		zoom(negative ? out : in, negative ? out : in);
	    calc.redoCalc(para);
	    resetGraphics();
	}
	public void mouseDragged(MouseEvent arg0) {
		int x = arg0.getX();
		int y = arg0.getY();
		translate(currentx - x, currenty - y);
//		System.out.println(x);
		currentx = x; currenty = y;
		calc.redoCalc(para);
		resetGraphics();
	}
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		calcx = currentx - ox;
		calcy = currenty - oy;
		calcx = Math.max(calcx, 0); calcx = Math.min(calcx, N-1);
		calcy = Math.max(calcy, 0); calcy = Math.min(calcy, N-1);
		refresh = false;
		repaint();
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
	
	public static class SliderPanel extends JPanel implements ChangeListener
	{
		private static final long serialVersionUID = 3900716018670779059L;
		FieldCalcDrawerP parent;
		JFrame frame;
		public JSlider s;
		public JSlider min, max;
		int parentmember;
		int oldvalue = 0;
		int oldminv = 0, oldmaxv = 999;
		static int npts = 1001;
		
		public SliderPanel(FieldCalcDrawerP parent, JFrame frame)
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
			frame.setVisible(true);
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
//			if (s.getValue() == oldvalue && min.getValue() == oldminv && max.getValue() == oldmaxv) return;
			if (min.getValue() != oldminv || max.getValue() != oldmaxv)
			{
				parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
				parent.formImage();
				parent.refresh = true;
				parent.repaint();
				oldminv = min.getValue();
				oldmaxv = max.getValue();
				return;
			}
//			else if (s.getValue() == oldvalue) return;
//			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
			parent.setPara(s.getValue(), npts);
			parent.doCalc();
			parent.resetGraphics();
			frame.setTitle(parent.calc.getInfo(parent.para));
			System.out.println(s.getValue());
			oldvalue = s.getValue();
		}
	}

	public static void main(String[] args)
	{
		String dir = "C:\\data\\lawlerhex\\run153topo15_830\\";
//		dir = "C:\\data\\lintrans\\run136topo6_620\\cutoff_pt4a\\";
//		dir = "C:\\data\\lawlerhex\\superm\\";
//		dir = "C:\\data\\lintrans\\run136topo6_620\\cutoff_pt4a\\";
		dir = "C:\\data\\analysis\\";
//		dir = "C:\\data\\lawlerhex\\flat8022010\\";
		//		Calculation2D c = new util.calc.UCalcStorer(dir, "topo", "topocorrect", "topoxL16cont", "topoyL16cont");
		//int i = 6;
		Calculation2D_1P c = null;
//		Calculation2D c = new util.calc.LinearTransBragg2(dir, "didv1", Math.PI/3, "bragg5", 1000, 4.3717*Math.sqrt(3)/2, "");
//		Calculation2D c = Translation.makeTranslation(dir + "topo.dat", 0, 1.1);
//		c = new util.calc.UFieldCalcFT(dir, "topo", 0);
//		c = GaussSmoothFT.getNew(dir, "8-27-2011--5k-0023", 10);
//		c = GaussSmoothFT.getNew(dir, "0023L85locavg");
//		new FieldCalcDrawerP(c, dir, new double [] {0, 40}, true, 200);
		
//		dir = "C:\\data\\lintrans\\run145tm_723\\";
//		c = TopomapView.getNew(dir + "topomap.txt", dir + "topo\\");
		c = GaussSmoothFT.getNew(dir);
//		int n = ((TopomapView)c).nlayers;
//		new FieldCalcDrawerP(c, dir, new double [] {0, n-1}, true, (n-1));
		new FieldCalcDrawerP(c, dir, new double [] {0, 100}, true, 100);
	}

}

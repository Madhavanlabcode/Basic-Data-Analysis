package drawing;

import image.ImageEditing;

import java.awt.*;

import javax.swing.*;

import main.SRAW;

import java.awt.event.*;

import util.color.ColorScale1D;
import util.fileops.FileOps;
import util.fileops.Topomap;
import util.*;
import util.scalar.Function;

//import graph.*;
import java.awt.image.BufferedImage;

public class GraphDrawerCart extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener
{
	private static final long serialVersionUID = -8088494909383992486L;
	static int DEFAULT_WIDTH = 1000, DEFAULT_HEIGHT = 1040;
	int WIDTH = 1000, HEIGHT = 1040;
	static int x0 = 5, y0 = 30;
	static int topSpace = 50, bottomSpace = 50, rightSpace = 40, leftSpace = 80;
	
	double aspect = (double)(WIDTH - rightSpace - leftSpace)/(HEIGHT - topSpace - bottomSpace);
	
	Point mousePos;
	Point topCorner = new Point(x0 + rightSpace, y0 + topSpace);
	
	static Color[] drawColors = {Color.RED, Color.BLUE, Color.MAGENTA, Color.CYAN, Color.ORANGE, Color.YELLOW, Color.GREEN};

	GraphObject[] plot = new GraphObject[7];
	
	BufferedImage image = null;
	ColorScale1D imageScale;
	double[][] yData;
	
	boolean usePlot3 = false;
	
	double minX = 0, maxX = WIDTH, minY = 0, maxY = HEIGHT;
	double xScale, yScale;
	boolean paintingArrays = true;
	boolean fixedDomain = true, fixedRange = true;
	
	boolean plot0Changed = false;
	
	public LineCutDrawer parent = null;
	public GraphDrawer_User user = null;
//	boolean drawTickMarks = true;
	
//	boolean drawMultipleYs = false;
	
	
	boolean drawSpecialMessage = false;
	String specialMessage;
	java.awt.Point specialMessagePoint = new java.awt.Point(WIDTH - (rightSpace + 250), topSpace);
	
	boolean drawAxisLabels = false;
	String xAxisLabel = "X-Axis";
	Point xAxisLabelPoint;
	String yAxisLabel = "Y-Axis";
	Point yAxisLabelPoint;

	String y2AxisLabel = "Y-Axis";
	Point y2AxisLabelPoint;

	String y3AxisLabel = "Y-Axis";
	Point y3AxisLabelPoint;
	
	boolean switchZooming = false;
	
	public Image dbimage;
	Graphics dbg;
	
	public int[] currentMouseScreen = new int [2];
	public double[] currentMouseData = new double [2];
	
	public KeyStrokeProccessor extra = null;
	
	public JFileChooser fc;
	public boolean drawingLayer = false;
	public util.fileops.Layer layer = null;
	public util.color.ColorScale1D cscale = null;
	public BufferedImage layerImage = null;
	
    public GraphDrawerCart(String title, boolean yesGraphNoImage)
    {
    	if (!yesGraphNoImage)
    		paintingArrays = false;
    	
    	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	setSize(WIDTH, HEIGHT);
    	setTitle(title);
    	this.setFont(new Font("monospaced", Font.ROMAN_BASELINE, 12));
//    	yLabel.setDrawingComponent(this);
//    	yLabel.setBackground(Color.WHITE);
    	this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	//////////////////////////////////////////
    	this.addComponentListener(new ComponentListener()    {     
        public void componentResized(ComponentEvent evt) {
            GraphDrawerCart  g = (GraphDrawerCart)evt.getSource();
            // Get new size
            Dimension d = g.getSize();
           	WIDTH = d.width;
        	HEIGHT = d.height;
        	aspect = (double)(WIDTH - rightSpace - leftSpace)/(HEIGHT - topSpace - bottomSpace);
        	g.setUpScale();
        	g.repaint();
    		}

		public void componentHidden(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentMoved(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentShown(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}
    		}
    	);
    ////////////////////////////////////////	
    	addMouseListener(this);
    	addMouseMotionListener(this);
    	addMouseWheelListener(this);
    	addKeyListener(this);
    }
    public GraphDrawerCart(String title, double[] x, double[] data)
    {
    	paintingArrays = true;
    	
    	plot[0] = new GraphObject(x, data);
    	
    	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	setSize(WIDTH, HEIGHT);
    	setTitle(title);
    	this.setFont(new Font("monospaced", Font.ROMAN_BASELINE, 12));
//    	yLabel.setDrawingComponent(this);
//    	yLabel.setBackground(Color.WHITE);
    	this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	this.addComponentListener(new ComponentListener()    {     
        public void componentResized(ComponentEvent evt) {
            GraphDrawerCart  g = (GraphDrawerCart)evt.getSource();
    
            // Get new size
            Dimension d = g.getSize();
           	WIDTH = d.width;
        	HEIGHT = d.height;
        	aspect = (double)(WIDTH - rightSpace - leftSpace)/(HEIGHT - topSpace - bottomSpace);
        	g.setUpScale();
        	g.repaint();
    		}
		public void componentHidden(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentMoved(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentShown(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}
    		}
    	);
    	
    	setDomain(x[0], x[x.length-1]);
    	setRange(ArrayOps.min(data), ArrayOps.max(data));
    	
    	addMouseListener(this);
    	addMouseMotionListener(this);
    	addMouseWheelListener(this);
    	addKeyListener(this);
    }
    public GraphDrawerCart(String title, double[] x, double[] data, double xmin, double xmax, double ymin, double ymax)
    {
    	paintingArrays = true;
    	
    	plot[0] = new GraphObject(x, data);
    	
    	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	setSize(WIDTH, HEIGHT);
    	setTitle(title);
    	this.setFont(new Font("monospaced", Font.ROMAN_BASELINE, 12));
//    	yLabel.setDrawingComponent(this);
//    	yLabel.setBackground(Color.WHITE);
    	this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	this.addComponentListener(new ComponentListener()    {     
        public void componentResized(ComponentEvent evt) {
            GraphDrawerCart  g = (GraphDrawerCart)evt.getSource();
    
            // Get new size
            Dimension d = g.getSize();
           	WIDTH = d.width;
        	HEIGHT = d.height;
        	aspect = (double)(WIDTH - rightSpace - leftSpace)/(HEIGHT - topSpace - bottomSpace);
        	g.setUpScale();
        	g.repaint();
    		}
		public void componentHidden(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentMoved(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentShown(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}
    		}
    	);
    	
    	setDomain(xmin, xmax);
    	setRange(ymin, ymax);
    	
    	addMouseListener(this);
    	addMouseMotionListener(this);
    	addMouseWheelListener(this);
    	addKeyListener(this);
    }
    public GraphDrawerCart(String title, double[][][] data)//data is [n][2][length] where data[i][0] is x and data[i][1] is y.
    {
    	paintingArrays = true;
    	double minx = Double.MAX_VALUE, miny = Double.MAX_VALUE;
    	double maxx = -Double.MAX_VALUE, maxy = -Double.MAX_VALUE;
    	plot = new GraphObject[data.length];
    	for (int i = 0; i < data.length; i++)
    	{
    		minx = Math.min(minx, ArrayOps.min(data[i][0]));
    		maxx = Math.max(maxx, ArrayOps.max(data[i][0]));
    		miny = Math.min(miny, ArrayOps.min(data[i][1]));
    		maxy = Math.max(maxy, ArrayOps.max(data[i][1]));
        	plot[i] = new GraphObject(data[i][0], data[i][1]);
        	plot[i].c = drawColors[i % drawColors.length];
    		
    	}
    	
    	setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	setSize(WIDTH, HEIGHT);
    	setTitle(title);
    	this.setFont(new Font("monospaced", Font.ROMAN_BASELINE, 12));
//    	yLabel.setDrawingComponent(this);
//    	yLabel.setBackground(Color.WHITE);
    	this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    	this.addComponentListener(new ComponentListener()    {     
        public void componentResized(ComponentEvent evt) {
            GraphDrawerCart  g = (GraphDrawerCart)evt.getSource();
    
            // Get new size
            Dimension d = g.getSize();
           	WIDTH = d.width;
        	HEIGHT = d.height;
        	aspect = (double)(WIDTH - rightSpace - leftSpace)/(HEIGHT - topSpace - bottomSpace);
        	g.setUpScale();
        	if (!paintingArrays){ resizeImage();
        	refreshImage();
        	}
        	g.repaint();
    		}
		public void componentHidden(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentMoved(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		public void componentShown(ComponentEvent arg0) {
			// TODO Auto-generated method stub
			
		}
    		}
    	);
    	
    	setDomain(minx, maxx);
    	setRange(miny, maxy);
    	
    	addMouseListener(this);
    	addMouseMotionListener(this);
    	addMouseWheelListener(this);
    	addKeyListener(this);
    }

    public static GraphDrawerCart getNew(double[] x, double[] y)
    {
    	return new GraphDrawerCart("Title", new double[][][] {{x, y}});
    }
    public void actionPerformed(ActionEvent ae)
    {
//    	Dimension d = getSize();
    }
   public void setAxisLabels(String x, String y)
   {
	   xAxisLabel = x;
	   yAxisLabel = y;
	   drawAxisLabels = true;
   }
   public void setY2Label(String y2)
   {
	   y2AxisLabel = y2;
   }
   public void setY3Label(String y3)
   {
	   y3AxisLabel = y3;
   }   
    public void paint(Graphics g)
    {
    	dbg.setColor(new Color(255, 255, 255));
    	dbg.fillRect(0, 0, WIDTH, HEIGHT);
    	
    	dbg.setClip(0, 0, WIDTH, HEIGHT);
    	drawCartPlane(dbg);
		if (drawingLayer)
			drawLayer(layer, cscale, dbg);
		drawAxesArrays(dbg);
//		int xpi;
		if (paintingArrays)
    		for (int i = 0; i < plot.length; i++)
    		{
    			if (plot[i] != null)
    				plot[i].draw(dbg, this);
    		}	
    	else{
    		if (image == null) resizeImage();
    		refreshImage();
//    		if (plot[0] != null) {xpi = pixel(plot[0].x[0], 0).x;
    		dbg.drawImage(image, leftSpace, topSpace, this);
//    		}
    	}
//    	if (drawAxisLabels)
//    		drawAxisLabels(g);
    	if (drawSpecialMessage)
    		drawSpecialMessage(dbg);
		
    	g.drawImage(dbimage, 0, 0, this);
    }
    
    /**
     * The image is established as follows. 
     */
    public void refreshImage()
    {
    	ImageEditing.putBufferedImage_FixedSizeIntTranspose(yData, imageScale, image);
    }
    
    public void resizeImage()
    {
    	image = new BufferedImage(WIDTH - rightSpace - leftSpace, HEIGHT - topSpace - bottomSpace, BufferedImage.TYPE_INT_RGB);
    }
    
    public void paintRaw(Graphics g)
    {
    	g.setColor(new Color(255, 255, 255));
    	g.fillRect(0, 0, WIDTH, HEIGHT);
    	
    	g.setClip(0, 0, WIDTH, HEIGHT);
    	drawCartPlane(g);
		if (drawingLayer)
			drawLayer(layer, cscale, g);
    	drawAxesArrays(g);

		if (paintingArrays)
		{
    		for (int i = 0; i < plot.length; i++)
    			if (plot[i] != null)
    			{
    	    		g.setColor(drawColors[i]);
    				plot[i].draw(g, this);
    			}
    	}
    	else
    		g.drawImage(image, topCorner.x, topCorner.y, this);
//    	if (drawAxisLabels)
//    		drawAxisLabels(g);
    	if (drawSpecialMessage)
    		drawSpecialMessage(g);
		
//    	g.drawImage(dbimage, 0, 0, this);
    }
    public void drawLayer(util.fileops.Layer layer, util.color.ColorScale1D scale, Graphics g)
    {
    	java.awt.Point origin = this.pixel(layer.x[0], layer.y[0]);
    	java.awt.Point corner = this.pixel(layer.x[layer.nx-1], layer.y[layer.ny-1]);
    	
    	java.awt.Point min = new java.awt.Point(Math.min(origin.x, corner.x), Math.min(origin.y, corner.y));
    	java.awt.Point max = new java.awt.Point(Math.max(origin.x, corner.x), Math.max(origin.y, corner.y));
    	
    	g.drawImage(layerImage, min.x, min.y, max.x-min.x, max.y-min.y, null);
    }
    public int[] restrictPt(java.awt.Point point)
    {
    	int[] ans = new int [2];
    	if (point.x < 0) ans[0] = 0;
    	if (point.x >= WIDTH) ans[0] = WIDTH;
    	else ans[0] = point.x;
    	if (point.y < 0) ans[1] = 0;
    	if (point.y >= HEIGHT) ans[1] = HEIGHT;
    	else ans[1] = point.y;
    	return ans;
    }
    public void translate(double dx, double dy)
    {
    	minX += dx; maxX += dx;
    	minY += dy; maxY += dy;
    }
    public void zoom(boolean in, double factor)
    {
    	
    	double xSize = maxX - minX,
    		ySize = maxY - minY;
    	
       	maxX += xSize*(factor-1)/(2);
    	minX -= xSize*(factor-1)/(2);
       	maxY += ySize*(factor-1)/2;
    	minY -= ySize*(factor-1)/2;
    	xScale *= factor;
    	yScale *= factor;
    }
    public void zoomX(boolean in, double factor)
    {
    	double xSize = maxX - minX;
       	maxX += xSize*(factor-1)/(2);
    	minX -= xSize*(factor-1)/(2);
    	xScale *= factor;
    }
    public void zoomY(boolean in, double factor)
    {
    	double ySize = maxY - minY;
       	maxY += ySize*(factor-1)/2;
    	minY -= ySize*(factor-1)/2;
    	yScale *= factor;
    }
    public void setUpScale()
    {
    	xScale = (maxX - minX)/(WIDTH - leftSpace - rightSpace);
    	yScale = (maxY - minY)/(HEIGHT - bottomSpace - topSpace);
    	specialMessagePoint = new java.awt.Point(WIDTH - (rightSpace + 350), topSpace);
    	xAxisLabelPoint = new Point(leftSpace/2 + WIDTH/2 - rightSpace, HEIGHT - bottomSpace + 3*topSpace/4);
    	yAxisLabelPoint = new Point(leftSpace/4, HEIGHT/2 + bottomSpace);
    }
    
    public void drawAxisLabels(Graphics g)
    {
    	g.drawString(xAxisLabel, xAxisLabelPoint.x, xAxisLabelPoint.y);
//    	g.drawString(yAxisLabel, yAxisLabelPoint.x, yAxisLabelPoint.y);

    }
    
    public void drawAxesArrays(Graphics g)
    {
    	g.setColor(Color.BLACK);
    	g.drawLine(leftSpace, HEIGHT - bottomSpace, WIDTH - rightSpace, HEIGHT - bottomSpace);
    	g.drawLine(leftSpace, topSpace, leftSpace, HEIGHT - bottomSpace);
    }
    
    public void setSpecialMessage(boolean isDrawing, String message)
    {
    	drawSpecialMessage = isDrawing;
    	this.specialMessage = message;
    }
    public void drawSpecialMessage(Graphics g)
    {
    	String[] lines = specialMessage.split("\n");
    	for (int i = 0; i < lines.length; i++)
    		g.drawString(lines[i], specialMessagePoint.x, specialMessagePoint.y + i*20);
    }
    
    public void drawTickMarks(Graphics g)
    {
    	double xInt = maxX - minX;
    	double log = (int)Math.log10(xInt);
    	double xSpace = Math.pow(10, log > 0 ? (int)log : (int)log - 1);
    	while (xInt/xSpace < 5)
    		xSpace /= 2;
    	
    	double xNow = minX;
    	Point point;
    	while (xNow <= maxX)
    	{
    		point = pixel(xNow, minY);
    		g.drawLine(point.x, point.y, point.x, point.y + 10);
    		g.drawString(String.format("%.2f", xNow), point.x - 10, point.y + 23);
    		xNow += xSpace;
    	}
    	//y
    	double yInt = maxY - minY;
    	log = (int)Math.log10(yInt);
    	double ySpace = Math.pow(10, log > 0 ? (int)log : (int)log - 1);
    	while (yInt/ySpace < 5)
    		ySpace /= 2;
    	
    	double yNow = minY;
    	while (yNow <= maxY)
    	{
    		point = pixel(minX, yNow);
    		g.drawLine(point.x, point.y, point.x - 10, point.y);
    		g.drawString(String.format("%.2f", yNow), point.x - 40, point.y + 5);
    		yNow += ySpace;
    	}
    }
    
    public int[] getPlaneSize()
    {
    	return new int[] {WIDTH - rightSpace - leftSpace - x0, HEIGHT - topSpace - bottomSpace - y0};
    }
    public void drawCartPlane(Graphics g)
    {
    	double xInt = maxX - minX;
    	double log = Math.log10(xInt/2);
    	log = log > 0 ? (int)log : (int)(log - 1);
//    	System.out.println(log);
    	double xSpace = Math.pow(10, log);
    	
    	Point point;
    	//First draw the y-axis:
    	Point top = pixel(0, maxY);
    	Point bottom = pixel(0, minY);
    	
//    	g.setColor(Color.black);
    	
    	double xNow = xSpace*((int)(minX/xSpace));
    	while (xNow <= maxX)
    	{
    		point = pixel(xNow, 0);
    		g.setColor(Color.DARK_GRAY);
    		g.drawLine(point.x, bottom.y, point.x, top.y);
//    		g.drawString(NumFormat.scientific(xNow, (int)Math.max(3, Math.abs(log)+2)), point.x - 10, bottom.y + 23);
    		if (Math.abs(xNow) < xInt/2000) xNow = 0;
    		g.drawString(NumFormat.scientific(xNow, 3), point.x - 10, bottom.y + 23);
//    		g.setColor(Color.LIGHT_GRAY);
//    		for (int i = 1; i < 10; i++)
//    		{
//    			xNow += xSpace/10;
//        		point = pixel(xNow, 0);
//        		g.drawLine(point.x, bottom.y, point.x, top.y);
//    		}
    			
    		xNow += xSpace;
    	}
    	//y
    	//First draw the y-axis:
    	double yInt = maxY - minY;
    	log = Math.log10(yInt/2);
    	log = log > 0 ? (int)log : (int)(log - 1);
//    	System.out.println(log);
    	double ySpace = Math.pow(10, log);
    	
    	Point left = pixel(minX, 0);
    	Point right = pixel(maxX, 0);
    	
    	double yNow = ySpace*((int)(minY/ySpace));
    	while (yNow <= maxY)
    	{
    		point = pixel(0, yNow);
    		g.setColor(Color.DARK_GRAY);
    		g.drawLine(left.x, point.y, right.x, point.y);
//    		g.drawString(NumFormat.scientific(yNow, (int)Math.max(3, Math.abs(log)+2)), left.x - 30, point.y + 2);
    		g.drawString(NumFormat.scientific(yNow, 3), left.x - 30, point.y + 2);
//    		g.setColor(Color.LIGHT_GRAY);
//    		for (int i = 1; i < 10; i++)
//    		{
//    			yNow += ySpace/10;
//        		point = pixel(0, yNow);
//        		g.drawLine(left.x, point.y, right.x, point.y);
//    		}
    			
    		yNow += ySpace;
    	}
    	g.setColor(Color.BLACK);
    	g.drawLine(left.x, left.y, right.x, right.y);
       	g.drawLine(bottom.x, bottom.y, top.x, top.y);
       	
//       	System.out.println(top + ", " + bottom + ", " + left + ", " + right);
    }
    
    public java.awt.Point pixel(double x, double y)
    {
    	return new Point((int)((x - minX)/xScale) + leftSpace, (int)(HEIGHT - bottomSpace - ((y - minY)/yScale)));
    }
    public void putDataPoint(int[] pixels, double[] dataPt)
    {
    	dataPt[0] = (pixels[0]-leftSpace)*xScale + minX;
    	dataPt[1] = ((HEIGHT - bottomSpace)-pixels[1])*yScale + minY;
    }
    
//    public double[] point(int x1, int y1)
//    {
//    	double x, y;
//    	x = (x1 - leftSpace)*xScale + minX;
//    	y = (y1 - leftSpace)*xScale + minY;
//    }
   
    public void setXY(GraphObject curve, boolean newDomain, boolean newRange, int i)
    {
    	plot[i] = curve;
    	
    	if (curve == null) return;
    	if (yData == null) yData = new double [][] {curve.y};

       	if (i >= yData.length)
    	{
    		double[][] temp = yData;
    		yData = new double [yData.length+1][];
    		for (int j = 0; i < temp.length; j++)
    			yData[j] = temp[j];
    	}
//    	yData[i] = curve.y;
    	if (i < 7) curve.c = drawColors[i]; else curve.c = null;
    	if (newDomain)
    	{
    		minX = ArrayOps.min(curve.x);
    		maxX = ArrayOps.max(curve.x);
    	}
    	if (newRange)
    	{
    		minY = ArrayOps.min(curve.y);
    		maxY = ArrayOps.max(curve.y);
    	}
		if (maxX - minX == 0 && maxY - minY == 0)
		{
			maxX += 0.5;
			minX -= 0.5;
			maxY += 0.5/aspect;
			minY -= 0.5/aspect;
		}
		else if (maxY - minY == 0)
		{
			maxY += 0.5*(maxX - minX)/aspect;
			minY -= 0.5*(maxX - minX)/aspect;
		}
		else if (maxX - minX == 0)
		{
			maxX += 0.5*(maxY - minY)*aspect;
			minX -= 0.5*(maxY - minY)*aspect;
		}
//		System.out.println(maxX + ", " + minX + ", " + maxY + ", " + minY + ", " + aspect);
    	if (newDomain || newRange)
    		setUpScale();
    }
    public void setXY2(double[] x, double[] y, boolean useIt)
    {
    	plot[2] = new GraphObject(x, y);
    	usePlot3 = useIt;
    }
    
    //sets the data in plot[0] without changing the range.
    public void setCurve(double[] data)
    {
    	plot[0].y = data;
    	repaint();
    }
    public void setCurve(double[] data, double[] x)
    {
    	plot[0].y = data;
    	plot[0].x = data;
    	repaint();
    }
    public void setCurve(double[] data, double[] x, int index)
    {
    	plot[index].y = data;
    	plot[index].x = x;
    	repaint();
    }
    public void setColor(Color c, int plotIndex)
    {
    	plot[plotIndex].c = c;
    }
    public void setThin(boolean isThin, int plotIndex)
    {
    	plot[plotIndex].drawThick = !isThin;
    }
    public void setNumPlots(int n)
    {
    	this.plot = new GraphObject [n];
    	yData = new double[n][];
    }

    public void setDomain(double minX, double maxX)
    {
    	this.minX = minX;
    	this.maxX = maxX;
    }
    public void setRange(double minY, double maxY)
    {
    	this.minY = minY;
    	this.maxY = maxY;
    	setUpScale();
    }
    public void setRange(double[] y)
    {
    	this.minY = ArrayOps.min(y);
    	this.maxY = ArrayOps.max(y);
    	setUpScale();
    }
    public void resetRange(int plotIndex)
    {
    	setRange(ArrayOps.min(plot[plotIndex].y), ArrayOps.max(plot[plotIndex].y));
    }
    
    public static void main(String[] args)
    {
    	Topomap.setStdDir();
    	
    	JFileChooser fc = new JFileChooser(Topomap.stddir);
    	double[][] data = FileOps.openTable(fc);
//    	plotGraph(data[0], data[1]);
    	plotGraph(data, 2);
    }
    
    public static void plotGraph(double[] x, double[] y)
    {
    	if (x == null)
    	{
    		x = ArrayOps.generateArrayNotInclUpper(0, y.length, y.length);
    	}
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", x, y);
    	drawer.setColor(Color.RED, 0);
		drawer.showWindow();
    }
    public static void plotGraph(double[] x, int[] y)
    {
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", x, ArrayOps.toDouble(y));
    	drawer.setColor(Color.RED, 0);
		drawer.showWindow();
    }
    public static void plotGraphs(double[] x, int[][] y)
    {
    	double[][][] data = new double[y.length][2][];
    	for (int i = 0; i < y.length; i++)
    	{
    		data[i][0] = x;
    		data[i][1] = ArrayOps.toDouble(y[i]);
    	}
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", data);
    	drawer.setColor(Color.RED, 0);
		drawer.showWindow();
    }
    public static void plotGraphs(double[] x, double[][] y)
    {
    	double[][][] data = new double[y.length][2][];
    	for (int i = 0; i < y.length; i++)
    	{
    		data[i][0] = x;
    		data[i][1] = y[i];
    	}
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", data);
    	drawer.setColor(Color.RED, 0);
		drawer.showWindow();
    }
    public static void plotGraphs(double[][] x, double[][] y)
    {
    	double[][][] data = new double[y.length][2][];
    	for (int i = 0; i < y.length; i++)
    	{
    		data[i][0] = x[i];
    		data[i][1] = y[i];
    	}
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", data);
    	drawer.setColor(Color.RED, 0);
		drawer.showWindow();
    }
    public static GraphDrawerCart getNew(double[][] x, double[][] y)
    {
    	double[][][] data = new double[y.length][2][];
    	for (int i = 0; i < y.length; i++)
    	{
    		data[i][0] = x[i];
    		data[i][1] = y[i];
    	}
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", data);
    	return drawer;
    }
   public static void plotGraphs(double[] x, int[][] y, Color[] c)
    {
    	double[][][] data = new double[y.length][2][];
    	for (int i = 0; i < y.length; i++)
    	{
    		data[i][0] = x;
    		data[i][1] = ArrayOps.toDouble(y[i]);
    	}
    	GraphDrawerCart drawer = new GraphDrawerCart("Graph", data);
    	for (int i = 0; i < c.length; i++)
    		drawer.setColor(c[i], i);
		drawer.showWindow();
    }
   public static void plotGraph(double[][] data_x_ys, int modulo)
    {
    	double[] x = data_x_ys[0];
    	GraphDrawerCart drawer = new GraphDrawerCart("", true);
    	if (modulo == 0) modulo = 6;
    	for (int i = 0; i < data_x_ys.length-1; i++)
    	{
    		drawer.setXY(new GraphObject(x, data_x_ys[i+1]), true, true, i);
        	drawer.setColor(drawColors[i%modulo], i);

    	}
		drawer.showWindow();
    }
    
    static void wait(int millis)
    {
    	try {
			Thread.sleep(millis);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
    }
    
    public static class GraphObject
    {
    	int length;
		double[] x, y;
		boolean drawThick = false;
		public Color c;
		boolean isDrawing = true;
		
		public GraphObject(double xmin, double xmax, double[] y) {
			super();
			this.y = y;
			length = y.length;
			
			this.x = new double [length];
			double dx = (xmin - xmax)/(length-1);
			for (int i = 0; i < length; i++)
				x[i] = xmin+i*dx;
		}
		public GraphObject(double[] x, double[] y) {
			super();
			this.x = x;
			this.y = y;
			length = x.length;
		}
		public GraphObject(double xmin, double xmax, Function f, int npts) {
			super();
			double dx = (xmax - xmin)/npts;
			x = new double [npts];
			y = new double [npts];
			for (int i = 0; i < npts; i++){
				x[i] = xmin + i*dx;
				y[i] = f.of(x[i]);
			}
			length = x.length;
		}
		
		public GraphObject(int length){
			x = new double[length];
			y = new double[length];
			this.length = length;
		}

	    public void draw(Graphics g, GraphDrawerCart drawer)
	    {
	    	if (!isDrawing)
	    		return;
	    	
	    	g.setColor(c);
	    	if(drawThick) {drawThick(g, drawer); return;}
	    	
	    	Point first, second;
	    	first = drawer.pixel(x[0], y[0]);
//	    	int start = 0;
//	    	while(first.x < GraphDrawer.leftSpace)
//	    	{
//	    		start++;
//	    		first = drawer.pixel(x[start], y[start]);
//	    	}
	    	
	    	for (int i = 0; i < length - 1; i++)
	    	{
	    		second = drawer.pixel(x[i + 1], y[i + 1]);
	    		if (x[i+1] <= drawer.maxX || x[i+1] >= drawer.minX || y[i+1] >= drawer.minY || y[i+1] <= drawer.maxY)
		    		g.drawLine(first.x, first.y, second.x, second.y);
	    		first = new Point(second.x, second.y);
	    	}
	    }
	    public void drawThick(Graphics g, GraphDrawerCart drawer)
	    {
	    	Point first, second;
	    	first = drawer.pixel(x[0], y[0]);
	    	
//	    	int start = 0;
//	    	while(first.x < GraphDrawer.leftSpace)
//	    	{
//	    		start++;
//	    		first = drawer.pixel(x[start], y[start]);
//	    	}
	    	for (int i = 0; i < length - 1; i++)
	    	{
	    		second = drawer.pixel(x[i + 1], y[i + 1]);
	    		g.drawLine(first.x, first.y, second.x, second.y);
	    		g.drawLine(first.x + 1, first.y, second.x + 1, second.y);
	    		g.drawLine(first.x, first.y + 1, second.x, second.y + 1);
	    		first = new Point(second.x, second.y);
	    	}
	    	for (int i = -1; i < 2; i++)
	    		g.drawLine(first.x + i, first.y - 1, first.x + i, first.y + 1);
	    	first = drawer.pixel(x[0], y[0]);
	    	
	    	for (int i = -1; i < 2; i++)
	    		g.drawLine(first.x + i, first.y - 1, first.x + i, first.y + 1);
	    }

    }
    static class TwoDDist
    {
    	public double x, y;
    	public double dx, dy;
    	
		public TwoDDist(double x, double y, double dx, double dy) {
			this.x = x;
			this.y = y;
			this.dx = dx;
			this.dy = dy;
		}
	    public void draw(Graphics g, GraphDrawerCart drawer)
	    {
	    	Point point = drawer.pixel(x, y), dPoint = new Point((int)(dx/drawer.xScale), (int)(dy/drawer.yScale));
	
	    	g.fillOval(point.x - 1, point.y - 1, 3, 3);
	    	g.drawOval(point.x - dPoint.x, point.y - dPoint.y, dPoint.x*2, dPoint.y*2);
	    }
    }

    static class ThickLine extends GraphObject
    {
		public Color lineColor, endColor;
		public boolean changed = false;
		public ThickLine(double[] x, double[] y, Color lineColor, Color endColor) {
			super(x, y);
			this.lineColor = lineColor;
			this.endColor = endColor;
			
		}
		
		public ThickLine(int length){
			super(length);
		}

	    public void draw(Graphics g, GraphDrawerCart drawer)
	    {
	    	Point first, second;
	    	first = drawer.pixel(x[0], y[0]);
	    	
//	    	int start = 0;
//	    	while(first.x < GraphDrawer.leftSpace)
//	    	{
//	    		start++;
//	    		first = drawer.pixel(x[start], y[start]);
//	    	}
	    	g.setColor(lineColor);
	    	for (int i = 0; i < length - 1; i++)
	    	{
	    		second = drawer.pixel(x[i + 1], y[i + 1]);
	    		g.drawLine(first.x, first.y, second.x, second.y);
	    		g.drawLine(first.x + 1, first.y, second.x + 1, second.y);
	    		g.drawLine(first.x, first.y + 1, second.x, second.y + 1);
	    		first = new Point(second.x, second.y);
	    	}
	    	g.setColor(endColor);
	    	for (int i = -1; i < 2; i++)
	    		g.drawLine(first.x + i, first.y - 1, first.x + i, first.y + 1);
	    	first = drawer.pixel(x[0], y[0]);
	    	
	    	g.setColor(endColor);
	    	for (int i = -1; i < 2; i++)
	    		g.drawLine(first.x + i, first.y - 1, first.x + i, first.y + 1);
	    }

    }

	public void mouseDragged(MouseEvent arg0) {
		// TODO Auto-generated method stub
//		System.out.println("Dragged");
		int dx = arg0.getX() - mousePos.x;
		int dy = arg0.getY() - mousePos.y;
		
		translate(-dx*xScale, dy*yScale);
		plot0Changed = false;
		repaint();
		mousePos.x += dx;
		mousePos.y += dy;
		return;
		
	}
	public void mouseMoved(MouseEvent arg0) {
		currentMouseScreen[0] = arg0.getX();
		currentMouseScreen[1] = arg0.getY();
		putDataPoint(currentMouseScreen, currentMouseData);
		setTitle(Printer.vectorP(currentMouseData));
	}
	public void mouseClicked(MouseEvent arg0) {
		if (user!= null) user.processMouseClick(arg0);
	}
	public void mouseEntered(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	public void mouseExited(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	public void mousePressed(MouseEvent arg0) {
		// TODO Auto-generated method stub
		mousePos = new Point(arg0.getX(), arg0.getY());
	}
	public void mouseReleased(MouseEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	public void mouseWheelMoved(MouseWheelEvent arg0) {
		// TODO Auto-generated method stub

		int notches = arg0.getWheelRotation();
	    boolean negative = notches < 0;
	    if (switchZooming)
	    	negative = !negative;
	    notches = Math.abs(notches);
	    if (arg0.isAltDown())
	    	for (int i = 0; i < notches; i++)
	    		zoomY(negative, negative ? 0.9 : 1.1);
	    else if (arg0.isControlDown())
	    	for (int i = 0; i < notches; i++)
	    		zoomX(negative, negative ? 0.9 : 1.1);
	    else
	    	for (int i = 0; i < notches; i++)
	    		zoom(negative, negative ? 0.9 : 1.1);
	    	
	    repaint();
	}
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub
		if (e.getKeyCode() == KeyEvent.VK_QUOTE)
			switchZooming = !switchZooming;
		
		if (user != null) user.processKeyStroke(e);
		
	}
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub
		
	}
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub
		if (e.getKeyChar() == 's' && extra == null)
		{
			SRAW.writeImage(FileOps.selectSave(fc).toString(), (BufferedImage)dbimage);
		}
		if (extra != null) extra.processKeyStroke(e.getKeyChar());
		else if (e.getKeyChar() == 'o')
		{
	    	double[][] data_x_ys = FileOps.openTable(fc);
	    	double[] x = data_x_ys[0];
	    	for (int i = 0; i < data_x_ys.length-1; i++)
	    	{
	    		setXY(new GraphObject(x, data_x_ys[i+1]), true, true, i);
	        	setColor(drawColors[i%2], i);

	    	}

		}
		if (e.getKeyChar() == 'C')
		{
			ImageEditing.copyToClipboard((BufferedImage)dbimage);
		}
		if (e.getKeyChar() == 'p')
		{
//			for (int i = 0; i < plot[0].length; i++)
//				System.out.println(plot[0].x[i] + "\t" + plot[0].y[i]);
			int np = 0;
			for (int i = 0; i < plot.length; i++)
    		{
    			if (plot[i] != null)
    				np++;
    		}	
			double[][] x = new double [np][];
			double[][] y = new double [np][];
			int n =0;
			for (int i = 0; i < plot.length; i++)
				if (plot[i] != null)
				{
					x[n] = plot[i].x;
					y[n++] = plot[i].y;
				}
			
			Printer.printColumnSeries(x, y);

		}
		if (e.getKeyChar() == 'x')
			System.out.println(currentMouseData[0]);
			
	}
	public void printAll(int zero, int limit)
	{
		int np = 0;
		for (int i = zero; i < limit; i++)
		{
			if (plot[i] != null)
				np++;
		}	
		double[][] x = new double [np][];
		double[][] y = new double [np][];
		int n =0;
		for (int i = zero; i < limit; i++)
			if (plot[i] != null)
			{
				x[n] = plot[i].x;
				y[n++] = plot[i].y;
			}
		
		Printer.printColumnSeries(x, y);
	}
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		setVisible(true);
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		dbimage = createImage(WIDTH, HEIGHT);
		dbg = dbimage.getGraphics();
		try {
			Thread.sleep(1000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		repaint();
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
	}
	
	public static BufferedImage drawPlot(double[] x, double[] y)
	{
		BufferedImage image = new BufferedImage(DEFAULT_WIDTH, DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
		GraphDrawerCart g = new GraphDrawerCart("", x, y);
		g.paintRaw(image.getGraphics());
		return image;
	}
	public static BufferedImage drawPlot(double[] x, double[] y, int width, int height)
	{
		BufferedImage image = new BufferedImage(DEFAULT_WIDTH, DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
		GraphDrawerCart g = new GraphDrawerCart("", x, y);
       	g.WIDTH = width;
    	g.HEIGHT = height;
    	g.aspect = (double)(g.WIDTH - rightSpace - leftSpace)/(g.HEIGHT - topSpace - bottomSpace);
    	g.setUpScale();
		g.paintRaw(image.getGraphics());
		return image;
	}
	//small number of plots, < 6 and with the first plot setting domain and range
	public static BufferedImage drawPlots(double[][] x, double[][] y, int width, int height)
	{
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		GraphDrawerCart g = new GraphDrawerCart("", x[0], y[0]);
		for (int i = 1; i < x.length; i++)
			g.setXY(new GraphDrawerCart.GraphObject(x[i], y[i]), false, false, 1);
       	g.WIDTH = width;
    	g.HEIGHT = height;
    	g.aspect = (double)(g.WIDTH - rightSpace - leftSpace)/(g.HEIGHT - topSpace - bottomSpace);
    	g.setUpScale();
		g.paintRaw(image.getGraphics());
		return image;
	}
	public static BufferedImage drawPlot(double[] x, double[] y, double xmin, double xmax, double ymin, double ymax)
	{
		BufferedImage image = new BufferedImage(DEFAULT_WIDTH, DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
		GraphDrawerCart g = new GraphDrawerCart("", x, y, xmin, xmax, ymin, ymax);
		g.paintRaw(image.getGraphics());
		return image;
	}
	public static BufferedImage drawPlot(double[] x, double[] y, double xmin, double xmax, double ymin, double ymax, int width, int height)
	{
		BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		GraphDrawerCart g = new GraphDrawerCart("", x, y, xmin, xmax, ymin, ymax);
       	g.WIDTH = width;
    	g.HEIGHT = height;
    	g.aspect = (double)(g.WIDTH - rightSpace - leftSpace)/(g.HEIGHT - topSpace - bottomSpace);
    	g.setUpScale();

		g.paintRaw(image.getGraphics());
		return image;
	}
	public static BufferedImage drawPlot(double[][][] plots)
	{
		BufferedImage image = new BufferedImage(DEFAULT_WIDTH, DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
		GraphDrawerCart g = new GraphDrawerCart("", plots);
		g.paintRaw(image.getGraphics());
		return image;
	}
	public static BufferedImage getPlot(double[][] data_x_ys, int modulo)
	{
    	double[] x = data_x_ys[0];
		BufferedImage image = new BufferedImage(DEFAULT_WIDTH, DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
    	GraphDrawerCart drawer = new GraphDrawerCart("", true);
    	if (modulo == 0) modulo = 6;
    	for (int i = 0; i < data_x_ys.length-1; i++)
    	{
    		drawer.setXY(new GraphObject(x, data_x_ys[i+1]), true, true, i);
        	drawer.setColor(drawColors[i%modulo], i);
    	}
		drawer.paintRaw(image.getGraphics());
		return image;
	}
}

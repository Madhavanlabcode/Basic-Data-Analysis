package drawing;

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

import main.CentroidField;
import main.DataManip;
import main.SRAW;
import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFT2DSmall;
import util.fourier.FFTOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.geom.Mask;
import util.matrix.Matrix;
import util.*;

//This will guide the user through the drift-correction procedure.
public class DriftCorrectionConsole extends JFrame implements ActionListener, MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{

	private static final long serialVersionUID = -4335309540571779626L;

	double[][] topograph; int N;
	
	double field[][]; double fmax, fmin, fdelta;
	double[][][] cField;
	boolean complex = false;
	boolean twoComp = false;
	BufferedImage image;
	ColorScale scale;
	ColorScales.MYC2d cScale;
	
	BraggVectorRegularizer reg;
	String braggInfoString = "";
	
	double cfmin, cfmax, cfdelta;
	
	boolean etaPlus = true;
	boolean drawingFT = false;
	
	int imsx = 1024, imsy = 1024;
	int ox = 20, oy = 40;
	int iftox = ox+imsx+10, iftoy = oy;
	
	int sizeRatio;	
	int[] writepoint = {ox, oy + imsy + 10};
	
	int nBragg = 2;
	int currentBragg = 0;
	double[][] bragg = new double [nBragg][2];
	double[][] braggTrue = new double [nBragg][2];
	int[] chosenL = new int [nBragg];
	boolean realSpaceCalc = true;
	
	int currentL = 0;
	
	//Now the good parts
	int nLengths = 100;
	double maxL;
	double[] lengths;
	double[][][][] phaseTopo;
	double[][] gaussMask;
	double[][][][][] expU;
	double[][][][] expUPhase;
	double[][][] uField; //we need to keep both signs //we know that the correct sign is plus.
	double[][] topoAfter; //same as ufield
	int[][][] phaseN;
	double[][][] phaseCont;
	int zoomLevel = 0, zoomFactor = 1;
	
	int WIDTH = 1060, HEIGHT = 1080;
	
	FFT2DSmall fft;
	double[][] ftmag;
	double[][] storedftmags;
	double[][] ftmagdraw;
	double[][][] fftz;
	double[][] iftfield, masterift;
	double[][][] iftz;
	double iftmin, iftmax;
	double fftminmag;
	boolean ftlog = true;
	int[] iftwritepoint = {iftox, iftoy};
	int ftsize;
	
	double ftmax, ftmin;
	double upfrac = 1, downfrac = 0;
	double ftscalemin, ftscalemax;
	
	int currentx, currenty;
	int currftx, currfty; //on screen.
	int currftxp, currftyp; //in the data set.
	int xdown, ydown, xup, yup;
	//The field must be 2^n by 2^n;
	
	int stage = 0;
	
	AtomicCoordinatesSet latt, lattReg; //these are "centered" coordinate sets with lattReg being with regularized Bragg vectors
	double[][][] uRegularization, uBasic;
	
	//For the boxes demarcating the FT:
	Box[] boxes = new Box[2];
	int selectedbox = -1;
	int nboxes = 0;
	Color[] boxcolors = {Color.BLUE, Color.GREEN};
	
	
	
	String dir;
	//I think this no longer works for sizes greater than 512.
	public DriftCorrectionConsole(double[][] field, String dir)
	{
		this.dir = dir;
		Layer temp = Layer.getFreeLayer(field);
		RHKFileOps.doUserFitting(temp);
		topograph = FieldOps.copy(field);
		this.field = FieldOps.copy(field);
		N = topograph.length;
		ftsize = N;
//		ftox = ox + sx + 10; ftoy = oy;
		writepoint = new int[] {ox, oy};
		imsx -= imsx % N;
		imsy = imsx;
		sizeRatio = imsx/N;
		int angDeg = Integer.parseInt(JOptionPane.showInputDialog(null, "Enter the angle for Bragg regularization, in degrees.", 90));
		System.out.println(Math.toRadians(angDeg));
		reg = new BraggVectorRegularizer(Math.toRadians(angDeg), this);
		nLengths = Integer.parseInt(JOptionPane.showInputDialog("User, how many different length scales do you want to try? \r\n" +
				"(Will spread linearly from 0 to a maximum length to be specified later.)"));
		
		
//		box = new Box(0, 0, field.length);
		System.out.println("image");
		fft = new FFT2DSmall(FieldOps.copy(topograph));
		ftmag = new double[N][N];
		fftz = new double [N][N][2];
		ftmagdraw = new double [N][N];
		iftfield = new double [N][N];
		masterift = new double [N][N];
		phaseTopo = new double [nBragg][N][N][2];
		expU = new double [nBragg][nLengths][N][N][2];
		expUPhase = new double [nBragg][nLengths][N][N];
		uField = new double [N][N][2];
		topoAfter = new double[N][N];
		phaseN = new int [nBragg][N][N];
		phaseCont = new double [nBragg][N][N];
		cField = expU[0][0];
		
		doFFT();
		fftminmag = FieldOps.magMin(fftz);
		storedftmags = FieldOps.magnitude(fftz);
		this.field = ftmagdraw;

		JOptionPane.showMessageDialog(null, "Please draw a box around the first Bragg peak. It doens't have to be too small. When satisfied, hit the Space bar.");
		
//		doIFFT();
		image = new BufferedImage(imsx, imsy, BufferedImage.TYPE_INT_RGB);
		setFieldInfo();
		resetColorScale(0, 1);
//		resetFRColorScale(0, 1);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		showWindow();
//		setTitle("Draw boxes with the mouse. x to exclude boxed regions from FT. i to include only boxed region. r to reset. s to save.");

		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		addMouseListener(this);
		JFrame f = new JFrame();
		JPanel p = new JPanel();
		p.setLayout(new BoxLayout(p, BoxLayout.PAGE_AXIS));
		BoundsPanel panel1 = new BoundsPanel(this, f, 0);
		p.add(panel1);
		f.setSize(400, 100);
		f.add(p);
		f.setVisible(true);
	}
	public void doFFT()
	{
		fft.doFFT();
		fft.fHat2(fftz);
		setFTDraw();
		ftmax = ArrayOps.max(ftmagdraw);
		ftmin = ArrayOps.min(ftmagdraw);
	}
	public void doIFFT()
	{
		fft.doIFFT();
		FieldOps.getIndex(fft.f, masterift, 0);
		iftz = FieldOps.copy(fft.f);
		this.iftfield = masterift;
		int i = 0;
		while (i < zoomLevel-1)
			{iftfield = FieldOps.reduce(2, iftfield); i++;}
		iftmax = ArrayOps.max(iftfield);
		iftmin = ArrayOps.min(iftfield);
	}
	
	public void setFTDraw()
	{
		FieldOps.magnitude(fftz, ftmag);
		ftmagdraw = FieldOps.copy(ftmag);
		if (ftlog) FieldOps.log(ftmagdraw);
	}
	
	public void setFieldInfo()
	{
		cfmax = FieldOps.magMax(cField); cfmin = FieldOps.magMin(cField);
		fmax = max(field); fmin = min(field);
		
		fdelta = fmax - fmin;
		cfdelta = cfmax - cfmin;
		System.out.println("[" + fmin + ", " + fmax + "]");
	}
	public void resetColorScale(double downnum, double upnum)
	{
		cScale = new ColorScales.MYC2d(cfmin + cfdelta*upnum, cfmin + cfdelta*downnum, 2*Math.PI); 
		scale = new ColorScales.LinearBRYW(fmin + fdelta*upnum, fmin + fdelta*downnum);
	}
	public void formImage()
	{
		if (complex)
			for (int i = 0; i < imsx; i++)
				for (int j = 0; j < imsy; j++)
					image.setRGB(i, j, cScale.of(cField[i/sizeRatio][j/sizeRatio]).getRGB());
		else
			for (int i = 0; i < imsx; i++)
				for (int j = 0; j < imsy; j++)
					image.setRGB(i, j, scale.of(field[i/sizeRatio][j/sizeRatio]).getRGB());
	}

	public int[] getFieldPoint(int x, int y)
	{
		return new int [] {x - ox, y - oy};
	}
	
	public int[] getQuad(int x, int y)
	{
		int[] r = getFieldPoint(x, y);
		r[0] *= 2; r[0] /= imsx;
		r[1] *= 2; r[1] /= imsy;
		return r;
	}
	public void drawImage(Graphics g)
	{
		g.drawImage(image, ox, oy, null);
	}
	public void paint(Graphics g)
	{
		g.clearRect(0, 0, 2000, 2000);
		drawImage(g);
		g.drawString(scale.toString(), writepoint[0], writepoint[1]);
		if (stage == 0)
			drawBoxes(g);
//		g.setColor(java.awt.Color.BLUE);
//		box.draw(g, this);
	}
	public void drawBoxes(Graphics g)
	{
		for (int i = 0; i < nboxes; i++){
			g.setColor(boxcolors[i]);
			boxes[i].draw(g, this);
//			boxes[i].draw(g, this);
		}
	}
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		repaint();
		setVisible(true);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	public void setLIndex(int value) {
		// TODO Auto-generated method stub
		cField = expU[currentBragg][value];
		field = expUPhase[currentBragg][value];
		currentL = value;
		setFieldInfo();
		resetColorScale(0, 1);
		formImage();
		repaint();
	}
	public void moveBox(int dx, int dy)
	{
		boxes[0].x += dx;
		boxes[0].y += dy;
		repaint();
	}
	public void resizeBox(int dx, int dy)
	{
		boxes[0].dx += dx;
		boxes[0].dy += dy;
		repaint();
	}
	public void doSummation(int braggIndex)
	{
		for (int j = 0; j < N; j++)
			for (int k = 0; k < N; k++){
				phaseTopo[braggIndex][j][k][0] = Math.cos(DataManip.dot(braggTrue[braggIndex], j, k))*topograph[j][k];
				phaseTopo[braggIndex][j][k][1] = -Math.sin(DataManip.dot(braggTrue[braggIndex], j, k))*topograph[j][k];
			}
		for (int j = 0; j < lengths.length; j++)
		{
			setTitle("On Bragg Vector # " + (braggIndex+1) + ", doing calculation " + (j+1) + " out of " + nLengths + ".");
			gaussMask = Mask.getGaussianMask(lengths[j]);
			DataManip.getDeviationGaussDefault(lengths[j], phaseTopo[braggIndex], expU[braggIndex][j], gaussMask);
		}
	}
	public void setBraggTrue(int i)
	{
		double[][] temp = new double [2][2];//for record-keeping.
		temp[0][0] = bragg[i][0]; temp[0][1] = bragg[i][1];
		bragg[i][0] = FieldOps.round(bragg[i][0]);
		bragg[i][1] = FieldOps.round(bragg[i][1]);
		temp[1][0] = bragg[i][0]; temp[1][1] = bragg[i][1];
		bragg[i][0] -= N/2;
		bragg[i][1] -= N/2;
		
		braggTrue[i][0] = bragg[i][0]*2*Math.PI/N;
		braggTrue[i][1] = bragg[i][1]*2*Math.PI/N;
		
		braggInfoString += "Bragg vetor # " + (i+1) + " is " + Printer.vectorP(bragg[i]) + " with respect to the center, which is " + Printer.vectorP(temp[1]) + " in the image.\r\n";
		braggInfoString += "(It was rounded from " + Printer.vectorP(temp[0]) + ".)\r\n";
		braggInfoString += "In inverse pixels its value is " + Printer.vectorP(braggTrue[i]) + ".\r\n";
	}
	public void executeBraggRegularization()
	{
		int ans = JOptionPane.showConfirmDialog(null, "User, do you want to force-regularize the Bragg peaks? (This will also equalize their magnitudes.)", "Force-regularization", JOptionPane.YES_NO_OPTION);
		
		if (ans == JOptionPane.NO_OPTION) {reg = null; return;}
		reg.regularize();
		uRegularization = new double [N][N][2];
		uBasic = new double [N][N][2];
		
	}
	public void doExpUCalc()
	{
		for (int i = 0; i < nBragg; i++)
			doSummation(i);
		
		setTitle("Populating the phase array....");
		for (int i = 0; i < nBragg; i++)
			for (int j = 0; j < nLengths; j++)
				for (int m = 0; m < N; m++)
					for (int n = 0; n < N; n++)
						expUPhase[i][j][m][n] = FieldOps.atan(expU[i][j][m][n][0], expU[i][j][m][n][1]);
		setTitle("Done.");
		
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
		if (stage == 0){
			if (arg0.getKeyChar() == 'a')// || arg0.getKeyChar() == 'A')
				moveBox(-1, 0);
			if (arg0.getKeyChar() == 'd')// || arg0.getKeyChar() == 'D')
				moveBox(+1, 0);
			if (arg0.getKeyChar() == 'w')// || arg0.getKeyChar() == 'W')
				moveBox(0, -1);
			if (arg0.getKeyChar() == 's')// || arg0.getKeyChar() == 's')
				moveBox(0, +1);
			if (arg0.getKeyChar() == 'A')// || arg0.getKeyChar() == 'A')
				resizeBox(-1, 0);
			if (arg0.getKeyChar() == 'D')// || arg0.getKeyChar() == 'D')
				resizeBox(+1, 0);
			if (arg0.getKeyChar() == 'W')// || arg0.getKeyChar() == 'W')
				resizeBox(0, -1);
			if (arg0.getKeyChar() == 'S')// || arg0.getKeyChar() == 's')
				resizeBox(0, +1);
		}
		if (arg0.getKeyChar() == 'c' && stage > 0)
		{
			complex = !complex;
			formImage();
			repaint();
		}
		if (arg0.getKeyChar() == 'f' && stage > 1)
		{
			if (field == ftmagdraw)
				field = topoAfter;
			else if (field == topoAfter)
				field = ftmagdraw;
			setFieldInfo();
			resetColorScale(0, 1);
			formImage();
			repaint();
		}
		if (arg0.getKeyChar() == 'o')
			setFieldInfo();
		if (arg0.getKeyChar() == ' ')
		{
			switch(stage)
			{
			case 0:
				bragg[currentBragg] = CentroidField.centroid(ftmag, boxes[0].x, boxes[0].dx, boxes[0].y, boxes[0].dy, true);
				currentBragg++;
				if (currentBragg == nBragg){
					stage++;
					maxL = Double.parseDouble(JOptionPane.showInputDialog("Please enter the maximum length scale to be used in the algorithm, in pixels (real space)."));
					lengths = ArrayOps.generateArrayNotInclLower(0, maxL, nLengths);
					setBraggTrue(currentBragg-1);
					latt = AtomicCoordinatesSet.generateCentered(bragg, N);
					reg.setSecondBraggPeak(braggTrue[1]);
					System.out.println("Initial Bragg vectors " + Printer.vectorP(braggTrue[0]) + " and " + Printer.vectorP(braggTrue[1]));
					executeBraggRegularization();
					System.out.println("Rounded Bragg vectors " + Printer.vectorP(braggTrue[0]) + " and " + Printer.vectorP(braggTrue[1]));
					long t0 = System.currentTimeMillis();
					doExpUCalc();
					currentBragg = 0;
					setLIndex(0);
					long t1 = System.currentTimeMillis();
					
					JOptionPane.showMessageDialog(null, "(The calculation took " + ((t1-t0)/1000) + " seconds.)\r\n" +
							"Now, vary the slider until the image before you is satisfactory. Then press Spacebar.\r\n" +
							"To swicth between the phase view and the complex view, press 'c'");
				}
				else{
					//-------------------------------------------------This is for the Bragg regularizer
					setBraggTrue(currentBragg-1);
					reg.setFirstBraggPeak(braggTrue[0]);
					boxes[0].move(FieldOps.round(reg.braggSetRot[0][1][0]/(2*Math.PI/N))+N/2, FieldOps.round(reg.braggSetRot[0][1][1]/(2*Math.PI/N))+N/2, this);
					repaint();
					//--------------------------------------------------------------------
					JOptionPane.showMessageDialog(null, "Now draw one aroung the second Bragg peak. Press Spacebar when done. \r\n" +
							"If you want to force-orthogonalize the Bragg peaks, draw your box around the one near the place the box has moved to.");
					
				}
				break;
			case 1:
				chosenL[currentBragg] = currentL;
				currentBragg++;
				if (currentBragg == nBragg)
				{
					stage++;
					currentBragg = 0;
//					JOptionPane.showMessageDialog(null, "Now the computation will finish. ");
					
					for (int i = 0; i < nBragg; i++)
					{
						FieldOps.putPhaseSteps(expUPhase[i][chosenL[currentBragg]], N/2, N/2, phaseN[i]);
						FieldOps.putAddedPhaseSteps(expUPhase[i][chosenL[currentBragg]], phaseN[i], 2*Math.PI, phaseCont[i]);
					}
					FieldOps.putU(phaseCont[0], phaseCont[1], braggTrue, 1, uField);
					FieldOps.subtractAvg(uField);
					//--------------------------------------------------This is for Bragg regularization
					if (reg != null)
					{
						uBasic = FieldOps.copy(uField);
						FieldOps.putUField(latt, lattReg, uRegularization);
						FieldOps.add(uBasic, uRegularization, uField);
					}
					//---------------------------------------------------------------------------------
					
					
					FieldOps.applyUFieldSmooth(topograph, uField, 2, 2, topoAfter);
					FieldOps.changeZeroToAverage(topoAfter);
					
					field = topoAfter;
					cField = uField;
					setFieldInfo();
					resetColorScale(0, 1);
					formImage();
					repaint();
					fft.setF(topoAfter);
					doFFT();
					JOptionPane.showMessageDialog(null, "Press the Space bar to save everything. In the save menu, just type the 'base' name. \r\n" +
							"To switch to viewing the U-field and back, press 'c'. \r\n" +
							"To view the FFT press 'f' in the real mode.");
				}
				else
				{
					JOptionPane.showMessageDialog(null, "Repeat, for the other component of the u-field");
					setLIndex(0);
				}
				break;
			case 2:
				//save everything then quit
				File f = FileOps.selectSave(new JFileChooser(dir));
				String s = f.toString();
				String braggWord, lengthWord;
				RHKFileOps.write3Files(s + "after", topoAfter);
				double[][][] fftz = new double [N][N][2];
				
				
				ColumnIO.writeString(braggInfoString, f.toString() + "Bragg Peak Info.txt");
				double[][][] fHat2 = new double [N][N][2];
				for (int i = 0; i < nBragg; i++)
				{
					braggWord = "b" + i;
					FileOps.writeBINandImage(s + "phaseSteps" + braggWord, phaseN[i]);
					FileOps.writeBINandImage(s + "phaseCont_" + braggWord, phaseCont[i]);
					for (int j = 0; j < nLengths; j++)
					{
						lengthWord = "L" + j;
						FileOps.writeBINandImage(s + "phaseRough_" + braggWord + "_" + lengthWord, expUPhase[i][j]);
						SRAW.writeImage(s + "localSum_" + braggWord + "_" + lengthWord, expU[i][j], true);
						FFTOps.putFFT(expU[i][j], fftz, false);
						FieldOps.shift(fftz, fHat2);
						ColumnIO.writeBinsPolar(fHat2, s + "localSumFFT_" + braggWord + "_" + lengthWord);
						SRAW.writeImage(s + "localSumFFT_" + braggWord + "_" + lengthWord, FFTOps.getImageCent(fftz, true, false, fftminmag));
					}
					String lengthInfo = "You picked the lengths " + lengths[chosenL[0]] + " (" + chosenL[0] + ") " + " and " + lengths[chosenL[1]] + " (" + chosenL[1] + ") " + " respectively. The calculation was done for the length scales \r\n";
					lengthInfo += Printer.arrayVertical(lengths);
					ColumnIO.writeString(lengthInfo, s + "length_info.txt");
				}
				if (reg == null)
				{
					ColumnIO.writeString(latt.toString(), f.toString() + "Final lattice.txt");
					ColumnIO.writeBinsPolar(uField, s + "uField");
					SRAW.writeImage(s + "uField", uField, true);
				}
				else
				{
					ColumnIO.writeString(lattReg.toString(), f.toString() + "Final lattice.txt");
					ColumnIO.writeString(latt.toString(), f.toString() + "Unregularized lattice.txt");
					//Now I want to write "after" shots without refugularization:
					FieldOps.applyUFieldSmooth(topograph, uBasic, 2, 2, topoAfter);
					FieldOps.changeZeroToAverage(topoAfter);
					RHKFileOps.write3Files(s + "_nonlinear_only", topoAfter);
					FieldOps.applyUFieldSmooth(topograph, uRegularization, 2, 2, topoAfter);
					FieldOps.changeZeroToAverage(topoAfter);
					RHKFileOps.write3Files(s + "_linear_only", topoAfter);
					/*
					 * ^^Save results with and without the linear regularizing U-field.
					 */
//					ColumnIO.writeBinsPolar(uRegularization, s + "uLinear");
					ColumnIO.writeTable(FieldOps.copy(uRegularization, 0), s + "uLinear_x.txt");
					ColumnIO.writeTable(FieldOps.copy(uRegularization, 1), s + "uLinear_y.txt");
					SRAW.writeImage(s + "uLinear", uRegularization, true);
//					ColumnIO.writeBinsPolar(uBasic, s + "uBasic");
					ColumnIO.writeTable(FieldOps.copy(uBasic, 0), s + "uBasic_x.txt");
					ColumnIO.writeTable(FieldOps.copy(uBasic, 1), s + "uBasic_y.txt");
					SRAW.writeImage(s + "uNonlinear", uBasic, true);
					ColumnIO.writeBinsPolar(uField, s + "uField");
//					ColumnIO.writeTable(FieldOps.copy(uField, 0), s + "uField_x.txt");
//					ColumnIO.writeTable(FieldOps.copy(uField, 1), s + "uField_y.txt");
					SRAW.writeImage(s + "uTotal", uField, true);
				}
				
				//Write the bicubic applications of U:
				File bicDir = new File(f.getParent() + "\\BiCubic\\");
				bicDir.mkdir();
				String sbc = bicDir.toString() + "\\" + f.getName();
				FieldOps.applyUFieldBiCubic(topograph, uField, topoAfter);
				RHKFileOps.write3Files(sbc + "after", topoAfter);
				
				FieldOps.applyUFieldBiCubic(topograph, uBasic, topoAfter);
				FieldOps.changeZeroToAverage(topoAfter);
				RHKFileOps.write3Files(sbc + "_nonlinear_only", topoAfter);
				FieldOps.applyUFieldBiCubic(topograph, uRegularization, topoAfter);
				FieldOps.changeZeroToAverage(topoAfter);
				RHKFileOps.write3Files(sbc + "_linear_only", topoAfter);
				

				
				
//				JOptionPane.showMessageDialog(null, "Done.");
				System.out.println();
				System.exit(0);
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
		xdown = arg0.getX();
		ydown = arg0.getY();
	}

	@Override
	public void mouseReleased(MouseEvent arg0) {
		// TODO Auto-generated method stub
		xup = arg0.getX();
		yup = arg0.getY();
		if(selectedbox == -1)
		{
			int[] dr = {(xup-xdown)/sizeRatio, (yup-ydown)/sizeRatio};
			createNewBox(this.ftrCent(xdown, ydown), dr);
			repaint();
		}
		
	}

	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	
	//returns the coordinates within the fft, with the center thereof being the origin.
	int[] ftrCent(int x, int y)
	{
		return new int[] {(x-ox)/sizeRatio - N/2, (y-oy)/sizeRatio - N/2};
	}
	
	//the coordinates r must be with respect to the center of the fft, in fft pixels
	public void createNewBox(int[] r, int [] dr)
	{
		//must create two boxes, one at r and one at (-r-dr)
		boxes[0] = new Box(r[0]+N/2, r[1]+N/2, dr[0], dr[1]);
		if (nboxes == 0)
			nboxes++;
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
	 public static void suppressTo(double[] z, double finalmag)
	 {
		 double mag = Complex.mag(z);
		 z[0] *= finalmag/mag;
		 z[1] *= finalmag/mag;
	 }
	 
	 public void reset()
	 {
		 fft.reset(FieldOps.copy(topograph));
		 doFFT();
		 doIFFT();
		 repaint();
		 nboxes = 0;
		 
	 }
	 public void save()
	 {
		JFileChooser fc = new JFileChooser(dir);
			if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION){
//				SRAW.writeImage(fc.getSelectedFile().toString() + "fft", ftmag);
//				SRAW.writeImage(fc.getSelectedFile().toString() + "fft_scaled", ftmag, ftscale);
//				SRAW.writeImage(fc.getSelectedFile().toString() + "ift", masterift);
//				SRAW.writeImage(fc.getSelectedFile().toString() + "ift_scaled", masterift, iftscale);
//				SRAW.writeImage(fc.getSelectedFile().toString() + "trans", newTopo);
//				ColumnIO.writeBin(masterift, fc.getSelectedFile().toString() + "ift.dat");
				ColumnIO.writeTable(masterift, fc.getSelectedFile().toString() + ".txt");
//				ColumnIO.writeString(info + "kmag used = " +  kmag + "\r\nscale factor = " + scale, fc.getSelectedFile().toString() + "transLog.txt");
				JOptionPane.showMessageDialog(null, "Done");
			}
	 }
	 
	 private static class Box
	 {
		 int x, y, dx, dy;

		public Box(int x, int y, int dx, int dy) {
			this.x = x;
			this.y = y;
			this.dx = dx;
			this.dy = dy;
		}
		
		public void draw(Graphics g, DriftCorrectionConsole parent)
		{
			g.drawRect(parent.ox + x*parent.sizeRatio, parent.oy + y*parent.sizeRatio, dx*parent.sizeRatio, dy*parent.sizeRatio);
		}
		public void move(int x, int y, DriftCorrectionConsole parent)
		{
			this.x = x > 0 ? x+dx < parent.N ? x : parent.N - dx : 0; this.y = y > 0 ? y+dy < parent.N ? y : parent.N - dy : 0;
		}
		@SuppressWarnings("unused")
		public void setSize(int sizex, int sizey)
		{
			this.dy = sizex;
			this.dy = sizey;
		}
		@SuppressWarnings("unused")
		public int[] getR()
		{
			return new int[] {x, y};
		}
		@SuppressWarnings("unused")
		public int[] getDR()
		{
			return new int[] {dx, dy};
		}
	 }
	 public static class BoundsPanel extends JPanel implements ChangeListener
	 {
		private static final long serialVersionUID = -8428501987212326257L;
		DriftCorrectionConsole parent;
		 JFrame frame;
		 public JSlider min, max;
		 public JSlider ls;
		 int parentmember;
		 public BoundsPanel(DriftCorrectionConsole parent, JFrame frame, int member)
		 {
				this.frame = frame;
				this.parent = parent;
				parentmember = member;
				this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
				this.setBorder(new LineBorder(Color.GRAY));
				min = new JSlider(0, 1000, 0);
				min.addChangeListener(this);
				max = new JSlider(0, 1000, 1000);
				max.addChangeListener(this);
				ls = new JSlider(0, parent.nLengths-1, 0);
				ls.addChangeListener(this);
				
				add(ls);
				add(min);
				add(max);
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
			if (ls.getValueIsAdjusting() && parent.stage == 1)
			{
				parent.setLIndex(ls.getValue());
			}
			
			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
			parent.formImage();
			parent.repaint();
			frame.setTitle("Min: " + parent.ftscalemin + "   Max: " + parent.ftscalemax);
		}
	}
	 
	 //This class is to regularize the Bragg vectors to a square lattice. This is to be done so as to approximately minimize the maximum distance of ANY bragg point from the actual initial Bragg peak.
	 public static class BraggVectorRegularizer{
		
		public double angle = Math.PI/2;
		public DriftCorrectionConsole parent;
		
		double[][] rotMat = new double [2][2];
		double[][] rotMatT = new double [2][2];
		
		double[][][] braggSetRot = new double [2][2][2];
		double[][][] unitBraggSetRot = new double [2][2][2];
		//[0] is the first user-defined bragg peak and its +angle cousin. [1] = second user-defined Bragg peak and its -angle cousin 
		public BraggVectorRegularizer(double angle,
				DriftCorrectionConsole parent) {
			super();
			this.angle = angle;
			this.parent = parent;
			Matrix.putRotationMatrix(angle, rotMat);
			Matrix.putRotationMatrix(-angle, rotMatT);
		}
		
		public void setFirstBraggPeak(double[] bragg)
		{
				for (int j = 0; j < 2; j++)
					braggSetRot[0][0][j] = bragg[j];
				braggSetRot[0][1] = Matrix.getProductWith(rotMat, bragg);
		}
		public void setSecondBraggPeak(double[] bragg)
		{
				for (int j = 0; j < 2; j++)
					braggSetRot[1][1][j] = bragg[j];
				braggSetRot[1][0] = Matrix.getProductWith(rotMatT, bragg);
		}
		
		public void regularize()
		{
			double mag1 = Complex.mag(braggSetRot[0][0]);
			double mag2 = Complex.mag(braggSetRot[1][1]);
			double magAvg = (mag1+mag2)/2;
			
			System.out.println(mag1 + "\t" + mag2 + "\t" + magAvg);
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					unitBraggSetRot[i][j] = Distance.unitVector(braggSetRot[i][j]);
			
			double angle1 = Complex.phase(unitBraggSetRot[0][0]);
			double angle2 = Complex.phase(unitBraggSetRot[1][0]);
			//In the case where angle1 or angle2 is close to 0, the other may be on the other side of 2pi. In that case special care must be taken:
			if (Math.abs(angle1-angle2) > Math.PI)
				if (angle1 > Math.PI) angle1 -= 2*Math.PI;
				else if (angle2 > Math.PI) angle2 -= 2*Math.PI;
			
			double angleAvg = (angle1+angle2)/2;
			System.out.println(angle1 +"\t"+ angle2 + "\t"+angleAvg);
			
			double[][] finalBragg = new double[2][2];
			finalBragg[0][0] = magAvg*Math.cos(angleAvg);
			finalBragg[0][1] = magAvg*Math.sin(angleAvg);
			finalBragg[1][0] = magAvg*Math.cos(angleAvg+angle);
			finalBragg[1][1] = magAvg*Math.sin(angleAvg+angle);
			double[][] temp = new double [2][2];
			
			int even = JOptionPane.showConfirmDialog(null, "Should we force the Bragg peak to have an even pixel number (useful for rt2 order)?", "Even or Odd", JOptionPane.YES_NO_OPTION);
			if (even == JOptionPane.YES_OPTION){
				temp[0] = roundBraggIntEven(finalBragg[0], parent.N);
				temp[1] = roundBraggIntEven(finalBragg[1], parent.N);
			}
			else
			{
				temp[0] = roundBraggInt(finalBragg[0], parent.N);
				temp[1] = roundBraggInt(finalBragg[1], parent.N);
			}
			for (int i = 0; i < 2; i++){
				parent.braggInfoString += "Upon regularization, Bragg vetor # " + (i+1) + " was converted to " + Printer.vectorP(temp[i]) + " with respect to the center.\r\n";
//				parent.braggInfoString += "In inverse pixels its value is " + Printer.vectorP(finalBragg[i]) + ".\r\n";
			}
			parent.lattReg = AtomicCoordinatesSet.generateCentered(finalBragg, parent.N);
//			parent.braggTrue = finalBragg;
		}
		
		//this returns the temp (pixel) variable for record-keeping purposes
		public static double[] roundBraggInt(double[] vec, int N)
		{
			double[] temp = new double[] {vec[0]/(2*Math.PI/N), vec[1]/(2*Math.PI/N)};
			temp[0] = FieldOps.round(temp[0]);
			temp[1] = FieldOps.round(temp[1]);
//			System.out.println(Printer.vectorP(temp));
			vec[0] = temp[0];//*(2*Math.PI/N);
			vec[1] = temp[1];//*(2*Math.PI/N);
			return temp; 
			
		}
		public static double[] roundBraggIntEven(double[] vec, int N)
		{
			double[] temp = new double[] {vec[0]/(2*Math.PI/N), vec[1]/(2*Math.PI/N)};
			temp[0] = FieldOps.roundEven(temp[0]);
			temp[1] = FieldOps.roundEven(temp[1]);
//			System.out.println(Printer.vectorP(temp));
			vec[0] = temp[0];//*(2*Math.PI/N);
			vec[1] = temp[1];//*(2*Math.PI/N);
			return temp; 
			
		}
	 }

	 public static void main(String[] args)
	 {
//	 	String dir = "C:\\data\\lawlerhex\\run153topo15_830\\";
//		dir = "C:\\data\\lintrans\\";
		 Topomap.setStdDir();
	 	JFileChooser fc = new JFileChooser(Topomap.stddir);
		new DriftCorrectionConsole(Layer.openFree(fc).data, fc.getCurrentDirectory().toString() + "\\");
	 }

}

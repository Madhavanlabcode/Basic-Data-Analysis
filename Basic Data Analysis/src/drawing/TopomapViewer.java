package drawing;

import image.ImageEditing;
import impurity.PointImp;

import java.awt.Toolkit;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import main.SRAW;
import schrodinger.MovieMaker;
import util.ArrayOps;
import util.Complex;
import util.ExportUtil;
import util.FieldOps;
import util.ImpurityUtil;
import util.NumFormat;
import util.Printer;
import util.SpectraUtil;
import util.TopomapUtil;
import util.color.ColorScale1D;
import util.color.ColorScale2d;
import util.color.ColorScaleHolder;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.PointSpectra;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFT2DSmall;
import util.fourier.FFT3D_Wrapper;
import util.fourier.FFTOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.regression.FFTCutFitter;
import util.regression.FFTCutFitter.EQPair;
import util.robot.Robo;

public class TopomapViewer extends JFrame implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	
	
	Topomap t;
	int N;
	double fmax, fmin, fdelta;
	double[][][] fftmag = null;
	boolean showFFT = false;
	double[][] drawField;
	boolean real = true;
	
	GraphDrawerCart spec;
	double[] spectrum;
	double[] twoV, twoVy;
	
	static JFileChooser fc;
	
	BufferedImage image;
	boolean refresh = true;
	
//	ColorScale scale;
	ColorScale2d cscale;
	double scalemin, scalemax;
	
	int defaultSize = 512;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 60;
	
	int[] writepoint = {ox + sx + 50, oy + 10};
	int linesize = 15;
	
	int zoomLevel = 0, zoomFactor = 1;
	int sizeratio = 1;
	
	int WIDTH = 1600, HEIGHT = 1080;
	
	Point mouseDownCompCoords=null;
	
	//para is the current index
	int para = 0;
	double imageOverField = 1;

	int currentx, currenty;
	int currenti, currentj;
	double currentid, currentjd;
	double currentr, currentphi;
	int calcx = 0, calcy = 0;
	String dir;
	String name = null;
	double[] pbounds;
	SliderPanel s;
	int snpts;
	
	//line cut feature
	LineCutDrawer lc = null;
	StripCutTool sct = null;
	double[][] line = null; //the two points
	public int nearestEndIndex = 0;
	double d1, d2;
	
	//color scale change feater.
//	int currentCScale = 0;
	ColorScaleHolder csh;
	private boolean changeScaleWithLayer=true;
	private boolean ftlog = true;
	private boolean refreshFFT = false;
	
	double[] min, max, delta;
	double[] realmin, realmax, realdelta;
	double[] fftmin, fftmax, fftdelta;
	
	ArrayList<PointImp> imps = null;
	
	AtomicCoordinatesSet latt = null;

	
	FFTCutFitter fitter;
	LayerViewer lv;
	
	Robo rob = new Robo();
	
	JMenuBar menuBar;
	
	public TopomapViewer(Topomap t, String dir,  int size)
	{
		//Setting up menu bar
		menuBar = new JMenuBar();
		
	    // File Menu, F - Mnemonic
	    JMenu fileMenu = new JMenu("File");
	    fileMenu.setMnemonic(KeyEvent.VK_F);
	    menuBar.add(fileMenu);
	    
	    //File->save map as a .bin
        JMenuItem saveMI = new JMenuItem("Save map as a .bin");
        fileMenu.add(saveMI);
        class saveAL implements ActionListener{
        	Topomap t;
        	saveAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		Topomap.writeBIN(t, fc);
    		}
        }
        saveMI.addActionListener(new saveAL(t));
        
        //File-> save current layer
        JMenuItem saveLayerMI = new JMenuItem("Save current layer as a layer .bin");
        fileMenu.add(saveLayerMI);
        class saveLayerAL implements ActionListener{
        	Topomap t;
        	saveLayerAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			int pixcoord = JOptionPane.showConfirmDialog(null, " YES: Use original length coordinates in layer \r\n NO: Use pixel length coordinates in layer", "Length Units", JOptionPane.YES_NO_OPTION);
    			if (pixcoord == JOptionPane.YES_OPTION)
    				Layer.writeBIN(t.getLayer(para));
    			else Layer.writeBIN(t.getLayerPixels(para));
    		}
        }
        saveLayerMI.addActionListener(new saveLayerAL(t));
        
        //File-> save all layers as .bin
        JMenuItem saveAllLayersMI = new JMenuItem("Save all layers as layer .bin's");
        fileMenu.add(saveAllLayersMI);
        class saveAllLayersAL implements ActionListener{
        	Topomap t;
        	saveAllLayersAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
				String base = FileOps.selectSave(fc).toString();
				String insert = JOptionPane.showInputDialog("Enter what you with to name each layer.", "layer");
				for (int i = 0; i < t.nlayers; i++)
					Layer.writeBIN(t.getLayer(i), base + "_"+insert+"_" + i + ".bin");
    		}
        }
        saveAllLayersMI.addActionListener(new saveAllLayersAL(t));
	    
        //File->Copy image to clipboard
        JMenuItem imageToClipMI = new JMenuItem("Paste image to clipboard");
        fileMenu.add(imageToClipMI);
        imageToClipMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
	      	  	BufferedImage export = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
				Graphics g = export.getGraphics();
				refresh = true;
				paint(g);
				ImageEditing.copyToClipboard(ImageEditing.getSubset(export, ox, sx*sizeratio, oy, sy*sizeratio));
      	  	}
        });
        
        //File->Copy spectrum to clipboard
        JMenuItem specToClipMI = new JMenuItem("Paste spectrum image to clipboard");
        fileMenu.add(specToClipMI);
        specToClipMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		ImageEditing.copyToClipboard((BufferedImage)spec.dbimage);
      	  	}
        });
        
        //File->Copy color scale to clipboard
        JMenuItem scaleToClipMI = new JMenuItem("Paste color scale to clipboard");
        fileMenu.add(scaleToClipMI);
        scaleToClipMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		ImageEditing.copyToClipboard(csh.getScale().getScaleImage(8));
      	  	}
        });
        
        //File->spectrum to clipboard as text
        JMenuItem specToTsvMI = new JMenuItem("Copy spectrum as tab-separated values");
        fileMenu.add(specToTsvMI);
        class specToTsvAL implements ActionListener{
        	Topomap t;
        	specToTsvAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			String[] plines = new String [t.nlayers];
    			for (int i = 0; i < t.nlayers; i++)
    			{
    				plines[i] = t.v[i] + "\t" + spectrum[i];
    				System.out.println(plines[i]);
    				Printer.copyToClipboard(plines);
    			}
    		}
        }
        specToTsvMI.addActionListener(new specToTsvAL(t));
        
        //File->fourier transform to tab-separated values
        JMenuItem fftToTsvMI = new JMenuItem("Copy Fourier transform as tab-separated table");
        fileMenu.add(fftToTsvMI);
        class fftToTsvAL implements ActionListener{
        	Topomap t;
        	fftToTsvAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
            	//Write a spectrum of the real and imaginary parts of the FFT as a function of energy, and copy it to the cliboard.
    			double[][][] fftz = new double [t.nlayers][t.nx][t.ny];
    			int i = currenti;
    			int j = currentj;
    			String[] lines = new String[t.nlayers];
    			for (int k = 0; k < t.nlayers; k++){
    				FFTOps.putFFT(t.data[k], fftz, true);
    				lines[k] = "" + t.v[k] + "\t"  + Complex.mag(fftz[i][j]) + "\t" + fftz[i][j][0] + "\t" + fftz[i][j][1] + "\t" + Complex.phase(fftz[i][j]);
    				System.out.println(lines[k]);
//    				System.out.println(i + "\t" + j + "\t" + currenti + "\t" + currentj);
    			}
    			Printer.copyToClipboard(lines);
    		}
        }
        fftToTsvMI.addActionListener(new fftToTsvAL(t));
        
        //File->Layers to movie
        JMenuItem layersToMovieMI = new JMenuItem("Layers to movie utility");
        fileMenu.add(layersToMovieMI);
        class layersToMovieAL implements ActionListener{
        	Topomap t;
        	layersToMovieAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
				int typeOfGIF = Integer.parseInt(JOptionPane.showInputDialog("Select the output type:\r\n" +
						" 0 - GIF with topo and FFT inset\r\n" +
						" 1 - GIF with FFT only (and voltage)\r\n" + 
						" 2 - Large BMP with layers as tiles\r\n" + 
						" 3 - Movie from BMP to AVI\r\n" + 
						" 4 - Tiled BMP with movie settings\r\n" +
						" 5 - Movie from BMP to AVI with both\r\n" +
						" 6 - Movie from BMP to AVI with another map\r\n" +
						" 7 - Movie from BMP to AVI about important FFT point\r\n" + 
						" 8 - Tiled BMP, movie settings, important FFT point\r\n"));
				File f;
				switch(typeOfGIF){
					case 0:
						f = FileOps.selectSave(fc);
						Layer topo = Layer.openFree(fc);
						RHKFileOps.doFitting(topo, RHKFileOps.getUserFittingChoice());
						if (showFFT)
							ExportUtil.exportGIFForPPTSlide(t, topo, csh.getScale(), null, false, true, true, 500, f.toString(), csh.getCurrentCS());
						else
							ExportUtil.exportGIFForPPTSlide(t, topo, null, csh.getScale(), false, true, true, 500, f.toString(), csh.getCurrentCS());
						break;
					case 1:
						f = FileOps.selectSave(fc);
						TopomapUtil.writeGIFWithVoltage(t, f.toString(), Integer.parseInt(JOptionPane.showInputDialog("How much to blow it up by?")), csh.getCurrentCS());
						break;
					case 2:
						f = FileOps.selectSave(fc);
						int nPerLine = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
//						BufferedImage image = ImageEditing.createTileImage(showFFT ? fftmag : t.data, csh.getScale(), nPerLine, t.v);
						BufferedImage image = ImageEditing.createTileImage(showFFT ? fftmag : t.data, csh.getScale(), nPerLine, t.v, csh.getCurrentCS());
						SRAW.writeImage(f.toString(), image);
						ImageEditing.copyToClipboard(image);
						break;
					case 3:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						boolean useCurrentScale = JOptionPane.showConfirmDialog(null, "Lock the color scale?") == JOptionPane.YES_OPTION; 
						double minperc = useCurrentScale ? 0 : Double.parseDouble(JOptionPane.showInputDialog("Minimum distribution cutoff?"));
						double maxperc = useCurrentScale ? 1 : Double.parseDouble(JOptionPane.showInputDialog("Maximum distribution cutoff?"));
						ColorScale1D scale = null;
						if (useCurrentScale) {
							String[] tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							scale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						int[] rF = new int [2];
						if (t.nx == t.ny) {
							int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
							rF[0] = resizeFactor; rF[1] = resizeFactor;
						}
						else {
							int resizeFactorX = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor? X", "" + 1));
							int resizeFactorY = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor? Y", "" + 1));
							rF[0] = resizeFactorX; rF[1] = resizeFactorY;
					
						}
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						BufferedImage[] stack;
						if (t.nx == t.ny)
							stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, rF[0], useCurrentScale ? scale : null);
						else
							stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, rF, useCurrentScale ? scale : null);
							
						BufferedImage[] substack = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++){
							substack[i-imin] = stack[i];
						}
						if (imps != null)
							for (int i = 0; i < substack.length; i++)
							{
								Graphics g = substack[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).pixelPos;
									drawPlus(g, (int)(r[0]*rF[0]), (int)(r[1]*rF[1]), 5);
								}
							}
						
						for (int i = 0; i < substack.length; i++)
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), substack[i]);
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substack.length-1, 5));
						break;
					}
					case 4:{
						int nPerLine1 = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
						boolean useCurrentScale = JOptionPane.showConfirmDialog(null, "Lock the color scale?") == JOptionPane.YES_OPTION; 
						double minperc = useCurrentScale ? 0 : Double.parseDouble(JOptionPane.showInputDialog("Minimum distribution cutoff?"));
						double maxperc = useCurrentScale ? 1 : Double.parseDouble(JOptionPane.showInputDialog("Maximum distribution cutoff?"));
						ColorScale1D scale = null;
						if (useCurrentScale) {
							String[] tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							scale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						
						BufferedImage[] stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, resizeFactor, useCurrentScale ? scale : null);
						BufferedImage[] substack = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++)
							substack[i-imin] = stack[i];
						if (imps != null)
							for (int i = 0; i < substack.length; i++)
							{
								Graphics g = substack[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).pixelPos;
									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
								}
							}
						BufferedImage tiles = ImageEditing.createTileImage(substack, nPerLine1, Color.WHITE);
						ImageEditing.copyToClipboard(tiles);
						SRAW.writeImage(FileOps.selectSave(fc).toString(), tiles);
						
						break;
					}
					case 5:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						
						String firstInput = JOptionPane.showInputDialog("Enter the half-length of the box you want to capture.", "" + t.nx/2);
						int width = Integer.parseInt(firstInput)*2;
						
						double minpercR = 0;
						double maxpercR = 1;
						double minpercK = 0;
						double maxpercK = 1;
						boolean useFixedScaleReal = JOptionPane.showConfirmDialog(null, "Use fixed color scale in real space?") == JOptionPane.YES_OPTION;
						ColorScale1D realScale = null;
						if (useFixedScaleReal){
							String[] tok;
							if (showFFT) tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.",  "" + ArrayOps.min(t.data) + "," + ArrayOps.max(t.data)).split(",");
							else tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							
							realScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						else{
						minpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Minimum distribution cutoff?"));
						maxpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Maximum distribution cutoff?"));
						}

						ColorScale1D kScale = null;
						boolean useFixedScaleK = JOptionPane.showConfirmDialog(null, "Use fixed color scale in k space?") == JOptionPane.YES_OPTION;
						if (useFixedScaleK){
							String[] tok;
							if (!showFFT) tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + ArrayOps.min(fftmag) + "," + ArrayOps.max(fftmag)).split(",");
							else tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							
							kScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						else{
						minpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Minimum distribution cutoff?"));
						maxpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Maximum distribution cutoff?"));
						}
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
//						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						
						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						BufferedImage[] stackR = ImageEditing.createImageArray(t.data, csh.getCurrentCS(), t.v, label, minpercR, maxpercR, resizeFactor, realScale);
						BufferedImage[] stackK = ImageEditing.createImageArray(fftmag, csh.getCurrentCS(), t.v, 0, minpercK, maxpercK, resizeFactor, kScale);
						BufferedImage[] substackR = new BufferedImage[imax-imin];
						BufferedImage[] substackK = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++){
							substackR[i-imin] = stackR[i];
							substackK[i-imin] = stackK[i];
							if (width < t.nx/2){
								substackR[i-imin] = ImageEditing.getSubset(substackR[i-imin], (t.nx-width)/2, width, (t.ny-width)/2, width);
								substackK[i-imin] = ImageEditing.getSubset(substackK[i-imin], (t.nx-width)/2, width, (t.ny-width)/2, width);;
							}
						}
						if (imps != null)
							for (int i = 0; i < substackR.length; i++)
							{
								Graphics g = substackK[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).pixelPos;
									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
								}
							}
						
						BufferedImage out;
						for (int i = 0; i < substackR.length; i++){
							if (showFFT)
								out = ImageEditing.weldHorizontal(substackK[i], substackR[i]);
							else
								out = ImageEditing.weldHorizontal(substackR[i], substackK[i]);
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), out);
						}
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substackR.length-1, 5));
						break;
					}
					case 6:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						
						double minpercR = 0;
						double maxpercR = 1;
						double minpercK = 0;
						double maxpercK = 1;
						Topomap other = Topomap.open(fc);
						boolean useFixedScaleReal = JOptionPane.showConfirmDialog(null, "Use fixed color scale in real space?") == JOptionPane.YES_OPTION;
						ColorScale1D realScale = null;
						if (useFixedScaleReal){
							String[] tok;
							tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							
							realScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						else{
						minpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Minimum distribution cutoff?"));
						maxpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Maximum distribution cutoff?"));
						}

						ColorScale1D kScale = null;
						boolean useFixedScaleK = JOptionPane.showConfirmDialog(null, "Use fixed color scale in other space?") == JOptionPane.YES_OPTION;
						int cschoice = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice of color scale index.", "" + csh.getCurrentCS()));
						if (useFixedScaleK){
							String[] tok;
							tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + ArrayOps.min(other.data) + "," + ArrayOps.max(other.data)).split(",");
							
							kScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), cschoice);
						}
						else{
						minpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Minimum distribution cutoff?"));
						maxpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Maximum distribution cutoff?"));
						}
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
//						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						
//						BufferedImage[] stackR = ImageEditing.createImageArray(t.data, csh.getCurrentCS(), t.v, 2, minpercR, maxpercR, resizeFactor, realScale);
//						BufferedImage[] stackK = ImageEditing.createImageArray(other.data, cschoice, t.v, 0, minpercK, maxpercK, resizeFactor, kScale);
//						BufferedImage[] substackR = new BufferedImage[imax-imin];
//						BufferedImage[] substackK = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++){
							BufferedImage imr = ImageEditing.createSingleImage(t.data, csh.getCurrentCS(), t.v, 2, minpercR, maxpercR, resizeFactor, realScale, i);
							BufferedImage imk = ImageEditing.createSingleImage(other.data, cschoice, t.v, 0, minpercK, maxpercK, resizeFactor, kScale, i);
//							if (imps != null)
//							{
//								Graphics g = imr.getGraphics();
//								g.setColor(csh.getUnusedColor());
//								double[] r;
//								for (int j = 0; j < imps.size(); j++)
//								{
//									r = imps.get(j).pixelPos;
//									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
//								}
//							}
						
							BufferedImage out = ImageEditing.weldHorizontal(imr, imk);
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), out);
							
						}
						
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, imax-imin-1, 5));
						break;
					}
					case 7:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						BufferedImage[] substack = getImageArrayHighSymPt();
						for (int i = 0; i < substack.length; i++)
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), substack[i]);
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substack.length-1, 5));
						break;
					}
					case 8:
					{
						int nPerLine1 = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
						BufferedImage[] substack = getImageArrayHighSymPt();
						BufferedImage tiles = ImageEditing.createTileImage(substack, nPerLine1, Color.WHITE);
						ImageEditing.copyToClipboard(tiles);
						SRAW.writeImage(FileOps.selectSave(fc).toString(), tiles);
						break;
					}
				}	
        	}
        }
        layersToMovieMI.addActionListener(new layersToMovieAL(t));
        
        //File-> save as comma-separated values
        JMenuItem saveAsCsvMI = new JMenuItem("Save all layers as comma-separated value files");
        fileMenu.add(saveAsCsvMI);
        class saveAsCsvAL implements ActionListener{
        	Topomap t;
        	saveAsCsvAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			String path = FileOps.selectSave(fc).toString();
    			for (int i = 0; i < t.nlayers; i++)
    				ColumnIO.writeTableCSV(t.data[i], path + "_" + i);
    		}
        }
        saveAsCsvMI.addActionListener(new saveAsCsvAL(t));
        
        //File-> Save the map as a PointSpectra obect
        JMenuItem saveAsPointSpecMI = new JMenuItem("Save map as a PointSpectra obect");
        fileMenu.add(saveAsPointSpecMI);
        class saveAsPointSpecAL implements ActionListener{
        	Topomap t;
        	saveAsPointSpecAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		PointSpectra ps = TopomapUtil.getTimeOrderedPointSpectra(t, 0);
				PointSpectra.writeBIN(ps, fc);
    		}
        }
        saveAsPointSpecMI.addActionListener(new saveAsPointSpecAL(t));
        
        //Display menu
	    JMenu dispMenu = new JMenu("Display");
	    menuBar.add(dispMenu);
	    
	    //Display->forward in color scale
        JMenuItem forwardScaleMI = new JMenuItem("Increment color scale (shortcut 'q')");
        dispMenu.add(forwardScaleMI);
        forwardScaleMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
    			csh.incrementScaleIndex();
    			formImage();
    			refresh = true;
    			repaint();
      	  	}
        });
	    
	    //Display->backward in color scale
        JMenuItem backwardScaleMI = new JMenuItem("Decrement color scale (shortcut 'e')");
        dispMenu.add(backwardScaleMI);
        backwardScaleMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
    			csh.decrementScaleIndex();
    			formImage();
    			refresh = true;
    			repaint();
      	  	}
        });
	    
	    //Display->set color scale to extremes
	    JMenuItem setColorExtremesMI = new JMenuItem("Set color scale by extremes");
	    dispMenu.add(setColorExtremesMI);
	    class setColorExtremesAL implements ActionListener{
	    	Topomap t;
	    	setColorExtremesAL(Topomap tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
				if (showFFT)
					csh.reBoundColorScale(ArrayOps.min(fftmag), ArrayOps.max(fftmag));
				else
					csh.reBoundColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data));
				resetGraphics(true, 0, 1);
	    	}
	    }
	    setColorExtremesMI.addActionListener(new setColorExtremesAL(t));
	    
	    //Display->toggle changing scale with layer
        JMenuItem autoscaleMI = new JMenuItem("Toggle changing scale with layer");
        dispMenu.add(autoscaleMI);
        autoscaleMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		changeScaleWithLayer = !changeScaleWithLayer;
      	  	}
        });
        
	    //Display->refresh graphics
        JMenuItem refreshMI = new JMenuItem("Refresh picture (shortcut spacebar)");
        dispMenu.add(refreshMI);
        refreshMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
    			refresh = true;
    			repaint();
      	  	}
        });
        
        //Display->higher layer
	    JMenuItem upLayerMI = new JMenuItem("Go up a layer (shortcut '=')");
	    dispMenu.add(upLayerMI);
	    class upLayerAL implements ActionListener{
	    	Topomap t;
	    	upLayerAL(Topomap tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
				if (para < t.nlayers-1){
					para++;
					twoV[0] = t.v[para];
					twoV[1] = t.v[para];
					spec.repaint();
					
					resetGraphics(changeScaleWithLayer, ((double)s.min.getValue()/1000), ((double)s.max.getValue()/1000));
					s.frame.setTitle("" + t.v[para] + "   (Layer " + para + ")");
				}
	    	}
	    }
	    upLayerMI.addActionListener(new upLayerAL(t));
        
        //Display->lower layer
	    JMenuItem downLayerMI = new JMenuItem("Go down a layer (shortcut '-')");
	    dispMenu.add(downLayerMI);
	    class downLayerAL implements ActionListener{
	    	Topomap t;
	    	downLayerAL(Topomap tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
				if (para > 0){ 
					para--;
					twoV[0] = t.v[para];
					twoV[1] = t.v[para];
					spec.repaint();
					
					resetGraphics(changeScaleWithLayer, ((double)s.min.getValue()/1000), ((double)s.max.getValue()/1000));
					s.frame.setTitle("" + t.v[para] + "   (Layer " + para + ")");
				}
	    	}
	    }
	    downLayerMI.addActionListener(new downLayerAL(t));
        
	    // Fourier Transform menu
	    JMenu fourierMenu = new JMenu("Fourier Transforms");
	    menuBar.add(fourierMenu);
	    
        //Fourier->take transform
        JMenuItem fftMI = new JMenuItem("Take/revert FFT (shortcut 'f')");
        fourierMenu.add(fftMI);
    	class fftAL implements ActionListener{
    		private Topomap t;
    		fftAL(Topomap tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    	    	//take Fourier transform.
      			if (fftmag == null)
      			{
      				refreshFFT = false;
      				fftmag = new double [t.nlayers][t.nx][t.ny];
      				for (int i = 0; i < t.nlayers; i++)
      				{
      					setTitle("Doing FFT: " + i);
      					FFTOps.obtainFFTmagCent(t.data[i], fftmag[i]);
      					if (ftlog)
      						FieldOps.log(fftmag[i]);
      				}
      				recalculateFFTBounds();
      				setTitle("Done");
      			}
      			else if (refreshFFT)
      			{
      				refreshFFT = false;
      				for (int i = 0; i < t.nlayers; i++)
      				{
      					setTitle("Doing FFT: " + i);
      					FFTOps.obtainFFTmagCent(t.data[i], fftmag[i]);
      					if (ftlog)
      						FieldOps.log(fftmag[i]);
      				}
      				recalculateFFTBounds();
      			}
      			showFFT = !showFFT;
      			if (showFFT) {
      				min = fftmin;
      				max = fftmax;
      				delta = fftdelta;
      				spec.setRange(ArrayOps.min(fftmag), ArrayOps.max(fftmag));
      			}
      			else{
      				min = realmin;
      				max = realmax;
      				delta = realdelta;
      				spec.setRange(ArrayOps.min(t.data), ArrayOps.max(t.data));
      			}
      			
      			if (lc != null)
      			{
      				if (showFFT)
      					lc.setDataset(fftmag);
      				else
      					lc.setDataset(t.data);
      			}
      			else if (sct != null)
      			{
      				if (showFFT)
      					sct.setData(fftmag);
      				else
      					sct.setData(t.data);
      			}
      			resetGraphics(true, 0, 1);
           }
        }
        fftMI.addActionListener(new fftAL(t));
        
        //Fourier->take complex transform
        JMenuItem complexFftMI = new JMenuItem("View complex Fourier transform");
        fourierMenu.add(complexFftMI);
        class complexFftAL implements ActionListener{
        	Topomap t;
        	String dir;
        	complexFftAL(Topomap tPrime, String dirPrime){
        		t=tPrime;
        		dir=dirPrime;
        	}
        	public void actionPerformed(ActionEvent e) {
      			double[][][][] results = new double [2][t.nlayers][t.nx][t.ny];
      			FFT2DSmall fft;
      			double mag, magmin = Double.MAX_VALUE, newmag, logmag;
      			for (int k = 0; k < t.nlayers; k++)
      			{
      				fft = FFTOps.obtainFFT(t.data[k]);
      				for (int i = 0; i < t.nx; i++)
      					for (int j = 0; j < t.ny; j++)
      					{
      						results[0][k][i][j] = fft.getFHat2Re(i, j);
      						results[1][k][i][j] = fft.getFHat2Im(i, j);
      						if (ftlog){
      							mag = Complex.mag(results[0][k][i][j], results[1][k][i][j]);
      							magmin = Math.min(mag, magmin);
      						}
      					}
      			}
      			if (ftlog)
      				for (int k = 0; k < t.nlayers; k++)
      					for (int i = 0; i < t.nx; i++)
      						for (int j = 0; j < t.ny; j++)
      						{
      							results[0][k][i][j] *= 1/magmin; //divide out magmin so the minimum magnitude is one.
      							results[1][k][i][j] *= 1/magmin;
      							newmag = Complex.mag(results[0][k][i][j], results[1][k][i][j]);
      							logmag = Math.log(newmag);
      							results[0][k][i][j] *= logmag/newmag; //divide out magmin so the minimum magnitude is one.
      							results[1][k][i][j] *= logmag/newmag;
      						}
      		
      			new TopomapViewer_complex2(Topomap.newTopomap(t, results[0]), Topomap.newTopomap(t, results[1]),dir, 512);
      		}	
        }
        complexFftMI.addActionListener(new complexFftAL(t,dir));
		
        //Fourier->toggle log scale in transform
        JMenuItem toggleLogMI = new JMenuItem("Toggle log scale in transform");
        fourierMenu.add(toggleLogMI);
        toggleLogMI.addActionListener(new ActionListener(){
        	public void actionPerformed(ActionEvent e){
    			ftlog = !ftlog;
    			refreshFFT = true;
    			refresh=true;
    			repaint();
        	}
        });
        
        //Fourier->show 3d FFT
        JMenuItem ThreeDFftMI = new JMenuItem("Replace map with 3D Fourier Transform");
        fourierMenu.add(ThreeDFftMI);
        class ThreeDFftAL implements ActionListener{
        	Topomap t;
        	ThreeDFftAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e) {
    			t.data = FFT3D_Wrapper.getFFTMagCentXY(t.data, JOptionPane.showConfirmDialog(null, "Take the log?") == JOptionPane.YES_OPTION);
      		}	
        }
        ThreeDFftMI.addActionListener(new ThreeDFftAL(t));

        //Fourier->Show spatial modulation of a Q-vector
        JMenuItem spatialLockinMI = new JMenuItem("Show spatial modulation of a q-vector");
        fourierMenu.add(spatialLockinMI);
        class spatialLockinAL implements ActionListener{
        	Topomap t;
        	String dir;
        	spatialLockinAL(Topomap tPrime, String dirPrime){
        		t=tPrime;
        		dir=dirPrime;
        	}
        	public void actionPerformed(ActionEvent e) {
      			double[][][] Rea = new double [t.nlayers][t.nx][t.ny];
      			double[][][] Ima = new double [t.nlayers][t.nx][t.ny];
      			double[][][] Mag = new double [t.nlayers][t.nx][t.ny];
      			
      			//get Q
      			String qString = JOptionPane.showInputDialog("Q-vector components (relative to center, in units of FFT pixels, comma-separated)");
				double Qx = Double.parseDouble(qString.split(",")[0].trim());
				double Qy = Double.parseDouble(qString.split(",")[1].trim());
      			
      			//multiply t by e^iQ.r
				System.out.println("Demodulating, with wave numbers " + Qx + ", " + Qy + ".");
				for(int i=0;i<t.nlayers;i++)
					for(int j=0;j<t.nx;j++)
						for(int k=0;k<t.ny;k++){
							Rea[i][j][k]=t.data[i][j][k]*Math.cos(Qx*j*2*Math.PI/t.nx+Qy*k*2*Math.PI/t.ny);
							Ima[i][j][k]=t.data[i][j][k]*Math.sin(Qx*j*2*Math.PI/t.nx+Qy*k*2*Math.PI/t.ny);
						}
      			
      			//get blur length
				double wavelength = t.nx/Math.sqrt(Math.pow(Qx,2)+Math.pow(Qy,2));
      			double blurLength = Double.parseDouble(JOptionPane.showInputDialog("Averaging length scale, in units of pixels. (Wavelength is " + wavelength + ")."));
      			//gaussian blur
      			System.out.print("Low-pass filtering with length scale " + blurLength + ". Be patient.");
      			for(int i=0;i<t.nlayers;i++){
      				Rea[i]=FieldOps.gaussSmooth(Rea[i], blurLength);
      				System.out.println("  Halfway there.");
      				Ima[i]=FieldOps.gaussSmooth(Ima[i], blurLength);
      				Mag[i]=FieldOps.magnitude(Rea[i], Ima[i]);
      			}
      			
      		
      			new TopomapViewer_complex2(Topomap.newTopomap(t, Rea), Topomap.newTopomap(t, Ima),dir, 512);
      			new TopomapViewer(Topomap.newTopomap(t, Mag),dir, 512);
      		}	
        }
        spatialLockinMI.addActionListener(new spatialLockinAL(t,dir));
        
        //Selecting data menu
        JMenu selectionMenu = new JMenu("Selection tools");
        menuBar.add(selectionMenu);
        
        //Selection->line cut toggle
        JMenuItem toggleLinecutMI = new JMenuItem("Toggle linecut tool");
        selectionMenu.add(toggleLinecutMI);
        class toggleLinecutAL implements ActionListener{
        	Topomap t;
        	toggleLinecutAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			if (lc == null)
    			{	
    				if (line == null) line = new double[][] {{0, 0}, {t.nx, t.ny}};
    				//2015.11.18 changed 2 points to more (2d argument)
    				//2nd argument controls how many cuts to use
    				lc = new LineCutDrawer(t, 16, line);
    				lc.pullLegend();
    				lc.setFC(fc);
    			}
    			else {
    				lc.dispose();
    				lc = null;
    			}
    			refresh = true;
    			repaint();
        	}
        }
        toggleLinecutMI.addActionListener(new toggleLinecutAL(t));

        //Selection->mirror data across linecut
        JMenuItem mirrorLinecutMI = new JMenuItem("Mirror data across linecut");
        selectionMenu.add(mirrorLinecutMI);
        class mirrorLinecutAL implements ActionListener{
        	Topomap t;
        	mirrorLinecutAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			if (line != null){
    				double[][][] result = new double[t.nlayers][t.nx][t.ny];
    				double[] mirror=new double[2];
    				mirror[0]=line[0][1]-line[1][1];
    				mirror[1]=line[1][0]-line[0][0];
    				double length = Math.pow(Math.pow(mirror[0],2)+Math.pow(mirror[1],2),0.5);
    				mirror[0] = mirror[0]/length;
    				mirror[1] = mirror[1]/length;
    				System.out.println("Mirror vector: " + mirror[0] + "," + mirror[1]);
    				for(int a=0;a<t.nlayers;a++){
    					for(int b=0;b<t.nx;b++){
    						for(int c=0;c<t.ny;c++){
    							double distance = (b-line[0][0])*mirror[0]+(c-line[0][1])*mirror[1];
    							int bPrime = (int)(-2*distance*mirror[0]+b);
    							int cPrime = (int)(-2*distance*mirror[1]+c);
    							if(bPrime>0 && bPrime<t.nx && cPrime>0 && cPrime<t.ny){
    								result[a][b][c] = (t.data[a][b][c]+t.data[a][bPrime][cPrime])/2;
    							}
    							else{
    								result[a][b][c]=t.data[a][b][c];
    							}
    						}
    					}
    				}
    				t.data = result;
    			}
    			resetGraphics(true, 0, 1);
        	}
        }
        mirrorLinecutMI.addActionListener(new mirrorLinecutAL(t));
        
        //Selection->how to mark a point
        selectionMenu.add(new JMenuItem("Use 'm' key to mark a point"));
        
        //Selection->read a coordinates file and mark them
        JMenuItem loadImpsMI = new JMenuItem("Load a coodinates file and mark them");
        selectionMenu.add(loadImpsMI);
        loadImpsMI.addActionListener(new ActionListener(){
        	public void actionPerformed(ActionEvent e){
    			if (imps == null && lc == null) {
    				PointImp[] der = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
    				imps = new ArrayList<PointImp>();
    				for (int i = 0; i < der.length; i++)
    					imps.add(der[i]);
    				refresh = false;
    				repaint();
    			}
    			else if(lc==null)
    				imps = null;
        	}
        });
        
        //Selection->save the average spectra at marked points
        JMenuItem saveImpsMI = new JMenuItem("Save spectra at marked points");
        selectionMenu.add(saveImpsMI);
        class saveImpsAL implements ActionListener{
        	Topomap t;
        	saveImpsAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		if(imps != null)
        		{	
        			PointSpectra spec = ImpurityUtil.getSpectraAt(ImpurityUtil.getFromList(imps), t, false, 1);
        			PointSpectra.writeBIN(spec, fc);
        			String[] lines = spec.toLines();
        			for (int i = 0; i < lines.length; i++)
        				System.out.println(lines[i]);
        		}
        	}
        }
        saveImpsMI.addActionListener(new saveImpsAL(t));
        
        //Selection->Save average spectra around marked points
        JMenuItem saveAvgImpsMI = new JMenuItem("Save average spectra around marked points");
        selectionMenu.add(saveAvgImpsMI);
        class saveAvgImpsAL implements ActionListener{
        	Topomap t;
        	saveAvgImpsAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		TopomapUtil.writeAverageSpectraAroundImps(t, imps, Double.parseDouble(JOptionPane.showInputDialog("Enter the gaussian radius in pixels")), fc);
        	}
        }
        saveAvgImpsMI.addActionListener(new saveAvgImpsAL(t));
        
        //Split the map into bins based on an existing layer
        JMenuItem splitMI = new JMenuItem("Split the map into bins based on an existing layer");
        selectionMenu.add(splitMI);
        class splitAL implements ActionListener{
        	Topomap t;
        	String dir;
        	splitAL(Topomap tPrime, String dirPrime){
        		t=tPrime;
        		dir=dirPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		boolean bool = JOptionPane.showConfirmDialog(null, "Apply a boolean filter first?") == JOptionPane.YES_OPTION;
    			if (!bool){
    				Topomap binshist = TopomapUtil.splitTopomapByBinsOfLayerFast_autoHist(t, dir, name, Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins")), fc);
//    				Topomap.writeBIN(binshist, fc);
    				if (t.nlayers > 1){
    				TopomapViewer tv = new TopomapViewer(binshist, dir, getGoodSize(binshist.nx));
    				tv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    				}
    			}
    			else
    			{
    				boolean[][] thebools = ImageEditing.loadFromImage(FileOps.selectOpen(fc).toString(), Color.WHITE);
    				int nb = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
    				Layer binstuff = TopomapUtil.splitTopomapByBinsOfLayer(t, dir, name, thebools, nb, fc);
    				ColorScale1D binStuff = ColorScales.getNew(0, nb-1, csh.getCurrentCS());
    				BufferedImage binBool = ImageEditing.getBufferedImage(binstuff.data, binStuff);
    				FieldOps.negate(thebools);
    				ImageEditing.setTrueToColor(thebools, binBool, csh.getUnusedColor());
    				ImageEditing.copyToClipboard(binBool);
    				
    			}
        	}
        }
        splitMI.addActionListener(new splitAL(t,dir));
        
        //Selection->Split the map into smoothed bins based on an existing layer
        JMenuItem splitSmoothMI = new JMenuItem("Split the map into smoothed bins based on an existing layer");
        selectionMenu.add(splitSmoothMI);
        class splitSmoothAL implements ActionListener{
        	Topomap t;
        	splitSmoothAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		{
    				int nbins = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
    				int[][] bins = FieldOps.getPercentileBinsForField(Layer.open(fc).data, nbins);
    				double[][][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins, Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels (to soften edges)")));
    				Topomap smw = new Topomap(smoothWeight, ArrayOps.generateArray(0, 1, nbins), t.x, t.y, null);
    				File save = FileOps.selectSave(fc);
    				for (int i = 0; i < nbins; i++)
    				{
    					Topomap.writeBIN(TopomapUtil.spatialFilter(t, smoothWeight[i]), save.toString() + "_" + i + ".bin");
    				}
    				Topomap.writeBIN(smw, save.toString() + "_WeightsForBins.bin");
    			}
        	}
        }
        splitSmoothMI.addActionListener(new splitSmoothAL(t));
        
        //Selection->Split the given layer against smoothed bins and write a map of the subsets
        JMenuItem splitSmoothWriteMI = new JMenuItem("Split the given layer against smoothed bins and write a map of the subsets");
        selectionMenu.add(splitSmoothWriteMI);
        class splitSmoothWriteAL implements ActionListener{
        	Topomap t;
        	splitSmoothWriteAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
				int nbins1 = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
				int[][] bins1 = FieldOps.getPercentileBinsForField(Layer.open(fc).data, nbins1);
				double[][][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins1, Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels (to soften edges)")));
				double[][][] smoothParts = new double [smoothWeight.length][t.nx][t.ny];
				for (int i = 0; i < nbins1; i++)
				{
					smoothParts[i] = FieldOps.spatialFilter(t.data[para], smoothWeight[i]);
				}
				double[] binsArray = ArrayOps.generateArrayNotInclUpper(0, nbins1, nbins1);
				Topomap result = new Topomap(smoothParts, binsArray, t.x, t.y, null);
				Topomap weights = new Topomap(smoothWeight, binsArray, t.x, t.y, null);
				
				Topomap.writeBIN(result, fc);
				Topomap.writeBIN(weights, fc);
        	}
        }
        splitSmoothWriteMI.addActionListener(new splitSmoothWriteAL(t));
        
        //Selection->Split the map layers against the layers of another 3D dataset
        JMenuItem splitByMapMI = new JMenuItem("Split the map layers against the smoothed layers of another 3D dataset");
        selectionMenu.add(splitByMapMI);
        class splitByMapAL implements ActionListener{
        	Topomap t;
        	splitByMapAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		int nbins1 = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
				Topomap binsource = Topomap.open(fc);
				int[][][] bins1 = new int[t.nlayers][][];
				double[] binsArray = ArrayOps.generateArrayNotInclUpper(0, nbins1, nbins1);
				double L1 = Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels (to soften edges)"));
				int savebin = JOptionPane.showConfirmDialog(null,  "Save the bins?");
				boolean writeEachLayerVsBin = JOptionPane.showConfirmDialog(null, "Would you also like to write a map for each layer, as a function of bin?") == JOptionPane.YES_OPTION;
				File fs = FileOps.selectSave(fc);
				double[][][] smoothWeight = new double [nbins1][t.nx][t.ny];
				double[][][][] smoothParts = new double [nbins1][t.nlayers][t.nx][t.ny];
				for (int j = 0; j < t.nlayers; j++){
					bins1[j] = FieldOps.getPercentileBinsForField(binsource.data[j], nbins1);
					smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins1[j], L1);
					for (int i = 0; i < nbins1; i++)
					{
						smoothParts[i][j] = FieldOps.spatialFilter(t.data[j], smoothWeight[i]);
					}
				}
				if (savebin == JOptionPane.YES_OPTION)
				{
					double[][][] bins2 = new double [t.nlayers][t.nx][t.ny];
					for (int i = 0; i < t.nlayers; i++) bins2[i] = ArrayOps.toDouble(bins1[i]);
					
					Topomap.writeBIN(Topomap.newTopomap(t,bins2), fs.toString() + "bins_unsmoothed.bin");
				}
				if (writeEachLayerVsBin)
				{
					for (int i = 0; i < t.nlayers; i++)
					{
						double[][][] layeri = new double [nbins1][t.nx][t.ny];
						for (int j = 0; j < nbins1; j++)
							layeri[j] = smoothParts[j][i];
						Topomap.writeBIN(new Topomap(layeri, binsArray, t.x, t.y, null), fs.toString() + "layer_" + i + ".bin");					
					}
				}
				
				
				for (int i = 0; i < nbins1; i++)
					Topomap.writeBIN(Topomap.newTopomap(t, smoothParts[i]), fs.toString() + "_" + i + ".bin");
				JOptionPane.showMessageDialog(null, "Done.");
        	}
        }
        splitByMapMI.addActionListener(new splitByMapAL(t));
        
        //Selection->Save a cropped copy of this map
        JMenuItem saveCroppedMI = new JMenuItem("Save a cropped copy of this map");
        selectionMenu.add(saveCroppedMI);
        class saveCroppedAL implements ActionListener{
        	Topomap t;
        	saveCroppedAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		TopomapUtil.saveCropped(t);
        	}
        }
        saveCroppedMI.addActionListener(new saveCroppedAL(t));

        //Selection->Mask using a mask map
        JMenuItem maskMaskMapMI = new JMenuItem("Spatially filter all layers using a mask map");
        selectionMenu.add(maskMaskMapMI);
        class maskMaskMapAL implements ActionListener{
        	Topomap t;
        	maskMaskMapAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		String[] options = {"Mask is a topomap", "Mask is a layer"};
        		switch(JOptionPane.showOptionDialog(null,"Choose mask type","TopomapViewer",JOptionPane.DEFAULT_OPTION,JOptionPane.PLAIN_MESSAGE,null,options,"Mask topomap")){
        		case 0:
            		Topomap spf = Topomap.open(fc);
    				for (int i = 0; i < t.nlayers; i++)
    					t.data[i] = FieldOps.spatialFilter(t.data[i], spf.data[i]);
        			break;
        		case 1:
            		Layer msk = Layer.openFree(fc);
    				for (int i = 0; i < t.nlayers; i++)
    					t.data[i] = FieldOps.spatialFilter(t.data[i], msk.data);
        			break;
        		}
        	}
        }
        maskMaskMapMI.addActionListener(new maskMaskMapAL(t));
        
        
        //Analysis menu
        JMenu analysisMenu = new JMenu("Analysis");
        menuBar.add(analysisMenu);
        
        //Analysis->fft fitter
        JMenuItem fftFitMI = new JMenuItem("Open FFT fitting tool");
        analysisMenu.add(fftFitMI);
        class fftFitAL implements ActionListener{
        	Topomap t;
        	fftFitAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		if(fitter == null)
        		{
        			if (latt == null) latt = new AtomicCoordinatesSet(FileOps.openText(fc));
        			fitter = new FFTCutFitter(t.data, latt);
        		}
        	}
        }
        fftFitMI.addActionListener(new fftFitAL(t));
        
        //Analysis-> work function fit
        JMenuItem wfFitMI = new JMenuItem("Do a work-function fitting");
        analysisMenu.add(wfFitMI);
        class wfFitAL implements ActionListener{
        	Topomap t;
        	wfFitAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		TopomapUtil.fitEachSpectrumExponential(t, "");
        	}
        }
        wfFitMI.addActionListener(new wfFitAL(t));
        
        //Analysis->Print or view correlation with another layer
        JMenuItem correlationMI = new JMenuItem("Print or view correlation with another layer");
        analysisMenu.add(correlationMI);
        class correlationAL implements ActionListener{
        	Topomap t;
        	correlationAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		if (JOptionPane.showConfirmDialog(null, "Open a map instead?") == JOptionPane.YES_OPTION){
    				Topomap ls = Topomap.open(fc);
    				String[] lines = new String [t.nlayers];
    				for (int i = 0; i < t.nlayers; i++){
    					lines[i] = "" + t.v[i];
    					for (int j = 0; j < ls.nlayers; j++)
    					{
    						lines[i] += "\t" + FieldOps.correlation(ls.data[j], t.data[i]);
    					}
    					System.out.println(lines[i]);
    				}
        		}
    			else{
    				Layer l = Layer.open(fc);
    				System.out.println();
    				for (int i = 0; i < t.nlayers; i++)
    					System.out.println("" + t.v[i] + "\t" + FieldOps.correlation(l.data, t.data[i]));
    			}
        	}
        }
        correlationMI.addActionListener(new correlationAL(t));
        
        //Analysis->Get a simple histogram of the spectra
        JMenuItem specHistMI = new JMenuItem("Get a simple histogram of the spectra");
        analysisMenu.add(specHistMI);
        class specHistAL implements ActionListener{
        	Topomap t;
        	String dir;
        	specHistAL(Topomap tPrime, String dirPrime){
        		t=tPrime;
        		dir=dirPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		boolean isDefault = JOptionPane.showConfirmDialog(null, "Use Default settings?") == JOptionPane.YES_OPTION;
    			
    			Layer hist = null;
    			if (isDefault)
    				hist = TopomapUtil.getSpectralDistributionBasic(1, t.nlayers, t);
    			else
    			{
    				String[] tok = JOptionPane.showInputDialog("Enter the limits of the histogram separated by commas.").split(",");
    				double lower = Double.parseDouble(tok[0]);
    				double upper = Double.parseDouble(tok[1]);
    				int nbins = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins", t.nlayers));
    				hist = TopomapUtil.getSpectralDistributionBasicLimited(1, nbins, lower, upper, t);
    			}
    			int size = t.nlayers;
    			while (size*((double)hist.ny/t.nlayers) < 512)
    				size*=2;
    			
    			lv = new LayerViewer(hist, dir, size);
    			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        	}
        }
        specHistMI.addActionListener(new specHistAL(t,dir));
        
        //Analysis->Replace each spectrum by its derivative
        JMenuItem specDerivMI = new JMenuItem("Replace each spectrum by its derivative");
        analysisMenu.add(specDerivMI);
        class specDerivAL implements ActionListener{
        	Topomap t;
        	specDerivAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			for (int i = 0; i < t.nx; i++)
    				for (int j = 0; j < t.ny; j++)
    				{
    					double[] tempSpec = ArrayOps.getDerivative(t.getSpectrum(i, j));
    					for (int k = 0; k < t.nlayers; k++)
    						t.data[k][i][j] = tempSpec[k];
    				}
        	}
        }
        specDerivMI.addActionListener(new specDerivAL(t));
        
        //Analysis->Get a histogram of all the values
        JMenuItem totalHistMI = new JMenuItem("Get a histogram of all the values");
        analysisMenu.add(totalHistMI);
        class totalHistAL implements ActionListener{
        	Topomap t;
        	totalHistAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			String input1 = JOptionPane.showInputDialog("Indices to include (exclude upper one)?", "" + 0 + "," + t.nlayers);
    			int[] is = new int [] {Integer.parseInt(input1.split(",")[0].trim()),Integer.parseInt(input1.split(",")[1].trim())};
    			int nbins = Integer.parseInt(JOptionPane.showInputDialog("How many bins? (Default = Sqrt(number of points)", "" + (int)Math.sqrt((is[1]-is[0])*t.nx*t.ny)));
    			double[] bins = ArrayOps.generateArrayInclBoth(ArrayOps.min(t.data), ArrayOps.max(t.data), nbins);
    			int[] values = ArrayOps.getHistogram(FieldOps.getArray(t.data, is[0], is[1]), bins);
    			GraphDrawerCart.plotGraph(bins, values);
        	}
        }
        totalHistMI.addActionListener(new totalHistAL(t));
        
        //Analysis->Get a simple gap map
        JMenuItem gapMapMI = new JMenuItem("Simple gap map");
        analysisMenu.add(gapMapMI);
        class gapMapAL implements ActionListener{
        	Topomap t;
        	gapMapAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			Layer gapmap = TopomapUtil.getSimplestGapMap(t, Double.parseDouble(JOptionPane.showInputDialog("Enter the cutoff.")));
    			LayerViewer.show(gapmap, 1024, false);
    			//plot now the histogram of the gap map properly.
    			double dv = (t.v[1]-t.v[0]);
    			double[] gapbins = new double [(int)((ArrayOps.max(gapmap.data) - ArrayOps.min(gapmap.data))/dv)];
    			for (int i = 0; i < gapbins.length; i++)
    				gapbins[i] = -dv/2+i*dv;
    			int[] gaphist = ArrayOps.getHistogram(FieldOps.getArray(gapmap.data), gapbins);
    			double[] gaphistRelative = new double [gaphist.length];
    			
    			for (int i = 0; i < gapbins.length; i++){
    				gaphistRelative[i] = ((double)gaphist[i]/(double)(t.nx*t.ny));
    				System.out.println("" + i + "\t" + (i*dv) + "\t" + gaphistRelative[i]);
    			}
    			GraphDrawerCart.plotGraph(gapbins, gaphistRelative);
        	}
        }
        gapMapMI.addActionListener(new gapMapAL(t));

        //Analysis->Take radial average
        JMenuItem radialAvgMI = new JMenuItem("Save the radial average of each layer as csv");
        analysisMenu.add(radialAvgMI);
        class radialAvgAL implements ActionListener{
        	Topomap t;
        	radialAvgAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		int nBins = Math.min(Integer.parseInt(JOptionPane.showInputDialog("Enter number of radial bins to put the data in")),(t.x.length*2/3));
        		double[][] total= new double[t.nlayers][nBins];
        		int[][] count = new int[t.nlayers][nBins];
        		for(int i=0;i<nBins;i++){
        			for(int j=0;j<t.nlayers;j++){
               			count[j][i]=0;
               			total[j][i]=0.0;
        			}
        		}
        		
        		double binSize = Math.sqrt(t.x.length * t.x.length + t.y.length * t.y.length)/nBins/2;
        		System.out.println("binSize: " + binSize);
        		
        		if(showFFT){
        			for(int i=0;i<t.nlayers;i++){
	        			for(int j=0;j<t.x.length;j++){
	        				for(int k=0;k<t.y.length;k++){
	        					total[i][Math.min((int)(Math.sqrt((j-t.x.length/2)*(j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize),nBins-1)]+=fftmag[i][j][k];
	        					count[i][Math.min((int)(Math.sqrt((j-t.x.length/2)*(j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize),nBins-1)]++;
	        				}
	        			}
        			}
        		}
        		else{
        			for(int i=0;i<t.nlayers;i++){
	        			for(int j=0;j<t.x.length;j++){
	        				for(int k=0;k<t.y.length;k++){
	        					total[i][Math.min((int)(Math.sqrt((j-t.x.length/2)*(j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize),nBins-1)]+=t.data[i][j][k];
	        					count[i][Math.min((int)(Math.sqrt((j-t.x.length/2)*(j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize),nBins-1)]++;
	        				}
	        			}
        			}
        		}
        		
        		double[][] radialAvg = new double[t.nlayers][nBins];
        		for(int i=0;i<nBins;i++){
        			for(int j=0;j<t.nlayers;j++){
        				radialAvg[j][i]=total[j][i]/count[j][i];
        			}
        		}
        		System.out.println("Took radial average");
        		ColumnIO.writeTableCSV(radialAvg, FileOps.selectSave(fc).toString());
        	}
        }
        radialAvgMI.addActionListener(new radialAvgAL(t));
        
        //Manipulation menu
        JMenu manipMenu = new JMenu("Manipulation");
        menuBar.add(manipMenu);
        
        //Manipulation->Do something to each layer
        JMenuItem eachLayerMI = new JMenuItem("Do something to each layer");
        manipMenu.add(eachLayerMI);
        class eachLayerAL implements ActionListener{
        	Topomap t;
        	eachLayerAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			int o3 = RHKFileOps.getUserFittingChoice();
    			RHKFileOps.doFittingTopo(t, o3);
        	}
        }
        eachLayerMI.addActionListener(new eachLayerAL(t));
        
        //Manipulation-> symmetrize
        JMenuItem symmMI = new JMenuItem("Symmetrize");
        manipMenu.add(symmMI);
        class symmAL implements ActionListener{
        	Robo rob;
        	symmAL(Robo robPrime){
        		rob=robPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		rob.typeChar('M');
        	}
        }
        symmMI.addActionListener(new symmAL(rob));
        
        //Manipulation->suppress marked points to the average
        JMenuItem suppressMI = new JMenuItem("Suppress masked areas to the average");
        manipMenu.add(suppressMI);
        class suppressAL implements ActionListener{
        	Topomap t;
        	suppressAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			Layer mask = Layer.open(fc);
    			for (int i = 0; i < t.nlayers; i++)
    				FieldOps.changeToWithoutMask(t.data[i], mask.data);
        	}
        }
        suppressMI.addActionListener(new suppressAL(t));
        
        //Manipulation->Sum over one of the dimensions
        JMenuItem dimSumMI = new JMenuItem("Sum over one of the dimensions");
        manipMenu.add(dimSumMI);
        class dimSumAL implements ActionListener{
        	Topomap t;
        	String dir;
        	dimSumAL(Topomap tPrime, String dirPrime){
        		t=tPrime;
        		dir=dirPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			lv = new LayerViewer(TopomapUtil.sumOneDimension(t, Integer.parseInt(JOptionPane.showInputDialog("Enter the code:\r\n\t 0 - Energy \r\n\t 1 - X \r\n\t 2 - Y"))), dir, 1024);
    			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        	}
        }
        dimSumMI.addActionListener(new dimSumAL(t,dir));
        
        //Manipulation->Smooth each spectrum with a gaussian
        JMenuItem smoothSpecMI = new JMenuItem("Smooth each spectrum with a gaussian");
        manipMenu.add(smoothSpecMI);
        class smoothSpecAL implements ActionListener{
        	Topomap t;
        	smoothSpecAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			double[] tempSpec;
    			double L = Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels."));
    			for (int i = 0; i < t.nx; i++)
    				for (int j = 0; j < t.ny; j++)
    				{
    					tempSpec = ArrayOps.gaussSmooth(t.getSpectrum(i, j), L);
    					for (int k = 0; k < t.nlayers; k++)
    						t.data[k][i][j] = tempSpec[k];
    				}
        	}
        }
        smoothSpecMI.addActionListener(new smoothSpecAL(t));
        
        //Manipulation->Subtract a polynomial from each spectrum
        JMenuItem subPolySpecMI = new JMenuItem("Subtract a polynomial from each spectrum");
        manipMenu.add(subPolySpecMI);
        class subPolySpecAL implements ActionListener{
        	Topomap t;
        	subPolySpecAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			int degree = Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial."));
    			if (JOptionPane.showConfirmDialog(null, "Save the fits before subtracting?") == JOptionPane.YES_OPTION)
    			{
    				double[][][] fits = new double[t.nlayers][t.nx][t.ny];
    				double[] tempSpec = new double[t.nlayers];
    				double[][] tempans;
    				for (int i = 0; i < t.nx; i++)
    					for (int j = 0; j < t.ny; j++)
    					{
    						t.putSpectrum(i, j, tempSpec);
    						tempans = ArrayOps.subtractPolynomialFit(t.v, tempSpec, degree);
    						for (int k = 0; k < t.nlayers; k++)
    						{
    							t.data[k][i][j] = tempans[0][k];
    							fits[k][i][j] = tempans[1][k];
    						}
    					}
    				
    				Topomap.writeBIN(Topomap.newTopomap(t, fits), fc);
    			}
    			else
    			{
    				TopomapUtil.subtractPolynomialFitEachSpectrum(t, degree);
    			}
        	}
        }
        subPolySpecMI.addActionListener(new subPolySpecAL(t));
        
        //Manipulation->Try to filter out noise
        JMenuItem tryNoiseFiltMI = new JMenuItem("Try to filter out noise");
        manipMenu.add(tryNoiseFiltMI);
        class tryNoiseFiltAL implements ActionListener{
        	Topomap t;
        	tryNoiseFiltAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			int degree = Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial."));
    			double[][][] fits = new double[t.nlayers][t.nx][t.ny];
    			double[] tempSpec = new double[t.nlayers];
    			double[][] tempans = null;
    			for (int i = 0; i < t.nx; i++){System.out.print("" + i + " ");
    				for (int j = 0; j < t.ny; j++)
    				{
    					t.putSpectrum(i, j, tempSpec);
    					tempans = ArrayOps.subtractPolynomialFit(t.v, tempSpec, degree);
    					for (int k = 0; k < t.nlayers; k++)
    					{
    						t.data[k][i][j] = tempans[0][k];
    						fits[k][i][j] = tempans[1][k];
    					}
    				}
    			}
    			boolean doEachSpectrum = JOptionPane.showConfirmDialog(null, "To fourier-filter each spectrum, click Yes.\r\n" +
    					"To Fourier-filter the map as a whole (i.e. vertical pillars of the 3D FFT) click No, or Cancel.") == JOptionPane.YES_OPTION;
    			if (doEachSpectrum){
    				
    				if (JOptionPane.showConfirmDialog(null, "Save the fits?") == JOptionPane.YES_OPTION)
    					Topomap.writeBIN(Topomap.newTopomap(t, fits), fc);
    					
    				String input = JOptionPane.showInputDialog("Enter upper left corner of the box, spearated by commas.");
    				int imin = Integer.parseInt(input.split(",")[0].trim());
    				int jmin = Integer.parseInt(input.split(",")[1].trim());
    	//			input = JOptionPane.showInputDialog("Enter width and height of the box, separated by commas.");
    	//			int width = Integer.parseInt(input.split(",")[0].trim());
    	//			int height = Integer.parseInt(input.split(",")[1].trim());
    	//			ArrayList<TopomapUtil.FilterBox> box = new ArrayList<TopomapUtil.FilterBox>();
    	//			box.add(new TopomapUtil.FilterBox(imin, jmin, width, height));
    	//			Topomap ans = TopomapUtil.fourierFilterTheMap(t, box);
    	//			for (int i = 0; i < t.nx; i++)
    	//				for (int j = 0; j < t.ny; j++)
    	//					for (int k = 0; k < t.nlayers; k++)
    	//						t.data[k][i][j] = ans.data[k][i][j] + fits[k][i][j];
    				for (int i = 0; i < t.nx; i++)
    					for (int j = 0; j < t.ny; j++){
    						tempans[0] = FFTOps.getFourierFiltered(t.getSpectrum(i, j), new int[] {imin}, new int[] {jmin});
    						for (int k = 0; k < t.nlayers; k++)
    							t.data[k][i][j] = tempans[0][k] + fits[k][i][j];
    					}
    			}
    			else
    			{
    				int nranges = Integer.parseInt(JOptionPane.showInputDialog("How many ranges shall we filter?"));
    				boolean getFromFile = JOptionPane.showConfirmDialog(null, "Get ranges from a file?") == JOptionPane.YES_OPTION;
    				ArrayList<TopomapUtil.FilterBox> boxes = new ArrayList<TopomapUtil.FilterBox>();
    				int[] minsX = new int [nranges];
    				int[] minsY = new int [nranges];
    				int[] widths = new int [nranges];
    				int[] heights = new int [nranges];
    				if (!getFromFile){
    					for (int i = 0; i < nranges; i++)
    					{
    						String imput2 = JOptionPane.showInputDialog("Enter coordinates of the upper-left corner of the box, comma separated.");
    						minsX[i] = Integer.parseInt(imput2.split(",")[0].trim());
    						minsY[i] = Integer.parseInt(imput2.split(",")[1].trim());
    						imput2 = JOptionPane.showInputDialog("Enter the width and height of the box, comma separated.");
    						widths[i] = Integer.parseInt(imput2.split(",")[0].trim());
    						heights[i] = Integer.parseInt(imput2.split(",")[1].trim());
    						boxes.add(new TopomapUtil.FilterBox(minsX[i], minsY[i], widths[i], heights[i]));
    					}
    				
    					String[] lines = new String[nranges+2];
    					lines[0] = "The filter ranges were ";
    					lines[1] = "MinX\tMinY\tWidth\tHeight\t";
    					for (int i = 0; i < nranges; i++)
    					{
    						lines[i+2] = minsX[i] + "\t" + minsY[i] + "\t" + widths[i] + "\t" + heights[i];
    					}
    					JOptionPane.showMessageDialog(null, "Save the record.");
    					FileOps.writeLines(fc, lines);
    				}
    				else
    				{
    					double[][] table = ColumnIO.readNColumns(FileOps.selectOpen(fc), 4, 2);
    					for (int i = 0; i < nranges; i++)
    					{
    						minsX[i] = (int)table[0][i];
    						minsY[i] = (int)table[1][i];
    						widths[i] = (int)table[2][i];
    						heights[i] = (int)table[3][i];
    						boxes.add(new TopomapUtil.FilterBox(minsX[i], minsY[i], widths[i], heights[i]));
    						
    					}
    				}
    				
    				Topomap temp = TopomapUtil.fourierFilterTheMap(t, boxes);
    				FieldOps.copy(temp.data, t.data);
    				FieldOps.plusEquals(t.data, fits);
    			}
        	}
        }
        tryNoiseFiltMI.addActionListener(new tryNoiseFiltAL(t));
        
        //Manipulation->Do 3D smoothing
        JMenuItem threeDSmoothMI = new JMenuItem("Do 3D smoothing");
        manipMenu.add(threeDSmoothMI);
        class threeDSmoothAL implements ActionListener{
        	Topomap t;
        	threeDSmoothAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
				double li = Double.parseDouble(JOptionPane.showInputDialog("Enter the real-space smoothing length in pixels."));
				double lk = Double.parseDouble(JOptionPane.showInputDialog("Enter the energy smoothing length in pixels."));
				boolean savenow = JOptionPane.showConfirmDialog(null, "Save immediately?") == JOptionPane.YES_OPTION;
				File save = null;
				if (savenow) save = FileOps.selectSave(fc);
				
				FieldOps.copy(FieldOps.getGaussianSmoothing3D(t.data, lk, li, li), t.data);
				if (savenow) Topomap.writeBIN(t, save.toString());
        	}
        }
        threeDSmoothMI.addActionListener(new threeDSmoothAL(t));
        
        //Manipulation->Subtract another map
        JMenuItem subMapMI = new JMenuItem("Subtract another map");
        manipMenu.add(subMapMI);
        class subMapAL implements ActionListener{
        	Topomap t;
        	subMapAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
				Topomap m = Topomap.open(fc);
				FieldOps.minusEquals(t.data, m.data);
        	}
        }
        subMapMI.addActionListener(new subMapAL(t));
        
        //Manipulation->Normalize spectra to some data
        JMenuItem normalizeSpectraMI = new JMenuItem("Normalize spectra by some region");
        manipMenu.add(normalizeSpectraMI);
        class normalizeSpectraAL implements ActionListener{
        	Topomap t;
        	normalizeSpectraAL(Topomap tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		String message = "Enter the voltages to integrate between (inclusive), separated by a comma.\nThe data range is " + t.v[0] + " to " + t.v[t.nlayers-1];
				String vLimits = JOptionPane.showInputDialog(message);
				String[] tokens = vLimits.split(",");
				for (int i = 0; i < tokens.length; i++)
					tokens[i] = tokens[i].trim();
				double first = Double.parseDouble(tokens[0]);
				double second = Double.parseDouble(tokens[1]);
				double[][] total = new double[t.x.length][t.y.length];
				FieldOps.zero(total);
				int flag=0;
				for(int x=0;x<t.nlayers;x++){
					if((first<=t.v[x] && t.v[x]<=second) || (second<=t.v[x] && t.v[x]<=first)){
						total = FieldOps.add(total, t.data[x]);
						flag=x;
					}
				}
				for(int a=0;a<t.nlayers;a++){
					for(int b=0;b<t.x.length;b++){
						for(int c=0;c<t.y.length;c++){
							t.data[a][b][c] = t.data[a][b][c]/total[b][c] * t.mean[flag];
						}
					}
				}
        	}
        }
        normalizeSpectraMI.addActionListener(new normalizeSpectraAL(t));
        
        
        
        
        setJMenuBar(menuBar);
		
		snpts = t.nlayers-1;
		this.t = t;
		this.dir = dir;
		this.defaultSize = size;
		
		if (fc != null && fc.getSelectedFile() != null) name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length());
		else if (dir == null) dir = Topomap.stddir;
		fc = new JFileChooser(dir);
		
		spectrum = new double [t.nlayers];
		t.putSpectrum(0, 0, spectrum);
		spec = new GraphDrawerCart("Spectrum", t.v, spectrum);
		twoV = new double[] {t.v[0], t.v[0]};
		double min = ArrayOps.min(t.data);
		double max = ArrayOps.max(t.data);
		this.realmin = new double [t.nlayers];
		this.realmax = new double [t.nlayers];
		realdelta = new double[t.nlayers];
		this.fftmin = new double [t.nlayers];
		this.fftmax = new double [t.nlayers];
		fftdelta = new double[t.nlayers];

		this.min = realmin;
		this.max = realmax;
		twoVy = new double[] {min, max};
		recalculateBounds();
		delta = realdelta;
		csh = new ColorScaleHolder(min, max);
		spec.setXY(new GraphDrawerCart.GraphObject(twoV, twoVy), false, false, 1);
		spec.setRange(min, max);
		spec.showWindow();
		
		N = t.nx;
		sx = N;
		drawField = t.data[0];
		if (defaultSize == 512)
		{
			WIDTH = 1000;
			HEIGHT = 562;
		}
		//System.out.println(drawField.length);
		while (drawField.length > defaultSize && drawField.length/2 >= defaultSize)
		{
			//System.out.println(drawField.length);
			this.drawField = FieldOps.reduce(2, this.drawField);
			sx = this.drawField.length;
			zoomLevel++;
			zoomFactor *= 2;
			imageOverField /= 2;
		}
		//Controls the scaling of the pixels
		while (sizeratio*sx < defaultSize)
		{
			sizeratio++;
			imageOverField++;
			//System.out.println(sizeratio);
			//System.out.println("sx"+sx);
		}
		
		writepoint[0] = ox + sizeratio*sx + 50;
		WIDTH = sizeratio*t.nx + 450;
		HEIGHT = sizeratio*t.ny + 200;
		
		sy = this.drawField[0].length;
		setFieldInfo();
		image = new BufferedImage(sx*sizeratio, sy*sizeratio, BufferedImage.TYPE_INT_RGB);
//		resetColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data), true);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		
		setTitle("Topomap");
		showWindow();

		s = new SliderPanel(this, new JFrame());
		s.show();
		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		addMouseListener(this);
	}
	
	public void setFieldInfo()
	{
//		double[] bounds = calc.getMinMax(real);
		double[] bounds = null;
		bounds = new double[] {ArrayOps.min(drawField), ArrayOps.max(drawField)};
			
		fmax = bounds[1]; fmin = bounds[0];
		fdelta = fmax - fmin;
		setTitle("Range: [" + fmin + ", " + fmax + "]" + "     " + "p = " + para);
	}
	public void resetCSHFromSlider(double downnum, double upnum)
	{
		scalemax = min[para] + delta[para]*upnum;
		scalemin = min[para] + delta[para]*downnum;
		
		setTitle("Range: [" + scalemin + ", " + scalemax + "]" + "     " + "p = " + para);
		csh.reBoundColorScale(scalemin, scalemax);
	}
//	public void resetColorScale()
//	{
//		if (real) scale = ColorScales.getNew(scalemin, scalemax, currentCScale);
//		else cscale = new ColorScales.MYC2d(scalemax, scalemin, 2*Math.PI);
//	}

	public void formImage()
	{
		SRAW.writeImage(image, drawField, csh.getScale(), sizeratio);
	}
	public void showWindow()
	{
		setMaximumSize(new Dimension(WIDTH, HEIGHT));
		setSize(WIDTH, HEIGHT);
		repaint();
		setVisible(true);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
	public void paint(Graphics g)
	{
		if (refresh){
			g.clearRect(100, 100, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			if (lc != null || sct != null) drawLine(g, csh.getUnusedColor());
			if (imps != null && !showFFT)
			{
				drawImpMarks(g, csh.getUnusedColor());
			}
			if (fitter != null)
			{
				double[] rOn = fitter.onLatt.getPixelCoords(new double[] {N/2, 0});
				double[] rOff = fitter.offLatt.getPixelCoords(new double[] {N/2, 0});
				g.setColor(csh.getUnusedColor());
				g.drawLine(screenX(N/2), screenY(N/2), screenX(rOn[0]), screenY(rOn[1]));
				g.setColor(Color.BLACK);
				g.drawLine(screenX(N/2), screenY(N/2), screenX(rOff[0]), screenY(rOff[1]));
				Printer.printlnHorizontal(rOn);
			
			}
			
			refresh = false;
			
		}
		g.setColor(Color.BLACK);
		drawText(g);
		menuBar.paintImmediately(0,0,WIDTH,50);
	}
	public void drawText(Graphics g)
	{
		g.clearRect(writepoint[0], writepoint[1]-linesize, 500, 800);
		String it = "";
		it += "Current mouse position in map: (" + currenti + ", " + currentj + ")" + "\r\n";
		it += "With respect to the center, that is (" + (currenti-t.nx/2) + ", " + (currentj-t.ny/2) + ") \r\n";
		if (!showFFT)
			it += "In metric units that is " + Printer.vectorPFormat(t.getMetricCoords(currenti, currentj)) + "\r\n"; 
		else
			it += "In metric units that is " + Printer.vectorPFormat(t.getFourierMetricCoords(currenti, currentj)) + "\r\n"; 
			
		it += "This vector has magnitude " + NumFormat.scientific(currentr,3) + " \r\n";
		if (!showFFT)
			it += "(metric " + NumFormat.scientific(Distance.distance(t.getMetricCoords(t.nx/2, t.ny/2), t.getMetricCoords(currenti, currentj)),3) + ")\r\n";
		else
			it += "(metric " + NumFormat.scientific(Complex.mag(t.getFourierMetricCoords(currenti, currentj)),3) + ")\r\n";
			
		it += "And corresponding angle " + NumFormat.scientific(Math.toDegrees(FieldOps.atan((currenti-t.nx/2), (currentj-t.ny/2))),3) + "degrees. \r\n";
		it += "The value of the drawn field there is " + NumFormat.scientific(drawField[currenti][currentj],3) + ".\r\n";
		if (lc != null)
			it += "The line is from " + Printer.vectorPFormat(line[0]) + " to " + Printer.vectorPFormat(line[1]) + "\r\n";
		String[] lines = it.split("\r\n");
		for (int i = 0; i < lines.length; i++)
			g.drawString(lines[i], writepoint[0], writepoint[1] + i*linesize);
	}
	public void drawImpMarks(Graphics g, Color c)
	{
		double[] r = new double[2];
		g.setColor(c);
		for (int i = 0; i < imps.size(); i++)
		{
			r = imps.get(i).pixelPos;
			drawPlus(g, (int)(r[0]*sizeratio)+ox, (int)(r[1]*sizeratio)+oy, 5);
		}
	}
	public void drawLine(Graphics g, Color c)
	{
		g.setColor(c);
		g.drawString("A", (int)(line[0][0]*sizeratio)+ox, (int)(line[0][1]*sizeratio)+oy);
		g.drawString("B", (int)(line[1][0]*sizeratio)+ox, (int)(line[1][1]*sizeratio)+oy);
		g.drawLine((int)(line[0][0]*sizeratio)+ox, (int)(line[0][1]*sizeratio)+oy, (int)(line[1][0]*sizeratio)+ox, (int)(line[1][1]*sizeratio)+oy);
	}
	public void drawPlus(Graphics g, int x, int y, int r)
	{
		g.drawLine(x-r, y, x+r, y);
		g.drawLine(x, y-r, x, y+r);
	}

	public void putFFTSpectrum()
	{
		for (int i = 0; i < t.nlayers; i++)
		{
			spectrum[i] = fftmag[i][currenti][currentj];
		}
	}
	public void moveLine(int dx, int dy)
	{
		line[nearestEndIndex][0] += dx/(double)sizeratio;
		line[nearestEndIndex][1] += dy/(double)sizeratio;
	}
	public void moveStripEnd(int dx, int dy)
	{
	
		
	}
	public void resetGraphics(boolean redoColorScale, double downNum, double upNum)
	{
		if (!showFFT) drawField = t.data[para];
		else drawField = fftmag[para];
		imageOverField = 1;
		//seems like zoomLevel, zoomFactor, etc. not used elsewhere
		while (drawField.length > defaultSize && drawField.length/2 >= defaultSize)
		{
			this.drawField = FieldOps.reduce(2, this.drawField);
			sx = this.drawField.length;
			zoomLevel++;
			zoomFactor *= 2;
			imageOverField /= 2;
		}
		//Changed to accept size ratios other than 2^n
		while (sizeratio*sx < defaultSize)
		{
			sizeratio++;
			imageOverField++;
		}
//		setFieldInfo();
		if (redoColorScale)
			csh.resetColorScaleChange(min[para] + delta[para]*downNum, min[para] + delta[para]*upNum, false);
	    this.formImage();
	    this.twoVy[0] = csh.getScale().getMin();
	    this.twoVy[1] = csh.getScale().getMax();
	    refresh = true;
	    repaint();
	}
	
	
	public void keyPressed(KeyEvent arg0) {
	}
	public void keyReleased(KeyEvent arg0) {
	}
	public void keyTyped(KeyEvent arg0) {
		System.out.println(arg0.getKeyChar());
		if (arg0.getKeyChar() == 'f')
		{
			//take Fourier transform.
			if (fftmag == null)
			{
				refreshFFT = false;
				fftmag = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
				{
					setTitle("Doing FFT: " + i);
					FFTOps.obtainFFTmagCent(t.data[i], fftmag[i]);
					if (ftlog)
						FieldOps.log(fftmag[i]);
				}
				recalculateFFTBounds();
				setTitle("Done");
			}
			else if (refreshFFT)
			{
				refreshFFT = false;
				for (int i = 0; i < t.nlayers; i++)
				{
					setTitle("Doing FFT: " + i);
					FFTOps.obtainFFTmagCent(t.data[i], fftmag[i]);
					if (ftlog)
						FieldOps.log(fftmag[i]);
				}
				recalculateFFTBounds();
			}
			showFFT = !showFFT;
			if (showFFT) {
				min = fftmin;
				max = fftmax;
				delta = fftdelta;
				spec.setRange(ArrayOps.min(fftmag), ArrayOps.max(fftmag));
			}
			else{
				min = realmin;
				max = realmax;
				delta = realdelta;
				spec.setRange(ArrayOps.min(t.data), ArrayOps.max(t.data));
			}
			
			if (lc != null)
			{
				if (showFFT)
					lc.setDataset(fftmag);
				else
					lc.setDataset(t.data);
			}
			else if (sct != null)
			{
				if (showFFT)
					sct.setData(fftmag);
				else
					sct.setData(t.data);
			}
			resetGraphics(true, 0, 1);
		}
		if (arg0.getKeyChar() == 'F') //take complex FFT and spawn a new viewer
		{
			double[][][][] results = new double [2][t.nlayers][t.nx][t.ny];
			FFT2DSmall fft;
			double mag, magmin = Double.MAX_VALUE, newmag, logmag;
			for (int k = 0; k < t.nlayers; k++)
			{
				fft = FFTOps.obtainFFT(t.data[k]);
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
					{
						results[0][k][i][j] = fft.getFHat2Re(i, j);
						results[1][k][i][j] = fft.getFHat2Im(i, j);
						if (ftlog){
							mag = Complex.mag(results[0][k][i][j], results[1][k][i][j]);
							magmin = Math.min(mag, magmin);
						}
					}
			}
			if (ftlog)
				for (int k = 0; k < t.nlayers; k++)
					for (int i = 0; i < t.nx; i++)
						for (int j = 0; j < t.ny; j++)
						{
							results[0][k][i][j] *= 1/magmin; //divide out magmin so the minimum magnitude is one.
							results[1][k][i][j] *= 1/magmin;
							newmag = Complex.mag(results[0][k][i][j], results[1][k][i][j]);
							logmag = Math.log(newmag);
							results[0][k][i][j] *= logmag/newmag; //divide out magmin so the minimum magnitude is one.
							results[1][k][i][j] *= logmag/newmag;
						}
		
			new TopomapViewer_complex2(Topomap.newTopomap(t, results[0]), Topomap.newTopomap(t, results[1]),dir, 512);
		}
		if (arg0.getKeyChar() == 'g')
		{
			ftlog = !ftlog;
			refreshFFT = true;
		}
		if (arg0.getKeyChar() == 'H')
		{
			ImageEditing.copyToClipboard((BufferedImage)spec.dbimage);
		}
		if (arg0.getKeyChar() == 'd')
			doAnything();
		if (arg0.getKeyChar() == ' '){
//			if (this.showFFT)
//				csh.reBoundColorScale(ArrayOps.min(fftmag), ArrayOps.max(fftmag));
//			else
//				csh.reBoundColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data));
//			resetGraphics(true, 0, 1);
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'p') //Put on the clipboard a custom image:
		{
			ColorScale1D scale = null;
			String[] tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
			scale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
			int[] r0 = Printer.getTwoInts("Upper left corner?");
			int[] dr = Printer.getTwoInts("Width and Height?");
			BufferedImage im = ImageEditing.getBufferedImage(showFFT ? fftmag[para] : t.data[para], scale);
			im = ImageEditing.getSubset(im, r0[0], dr[0], r0[1], dr[1]);
			int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
			if (resizeFactor > 1) im = ImageEditing.getEnlarged(im, resizeFactor);
			ImageEditing.copyToClipboard(im);
//			File f = FileOps.selectSave(fc);
//			if (f != null)
//				SRAW.writeImage(f.toString(), image);
		}
		if (arg0.getKeyChar() == 'G')
		{
				int typeOfGIF = Integer.parseInt(JOptionPane.showInputDialog("Select the output type:\r\n" +
						" 0 - GIF with topo and FFT inset\r\n" +
						" 1 - GIF with FFT only (and voltage)\r\n" + 
						" 2 - Large BMP with layers as tiles\r\n" + 
						" 3 - Movie from BMP to AVI\r\n" + 
						" 4 - Tiled BMP with movie settings\r\n" +
						" 5 - Movie from BMP to AVI with both\r\n" +
						" 6 - Movie from BMP to AVI with another map\r\n" +
						" 7 - Movie from BMP to AVI about important FFT point\r\n" + 
						" 8 - Tiled BMP, movie settings, important FFT point\r\n"));
				File f;
				switch(typeOfGIF){
					case 0:
						f = FileOps.selectSave(fc);
						Layer topo = Layer.openFree(fc);
						RHKFileOps.doFitting(topo, RHKFileOps.getUserFittingChoice());
						if (showFFT)
							ExportUtil.exportGIFForPPTSlide(t, topo, csh.getScale(), null, false, true, true, 500, f.toString(), csh.getCurrentCS());
						else
							ExportUtil.exportGIFForPPTSlide(t, topo, null, csh.getScale(), false, true, true, 500, f.toString(), csh.getCurrentCS());
						break;
					case 1:
						f = FileOps.selectSave(fc);
						TopomapUtil.writeGIFWithVoltage(t, f.toString(), Integer.parseInt(JOptionPane.showInputDialog("How much to blow it up by?")), csh.getCurrentCS());
						break;
					case 2:
						f = FileOps.selectSave(fc);
						int nPerLine = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
//						BufferedImage image = ImageEditing.createTileImage(showFFT ? fftmag : t.data, csh.getScale(), nPerLine, t.v);
						BufferedImage image = ImageEditing.createTileImage(showFFT ? fftmag : t.data, csh.getScale(), nPerLine, t.v, csh.getCurrentCS());
						SRAW.writeImage(f.toString(), image);
						ImageEditing.copyToClipboard(image);
						break;
					case 3:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						boolean useCurrentScale = JOptionPane.showConfirmDialog(null, "Lock the color scale?") == JOptionPane.YES_OPTION; 
						double minperc = useCurrentScale ? 0 : Double.parseDouble(JOptionPane.showInputDialog("Minimum distribution cutoff?"));
						double maxperc = useCurrentScale ? 1 : Double.parseDouble(JOptionPane.showInputDialog("Maximum distribution cutoff?"));
						ColorScale1D scale = null;
						if (useCurrentScale) {
							String[] tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							scale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						int[] rF = new int [2];
						if (t.nx == t.ny) {
							int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
							rF[0] = resizeFactor; rF[1] = resizeFactor;
						}
						else {
							int resizeFactorX = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor? X", "" + 1));
							int resizeFactorY = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor? Y", "" + 1));
							rF[0] = resizeFactorX; rF[1] = resizeFactorY;
					
						}
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						BufferedImage[] stack;
						if (t.nx == t.ny)
							stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, rF[0], useCurrentScale ? scale : null);
						else
							stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, rF, useCurrentScale ? scale : null);
							
						BufferedImage[] substack = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++){
							substack[i-imin] = stack[i];
						}
						if (imps != null)
							for (int i = 0; i < substack.length; i++)
							{
								Graphics g = substack[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).pixelPos;
									drawPlus(g, (int)(r[0]*rF[0]), (int)(r[1]*rF[1]), 5);
								}
							}
						
						for (int i = 0; i < substack.length; i++)
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), substack[i]);
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substack.length-1, 5));
						break;
					}
					case 4:{
						String shortname = name.contains(" ") ? name.split(" ")[0] : name;
						int nPerLine1 = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
						boolean useCurrentScale = JOptionPane.showConfirmDialog(null, "Lock the color scale?") == JOptionPane.YES_OPTION; 
						double minperc = useCurrentScale ? 0 : Double.parseDouble(JOptionPane.showInputDialog("Minimum distribution cutoff?"));
						double maxperc = useCurrentScale ? 1 : Double.parseDouble(JOptionPane.showInputDialog("Maximum distribution cutoff?"));
						ColorScale1D scale = null;
						if (useCurrentScale) {
							String[] tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							scale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						
						BufferedImage[] stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, resizeFactor, useCurrentScale ? scale : null);
						BufferedImage[] substack = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++)
							substack[i-imin] = stack[i];
						if (imps != null)
							for (int i = 0; i < substack.length; i++)
							{
								Graphics g = substack[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).pixelPos;
									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
								}
							}
						BufferedImage tiles = ImageEditing.createTileImage(substack, nPerLine1, Color.WHITE);
						ImageEditing.copyToClipboard(tiles);
						SRAW.writeImage(FileOps.selectSave(fc).toString(), tiles);
						
						break;
					}
					case 5:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						
						String firstInput = JOptionPane.showInputDialog("Enter the half-length of the box you want to capture.", "" + t.nx/2);
						int width = Integer.parseInt(firstInput)*2;
						
						double minpercR = 0;
						double maxpercR = 1;
						double minpercK = 0;
						double maxpercK = 1;
						boolean useFixedScaleReal = JOptionPane.showConfirmDialog(null, "Use fixed color scale in real space?") == JOptionPane.YES_OPTION;
						ColorScale1D realScale = null;
						if (useFixedScaleReal){
							String[] tok;
							if (showFFT) tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.",  "" + ArrayOps.min(t.data) + "," + ArrayOps.max(t.data)).split(",");
							else tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							
							realScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						else{
						minpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Minimum distribution cutoff?"));
						maxpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Maximum distribution cutoff?"));
						}

						ColorScale1D kScale = null;
						boolean useFixedScaleK = JOptionPane.showConfirmDialog(null, "Use fixed color scale in k space?") == JOptionPane.YES_OPTION;
						if (useFixedScaleK){
							String[] tok;
							if (!showFFT) tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + ArrayOps.min(fftmag) + "," + ArrayOps.max(fftmag)).split(",");
							else tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							
							kScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						else{
						minpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Minimum distribution cutoff?"));
						maxpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Maximum distribution cutoff?"));
						}
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
//						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						
						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						BufferedImage[] stackR = ImageEditing.createImageArray(t.data, csh.getCurrentCS(), t.v, label, minpercR, maxpercR, resizeFactor, realScale);
						BufferedImage[] stackK = ImageEditing.createImageArray(fftmag, csh.getCurrentCS(), t.v, 0, minpercK, maxpercK, resizeFactor, kScale);
						BufferedImage[] substackR = new BufferedImage[imax-imin];
						BufferedImage[] substackK = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++){
							substackR[i-imin] = stackR[i];
							substackK[i-imin] = stackK[i];
							if (width < t.nx/2){
								substackR[i-imin] = ImageEditing.getSubset(substackR[i-imin], (t.nx-width)/2, width, (t.ny-width)/2, width);
								substackK[i-imin] = ImageEditing.getSubset(substackK[i-imin], (t.nx-width)/2, width, (t.ny-width)/2, width);;
							}
						}
						if (imps != null)
							for (int i = 0; i < substackR.length; i++)
							{
								Graphics g = substackK[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).pixelPos;
									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
								}
							}
						
						BufferedImage out;
						for (int i = 0; i < substackR.length; i++){
							if (showFFT)
								out = ImageEditing.weldHorizontal(substackK[i], substackR[i]);
							else
								out = ImageEditing.weldHorizontal(substackR[i], substackK[i]);
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), out);
						}
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substackR.length-1, 5));
						break;
					}
					case 6:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						
						double minpercR = 0;
						double maxpercR = 1;
						double minpercK = 0;
						double maxpercK = 1;
						Topomap other = Topomap.open(fc);
						boolean useFixedScaleReal = JOptionPane.showConfirmDialog(null, "Use fixed color scale in real space?") == JOptionPane.YES_OPTION;
						ColorScale1D realScale = null;
						if (useFixedScaleReal){
							String[] tok;
							tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
							
							realScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
						}
						else{
						minpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Minimum distribution cutoff?"));
						maxpercR = Double.parseDouble(JOptionPane.showInputDialog("Real space Maximum distribution cutoff?"));
						}

						ColorScale1D kScale = null;
						boolean useFixedScaleK = JOptionPane.showConfirmDialog(null, "Use fixed color scale in other space?") == JOptionPane.YES_OPTION;
						int cschoice = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice of color scale index.", "" + csh.getCurrentCS()));
						if (useFixedScaleK){
							String[] tok;
							tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + ArrayOps.min(other.data) + "," + ArrayOps.max(other.data)).split(",");
							
							kScale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), cschoice);
						}
						else{
						minpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Minimum distribution cutoff?"));
						maxpercK = Double.parseDouble(JOptionPane.showInputDialog("K-space Maximum distribution cutoff?"));
						}
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());
						
//						int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
						
//						BufferedImage[] stackR = ImageEditing.createImageArray(t.data, csh.getCurrentCS(), t.v, 2, minpercR, maxpercR, resizeFactor, realScale);
//						BufferedImage[] stackK = ImageEditing.createImageArray(other.data, cschoice, t.v, 0, minpercK, maxpercK, resizeFactor, kScale);
//						BufferedImage[] substackR = new BufferedImage[imax-imin];
//						BufferedImage[] substackK = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++){
							BufferedImage imr = ImageEditing.createSingleImage(t.data, csh.getCurrentCS(), t.v, 2, minpercR, maxpercR, resizeFactor, realScale, i);
							BufferedImage imk = ImageEditing.createSingleImage(other.data, cschoice, t.v, 0, minpercK, maxpercK, resizeFactor, kScale, i);
//							if (imps != null)
//							{
//								Graphics g = imr.getGraphics();
//								g.setColor(csh.getUnusedColor());
//								double[] r;
//								for (int j = 0; j < imps.size(); j++)
//								{
//									r = imps.get(j).pixelPos;
//									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
//								}
//							}
						
							BufferedImage out = ImageEditing.weldHorizontal(imr, imk);
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), out);
							
						}
						
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, imax-imin-1, 5));
						break;
					}
					case 7:{
						String shortname = name != null ? (name.contains(" ") ? name.split(" ")[0] : name) : "mfv";
						BufferedImage[] substack = getImageArrayHighSymPt();
						for (int i = 0; i < substack.length; i++)
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), substack[i]);
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substack.length-1, 5));
						break;
					}
					case 8:
					{
						int nPerLine1 = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
						BufferedImage[] substack = getImageArrayHighSymPt();
						BufferedImage tiles = ImageEditing.createTileImage(substack, nPerLine1, Color.WHITE);
						ImageEditing.copyToClipboard(tiles);
						SRAW.writeImage(FileOps.selectSave(fc).toString(), tiles);
						break;
					}
				}	
		}
		if (arg0.getKeyChar() == 'b'){
			double[] avg = t.getAverageSpectrum();
			System.out.println("Average spectrum:");
			Printer.printlnVertical(avg);
		}
		if (arg0.getKeyChar() == 'L') //put a layer onto the spectra windo
			//Actually, load an atomicCoordinates Set.
		{
//			if (!spec.drawingLayer){
//				Layer layer = Layer.open(fc);
//				ColorScale1D scale = ColorScales.getNew(layer.data, csh.getCurrentCS());
//				spec.drawingLayer = true;
//				spec.layer = layer;
//				spec.cscale = scale;
//				spec.layerImage = ImageEditing.getBufferedImage(layer.data, scale);
//				spec.plot[0].c = ColorScales.getUnusedColor(scale);
//				spec.repaint();
//			}
//			else
//			{
//				spec.drawingLayer = false;
//			}
			
		}
		if (arg0.getKeyChar() == 'c')
		{
			if (sct == null)
			{
				if (line == null) line = new double[][] {{t.nx/2, t.ny/2}, {t.nx, t.ny}};
				sct = new StripCutTool(t, 2, 0, t.nx/2, 1);
				sct.putLine(line);
				refresh = true;
				repaint();

			}
		}
		if (arg0.getKeyChar() == 'C') //copy the current image to clipboard
		{
//			ImageEditing.copyToClipboard(image);
			BufferedImage export = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
			Graphics g = export.getGraphics();
			refresh = true;
			paint(g);
			ImageEditing.copyToClipboard(ImageEditing.getSubset(export, ox, sx*sizeratio, oy, sy*sizeratio));
		}
		if (arg0.getKeyChar() == 'V')
		{
			ImageEditing.copyToClipboard(csh.getScale().getScaleImage(8));
		}
		if (arg0.getKeyChar() == 'v')
			if (JOptionPane.showConfirmDialog(null, "Save the Fourier transform?") == JOptionPane.YES_OPTION)
			{
				boolean log = JOptionPane.showConfirmDialog(null, "Log scale?") == JOptionPane.YES_OPTION;
				double[][][] fftmagtemp = new double[t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
				{
					setTitle("Doing FFT: " + i);
					FFTOps.obtainFFTmagCent(t.data[i], fftmagtemp[i]);
					if (log)
						FieldOps.log(fftmag[i]);
				}
				double[] kx = new double[t.nx];
				double[] ky = new double[t.ny];
				for (int i = 0; i < t.nx; i++)
					kx[i] = t.getFourierMetricCoords(i, 0)[0];
				for (int i = 0; i < t.ny; i++)
					ky[i] = t.getFourierMetricCoords(0, i)[1];
				Topomap.writeBIN(new Topomap(fftmagtemp, t.v, kx, ky, null), fc);
				
	
			}
		if (arg0.getKeyChar() == 'y')
		{
			//Write a spectrum of the real and imaginary parts of the FFT as a function of energy, and copy it to the cliboard.
			double[][][] fftz = new double [t.nlayers][t.nx][t.ny];
			int i = currenti;
			int j = currentj;
			String[] lines = new String[t.nlayers];
			for (int k = 0; k < t.nlayers; k++){
				FFTOps.putFFT(t.data[k], fftz, true);
				lines[k] = "" + t.v[k] + "\t"  + Complex.mag(fftz[i][j]) + "\t" + fftz[i][j][0] + "\t" + fftz[i][j][1] + "\t" + Complex.phase(fftz[i][j]);
				System.out.println(lines[k]);
//				System.out.println(i + "\t" + j + "\t" + currenti + "\t" + currentj);
			}
			Printer.copyToClipboard(lines);
		}
		if (arg0.getKeyChar() == 'q')
		{
			csh.incrementScaleIndex();
			formImage();
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'e')
		{
			csh.decrementScaleIndex();
			formImage();
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'E') //export csv
		{
			String path = FileOps.selectSave(fc).toString();
			for (int i = 0; i < t.nlayers; i++)
				ColumnIO.writeTableCSV(t.data[i], path + "_" + i);
		}

		if (arg0.getKeyChar() == 'r')
		{
			changeScaleWithLayer = !changeScaleWithLayer;
		}
		if (arg0.getKeyChar() == 'R'){
			String menu = "What would you like to do?\r\n" + 
					"0 - smooth each line with a Gaussian\r\n" + 
					"1 - replace each line with -1* its second derivative\r\n" + 
					"2 - Replace each line with its derivative\r\n" + 
					"3 - adaptively fit some lines to parabolas and print the results\r\n";
			int choice = Printer.getAnInt(menu);
			double parameter = 0;
			if (choice == 0) parameter = Printer.getADouble("Smoothing length?");
//			backup = Layer.newLayer(t, FieldOps.copy(t.data));
			for (int k = 0; k< t.nlayers; k++){
				PointSpectra ps = t.getLayer(k).toSpectra();
				if (choice == 0){
					ps = SpectraUtil.getGaussSmoothed(ps, parameter);
					t.data[k] = FieldOps.transpose(ps.data);
					resetGraphics(true, 0, 1);
				}
				if (choice == 1){
					for (int i = 0; i < ps.nspec; i++)
						ps.data[i] = ArrayOps.getDerivative(ArrayOps.getDerivative(ps.data[i]));
					FieldOps.negate(ps.data);
					t.data[k] = FieldOps.transpose(ps.data);
					resetGraphics(true, 0, 1);
				}
				if (choice == 2){
					for (int i = 0; i < ps.nspec; i++)
						ps.data[i] = ArrayOps.getDerivative(ps.data[i]);
	//				FieldOps.negate(ps.data);
					t.data[k] = FieldOps.transpose(ps.data);
					resetGraphics(true, 0, 1);
				}
			}
		}
		if (arg0.getKeyChar() == 'x'){

			int pixcoord = JOptionPane.showConfirmDialog(null, " YES: Use original length coordinates in layer \r\n NO: Use pixel length coordinates in layer", "Length Units", JOptionPane.YES_NO_OPTION);
			if (pixcoord == JOptionPane.YES_OPTION)
				Layer.writeBIN(t.getLayer(para));
			else Layer.writeBIN(t.getLayerPixels(para));
		}
		if (arg0.getKeyChar() == 's')
		{
			String[] plines = new String [t.nlayers];
			for (int i = 0; i < t.nlayers; i++)
			{
				plines[i] = t.v[i] + "\t" + spectrum[i];
				System.out.println(plines[i]);
				Printer.copyToClipboard(plines);
			}
		}
		if (arg0.getKeyChar() == 'm') //mark a point
		{
			double[] coords = {currenti, currentj};
			PointImp mark = new PointImp(coords);
			if (imps == null && lc == null) {
				imps = new ArrayList<PointImp>();
				imps.add(mark);
			}
			else if(lc==null)
				imps.add(mark);
			
			refresh=true;
			repaint();
		}
		if (arg0.getKeyChar() == 'S' && imps == null)
		{
			Topomap.writeBIN(t, fc);
		}
		if (arg0.getKeyChar() == 'S' && imps != null)
		{	
			PointSpectra spec = ImpurityUtil.getSpectraAt(ImpurityUtil.getFromList(imps), t, false, 1);
			PointSpectra.writeBIN(spec, fc);
			String[] lines = spec.toLines();
			for (int i = 0; i < lines.length; i++)
				System.out.println(lines[i]);
		}
		if (arg0.getKeyChar() == 'M')
		{
			String[] options = {"6-fold","4-fold","2-fold"};
			int fold = JOptionPane.showOptionDialog(null,"Choose symmetry. You may need lattice parameters from drift correction.","Symmetrization",JOptionPane.DEFAULT_OPTION,JOptionPane.PLAIN_MESSAGE,null,options,"");
			switch(fold){
			case 0:
				if (latt == null) latt = new AtomicCoordinatesSet(FileOps.openText(fc));
				t = TopomapUtil.makeSymmetrizedFFTs(t, latt, 4);
				break;
			case 1:
				if (latt == null) latt = new AtomicCoordinatesSet(FileOps.openText(fc));
				t = TopomapUtil.makeSymmetrizedFFTs(t, latt, 6);
				break;
			case 2:
				t = TopomapUtil.makeSymmetrizedFFTs(t, latt, 2);
				break;
			}

			Topomap.writeBIN(t, fc);
			resetGraphics(true, 0, 1);
		}
		if (arg0.getKeyChar() == 'P') //print entire stack of images with custom color scale
		{
//			File f = FileOps.selectSave(fc);
//			String name = f.toString();
//			BufferedImage im = new BufferedImage(image.getWidth(), image.getHeight(), BufferedImage.TYPE_INT_RGB);
			String firstInput = JOptionPane.showInputDialog("Enter the half-length of the box you want to capture.", "" + t.nx/2);
			int width = Integer.parseInt(firstInput)*2;
			BufferedImage im = new BufferedImage(width, width, BufferedImage.TYPE_INT_RGB);
			
			String input = JOptionPane.showInputDialog("Enter the color scale range (comma separated).", "" + csh.getScale().getMin() + "," + csh.getScale().getMax());
			double min = Double.parseDouble(input.split(",")[0].trim());
			double max = Double.parseDouble(input.split(",")[1].trim());
			ColorScale1D cs = ColorScales.getNew(min, max, csh.getCurrentCS());
			int i = para;
			//			for (int i = 0; i < t.nlayers; i++)
//			{
//				SRAW.writeImage(name + MovieMaker.fromInt(i), t.data[i], csh.getScale());  //To print the bare image:
//				if (!showFFT)
//					SRAW.writeImage(im, t.data[i], cs);
//				else
//					SRAW.writeImage(im, fftmag[i], cs);
			BufferedImage fullImage = ImageEditing.getBufferedImage(showFFT ? fftmag[i] : t.data[i], cs);
			im = ImageEditing.getSubset(fullImage, (t.nx-width)/2, width, (t.ny-width)/2, width);
			ImageEditing.copyToClipboard(im);
			//				int wx = im.getWidth() - 150;
//				int wy = im.getWidth() -20;
//				ImageEditing.writeLine(im, "V = " + t.v[i], wx, wy, csh.getUnusedColor(), Font.SERIF, Font.PLAIN, 32);
//				SRAW.writeImage(name + MovieMaker.fromInt(i), im);
//			}
		}
		if (arg0.getKeyChar() == 'I' && lc == null)
		{
			if (imps == null) {
				PointImp[] der = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
				imps = new ArrayList<PointImp>();
				for (int i = 0; i < der.length; i++)
					imps.add(der[i]);
				refresh = false; repaint();}
			else imps = null;
		}
		if (arg0.getKeyChar() == '=')
		{
			if (para < t.nlayers-1){ para++;
				twoV[0] = t.v[para];
				twoV[1] = t.v[para];
				spec.repaint();
				
				resetGraphics(changeScaleWithLayer, ((double)s.min.getValue()/1000), ((double)s.max.getValue()/1000));
				s.frame.setTitle("" + t.v[para] + "   (Layer " + para + ")");
			}
		}
		if (arg0.getKeyChar() == '-')
		{
			if (para > 0){ 
				para--;
				twoV[0] = t.v[para];
				twoV[1] = t.v[para];
				spec.repaint();
				
				resetGraphics(changeScaleWithLayer, ((double)s.min.getValue()/1000), ((double)s.max.getValue()/1000));
				s.frame.setTitle("" + t.v[para] + "   (Layer " + para + ")");
			}
		}
		if (arg0.getKeyChar() == 'N')
		{
			writeCustomPicutres();
		}
		if (arg0.getKeyChar() == 'T' && fitter == null)
		{
			if (latt == null) latt = new AtomicCoordinatesSet(FileOps.openText(fc));
			fitter = new FFTCutFitter(t.data, latt);
		}
		if (lc != null) lc.processKeyStroke(arg0.getKeyChar());
		else if (sct != null) sct.processKeyStroke(arg0.getKeyChar());		
		
		if (fitter != null){
			if (arg0.getKeyChar() == 'j') //construct the e, q point at the present mouse position and add it to the list
			{
				//first, determine whether it is on- or off- lattice.
				double[] r =new double[] {currentid-t.nx/2, currentjd-t.nx/2};
				double angon = fitter.onLatt.getAngleBetween(r, 0);
				double angoff = fitter.offLatt.getAngleBetween(r, 0);
				double dist = showFFT ? Complex.mag(t.getFourierMetricCoords(currentid, currentjd)) : Distance.distance(t.getMetricCoords(t.nx/2, t.ny/2), t.getMetricCoords(currentid, currentjd));
				boolean onit = angon < angoff;
				System.out.println("" + onit + "\t" + angon + "\t" + angoff);
				double distproj = onit ? dist*Math.cos(angon) : dist*Math.cos(angoff);
				fitter.user.addEQPair(new EQPair(t.v[para], distproj, onit));
			}
			if (arg0.getKeyChar() == 'J')
			{
				fitter.user.removeLastEQPair();
			}
			fitter.processKeyStroke(arg0.getKeyChar());
		}
		if (arg0.getKeyChar() == 'U')
		{
			if (latt == null) latt = new AtomicCoordinatesSet(ColumnIO.getASCII(FileOps.selectOpen(fc).toString()));
//			t = TopomapUtil.makeSymmetrizedFFT(t, latt, JOptionPane.showConfirmDialog(null, "Lattice is square?") == JOptionPane.YES_OPTION);
			int index = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice for the second derivative:\r\n"+
					"0 - faa\r\n"+
					"1 - fab\r\n"+
					"2 - fba\r\n"+
					"3 - fbb\r\n"+
					"4 - 1/2(fab+fba)\r\n"+
					"5 - 1/2(fab-fba) which is zero identically\r\n"+
					"6 - 1/2(faa+fbb)\r\n"+
					"7 - 1/2(faa-fbb)\r\n"));
			
			AtomicCoordinatesSet recip = latt.getReciprocalLattice();
			double[][] units = FieldOps.getGrahamSchmidtUnitVectors(recip.getA(), recip.getB());
			for (int i = 0; i < t.nlayers; i++) t.data[i] = FieldOps.getDirectional2ndDerivativeStuff(t.data[i], units[0], units[1])[index];
			this.recalculateBounds();
			this.resetCSHFromSlider(0, 1);
			this.formImage();
			this.repaint();
		}

		
	}
	/**
	 * The picture here is of the two bragg peaks in their entirety
	 */
	private void writeCustomPicutres() {
		// TODO Auto-generated method stub
		int sizeUnit = t.nx;
		int miniSizeUnit = 32;
		int b1x = 216, b1y = -218;
		int b2x =  218, b2y = 216;
		int blowUpFactor = (2*sizeUnit)/miniSizeUnit;
		ColorScale1D temp = ColorScales.getNew(0, 10000, csh.getCurrentCS());
		for (int i = 0; i < t.nlayers; i++){
			BufferedImage fft = ImageEditing.getBufferedImage(fftmag[i], temp);
			BufferedImage real = ImageEditing.getBufferedImage(t.data[i]);
			BufferedImage braggOne = ImageEditing.getSubset(fft, b1x - miniSizeUnit/2 + t.nx/2, miniSizeUnit, b1y - miniSizeUnit/2 + t.ny/2, miniSizeUnit);
			BufferedImage braggTwo = ImageEditing.getSubset(fft, b2x - miniSizeUnit/2 + t.nx/2, miniSizeUnit, b2y - miniSizeUnit/2 + t.ny/2, miniSizeUnit);
			braggOne = ImageEditing.getEnlarged(braggOne, blowUpFactor);
			braggTwo = ImageEditing.getEnlarged(braggTwo, blowUpFactor);
			BufferedImage ans = new BufferedImage(4*sizeUnit, 3*sizeUnit, BufferedImage.TYPE_INT_RGB);
			ImageEditing.setColor(ans, Color.WHITE);
			ImageEditing.copyInto(fft, ans, 0, 0);
			ImageEditing.copyInto(real, ans, sizeUnit, 0);
			ImageEditing.copyInto(braggOne, ans, 0, sizeUnit);
			ImageEditing.copyInto(braggTwo, ans, 2*sizeUnit, sizeUnit);
			SRAW.writeImage(MovieMaker.avidir + "Custom_" + i, ans);
		}
	}
	public void mouseWheelMoved(MouseWheelEvent arg0) {
	}
	public void mouseDragged(MouseEvent arg0) {
		int x = arg0.getX();
		int y = arg0.getY();
		
		double r, phi;
		double[] pos = new double [] {((x-ox)/(double)sizeratio) - t.nx/2, ((y-oy)/(double)sizeratio) - t.ny/2};
		r = Complex.mag(pos);
		phi = Complex.phase(pos);
		
		if (lc != null)
		{
			moveLine((x - currentx), (y - currenty));
			lc.putPoints(line);
			lc.refresh();
			refresh = true;
			repaint();
		}
		else if (sct != null)
		{
			
			double dr, dphi;
			dr = r - currentr;
			dphi = phi-currentphi;
			System.out.println(dr + "\t" + dphi + "\t" + nearestEndIndex);
			if (nearestEndIndex == 0 && sct.cut[0].r0 + dr > 0)
				sct.setR0(sct.cut[0].r0 + dr);
			else if (sct.cut[0].r1 + dr > 0)
				sct.setR1(sct.cut[0].r1 + dr);
			sct.setPhi(sct.cut[0].phi + dphi);
			sct.putLine(line);
			sct.refresh();
			refresh = true;
			repaint();
		}
		else
		{
			Point currCoords = arg0.getLocationOnScreen();
            this.setLocation(currCoords.x - mouseDownCompCoords.x, currCoords.y - mouseDownCompCoords.y);
		}
//		System.out.println(x);
		currentx = x; currenty = y;
		currentr = r;
		currentphi = phi;
		}
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		calcx = currentx - ox;
		calcy = currenty - oy;
		currentid = calcx/(double)sizeratio;
		currentjd = calcy/(double)sizeratio;
		currenti = (int)currentid;
		currentj = (int)currentjd;
		currentr = Distance.distance((currenti-t.nx/2), (currentj-t.ny/2));
		currentphi = FieldOps.atan(currentid-t.nx/2, currentjd-t.ny/2);
		if (currenti < 0) currenti = 0;
		if (currenti >= t.nx) currenti = t.nx-1;
		if (currentj < 0) currentj = 0;
		if (currentj >= t.ny) currentj = t.ny-1;
		refresh = false;
		repaint();
		
		spec.setTitle("Spectrum [" + currenti + ", " + currentj + "]");
		if (withinBounds(currenti, currentj, t.data[0]) && !showFFT)
			t.putSpectrum(currenti, currentj, spectrum);
		else if (showFFT){
			putFFTSpectrum();
			//			spec.setRange(ArrayOps.min(spectrum), ArrayOps.max(spectrum));
		}
		spec.repaint();
		
		if (lc != null || sct != null)
		{
			d1 = Distance.distance(currentid - line[0][0], currentjd - line[0][1]);
			d2 = Distance.distance(currentid - line[1][0], currentjd - line[1][1]);
			if (d1 < d2) nearestEndIndex = 0;
			else nearestEndIndex = 1;
//			System.out.println(nearestEndIndex);
		}		
//		refresh = false;
//		repaint();
	}
	public int screenX(double currentid)
	{
		return (int)(currentid*sizeratio + ox);
	}
	public int screenY(double currentjd)
	{
		return (int)(currentjd*sizeratio + oy);
	}
	public void mouseClicked(MouseEvent arg0) {
		
		if (fitter != null)
		{
			
		}
		
	}
	public void mouseEntered(MouseEvent arg0) {
	}
	public void mouseExited(MouseEvent arg0) {
	}
	public void mousePressed(MouseEvent arg0) {
		mouseDownCompCoords = arg0.getPoint();
	}
	public void mouseReleased(MouseEvent arg0) {
		mouseDownCompCoords = null;
	}
	
	
	public static boolean withinBounds(int i, int j, double[][] array)
	{
		return i >= 0 && j >= 0 && i < array.length && j < array[i].length;
	}
	public void recalculateBounds()
	{
		for (int i = 0; i < t.nlayers; i++)
		{
			this.realmin[i] = FieldOps.min(t.data[i]);
			this.realmax[i] = FieldOps.max(t.data[i]);
			this.realdelta[i] = this.realmax[i]-this.realmin[i];
		}
	}
	public void recalculateFFTBounds()
	{
		for (int i = 0; i < t.nlayers; i++)
		{
			this.fftmin[i] = FieldOps.min(fftmag[i]);
			this.fftmax[i] = FieldOps.max(fftmag[i]);
			this.fftdelta[i] = this.fftmax[i]-this.fftmin[i];
		}
	}
	public void doAnything()
	{
		int o = Integer.parseInt(JOptionPane.showInputDialog(null, "What would you like to do?\r\n"+
				"1 - Do a work-function fitting\r\n" +
				"2 - Symmetrize\r\n" +
				"3 - Do some data processing technique to each layer\r\n"+
				"4 - Rescale the bias voltages by a constant factor (for Z-map)\r\n"+ 
				"5 - print a spectrum of the correlation with a layer (opens a layer)\r\n"+
				"6 - Split the map into bins based on an existing layer\r\n" +
				"7 - Split the map into smoothed bins based on an existing layer\r\n" + 
				"8 - Get locally-averaged spectra at impurities\r\n" +
				"9 - Get a simple histogram of the spectra\r\n" +
				"10 - Smoothly suppress features to the else average\r\n" +
				"11 - Replace the map with its full 3D FFT\r\n" + 
				"12 - Sum over one of the dimensions\r\n" +
				"13 - Smooth each spectrum with a gaussian\r\n" + 
				"14 - Subtract a polynomial from each spectrum\r\n" +
				"15 - Graph the average FFT of all the spectra\r\n" + 
				"16 - Try to filter out the noise\r\n" + 
				"17 - Replace each spectrum by its FFT\r\n" +
				"18 - Replace each spectrum by its derivative\r\n" + 
				"19 - Get a histogram of all the values\r\n" + 
				"20 - Get a simple gap map.\r\n" +
				"21 - Use 3D smoothing\r\n" + 
				"22 - Get layer showing the sum of each spectrum\r\n" +
				"23 - Split the given layer against smoothed bins and write a map of the subsets\r\n" +
				"24 - Split the map layers against the layers of another 3D dataset (smoothed) \r\n" + 
				"25 - Replace with 3D laplacian\r\n" + 
				"26 - Subtract another map\r\n" + 
				"27 - Save a cropped version of this map\r\n" + 
				"28 - Copy to the clipboard the dI/dV at this voltage as a function of time\r\n" +
				"29 - Save the map as a PointSpectra obect\r\n" + 
				"30 - Export all layers as .bin files\r\n" +
				"31 - Do an adptive peak fitting to a parabola starting from the center pixel\r\n" +
				"32 - Get a complex map of the smoothed, shifted Fourier transform\r\n"+
				"33 - Spatially filter each layer against a map of spatial filters(open the map)\r\n"
				));
		
		
		switch(o)
		{
		case 1:
			TopomapUtil.fitEachSpectrumExponential(t, "");
			break;
		case 2:
			int o2 = JOptionPane.showConfirmDialog(null, "Yes for 6-fold\r\nNo for 4-fold");
			if (latt == null) latt = new AtomicCoordinatesSet(FileOps.openText(fc));

			t = TopomapUtil.makeSymmetrizedFFTs(t, latt, o2 == JOptionPane.NO_OPTION);
			Topomap.writeBIN(t, fc);
			resetGraphics(true, 0, 1);
			break;
		case 3:
			int o3 = RHKFileOps.getUserFittingChoice();
			RHKFileOps.doFittingTopo(t, o3);
//			recalculateBounds();
			break;
		case 4:
			double fz = Double.parseDouble(JOptionPane.showInputDialog("Enter the factor to multiply Z by."));
			for (int i = 0; i < t.nlayers; i++)
				t.v[i] *= fz;
		break;
		case 5:
			if (JOptionPane.showConfirmDialog(null, "Open a map instead?") == JOptionPane.YES_OPTION){
				Topomap ls = Topomap.open(fc);
				String[] lines = new String [t.nlayers];
				for (int i = 0; i < t.nlayers; i++){
					lines[i] = "" + t.v[i];
					for (int j = 0; j < ls.nlayers; j++)
					{
						lines[i] += "\t" + FieldOps.correlation(ls.data[j], t.data[i]);
					}
					System.out.println(lines[i]);
				}
			}
			else{
				Layer l = Layer.open(fc);
				System.out.println();
				for (int i = 0; i < t.nlayers; i++)
					System.out.println("" + t.v[i] + "\t" + FieldOps.correlation(l.data, t.data[i]));
			}
			break;
		case 6:
			boolean bool = JOptionPane.showConfirmDialog(null, "Apply a boolean filter first?") == JOptionPane.YES_OPTION;
			if (!bool){
				Topomap binshist = TopomapUtil.splitTopomapByBinsOfLayerFast_autoHist(t, dir, name, Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins")), fc);
//				Topomap.writeBIN(binshist, fc);
				if (t.nlayers > 1){
				TopomapViewer tv = new TopomapViewer(binshist, dir, getGoodSize(binshist.nx));
				tv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
				}
			}
			else
			{
				boolean[][] thebools = ImageEditing.loadFromImage(FileOps.selectOpen(fc).toString(), Color.WHITE);
				int nb = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
				Layer binstuff = TopomapUtil.splitTopomapByBinsOfLayer(t, dir, name, thebools, nb, fc);
				ColorScale1D binStuff = ColorScales.getNew(0, nb-1, csh.getCurrentCS());
				BufferedImage binBool = ImageEditing.getBufferedImage(binstuff.data, binStuff);
				FieldOps.negate(thebools);
				ImageEditing.setTrueToColor(thebools, binBool, csh.getUnusedColor());
				ImageEditing.copyToClipboard(binBool);
				
			}
			break;
		case 7:
//			TopomapUtil.subtractPolynomialFitEachSpectrum(t, Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial")));
			{
				int nbins = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
				int[][] bins = FieldOps.getPercentileBinsForField(Layer.open(fc).data, nbins);
				double[][][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins, Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels (to soften edges)")));
				Topomap smw = new Topomap(smoothWeight, ArrayOps.generateArray(0, 1, nbins), t.x, t.y, null);
				File save = FileOps.selectSave(fc);
				for (int i = 0; i < nbins; i++)
				{
					Topomap.writeBIN(TopomapUtil.spatialFilter(t, smoothWeight[i]), save.toString() + "_" + i + ".bin");
				}
				Topomap.writeBIN(smw, save.toString() + "_WeightsForBins.bin");
			}
			break;
		case 8:
			TopomapUtil.writeAverageSpectraAroundImps(t, imps, Double.parseDouble(JOptionPane.showInputDialog("Enter the gaussian radius in pixels")), fc);
			break;
		case 9:
		{ 	boolean isDefault = JOptionPane.showConfirmDialog(null, "Use Default settings?") == JOptionPane.YES_OPTION;
			
			Layer hist = null;
			if (isDefault)
				hist = TopomapUtil.getSpectralDistributionBasic(1, t.nlayers, t);
			else
			{
				String[] tok = JOptionPane.showInputDialog("Enter the limits of the histogram separated by commas.").split(",");
				double lower = Double.parseDouble(tok[0]);
				double upper = Double.parseDouble(tok[1]);
				int nbins = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins", t.nlayers));
				hist = TopomapUtil.getSpectralDistributionBasicLimited(1, nbins, lower, upper, t);
			}
			int size = t.nlayers;
			while (size*((double)hist.ny/t.nlayers) < 512)
				size*=2;
			
			lv = new LayerViewer(hist, dir, size);
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			break;
		}
		case 10:
			Layer mask = Layer.open(fc);
			for (int i = 0; i < t.nlayers; i++)
				FieldOps.changeToWithoutMask(t.data[i], mask.data);
			break;
		case 11:
			t.data = FFT3D_Wrapper.getFFTMagCentXY(t.data, JOptionPane.showConfirmDialog(null, "Take the log?") == JOptionPane.YES_OPTION);
			break;
		case 12:
			lv = new LayerViewer(TopomapUtil.sumOneDimension(t, Integer.parseInt(JOptionPane.showInputDialog("Enter the code:\r\n\t 0 - Energy \r\n\t 1 - X \r\n\t 2 - Y"))), dir, 1024);
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			break;
		case 13:
			double[] tempSpec;
			double L = Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels."));
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					tempSpec = ArrayOps.gaussSmooth(t.getSpectrum(i, j), L);
					for (int k = 0; k < t.nlayers; k++)
						t.data[k][i][j] = tempSpec[k];
				}
			break;
		case 14:
			int degree = Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial."));
			if (JOptionPane.showConfirmDialog(null, "Save the fits before subtracting?") == JOptionPane.YES_OPTION)
			{
				double[][][] fits = new double[t.nlayers][t.nx][t.ny];
				tempSpec = new double[t.nlayers];
				double[][] tempans;
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
					{
						t.putSpectrum(i, j, tempSpec);
						tempans = ArrayOps.subtractPolynomialFit(t.v, tempSpec, degree);
						for (int k = 0; k < t.nlayers; k++)
						{
							t.data[k][i][j] = tempans[0][k];
							fits[k][i][j] = tempans[1][k];
						}
					}
				
				Topomap.writeBIN(Topomap.newTopomap(t, fits), fc);
			}
			else
			{
				TopomapUtil.subtractPolynomialFitEachSpectrum(t, degree);
			}
			break;
		case 15:
			if (JOptionPane.showConfirmDialog(null, "Take the log?") == JOptionPane.NO_OPTION)
				GraphDrawerCart.plotGraph(ArrayOps.generateArrayNotInclUpper(0, t.nlayers, t.nlayers), TopomapUtil.getAverageFFTOfSpectra(t));
			else{
				double[] temp = TopomapUtil.getAverageFFTOfSpectra(t);
				FieldOps.log(temp);
				GraphDrawerCart.plotGraph(ArrayOps.generateArrayNotInclUpper(0, t.nlayers, t.nlayers), temp);
				
			}
			break;
		case 16:{
			degree = Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial."));
			double[][][] fits = new double[t.nlayers][t.nx][t.ny];
			tempSpec = new double[t.nlayers];
			double[][] tempans = null;
			for (int i = 0; i < t.nx; i++){System.out.print("" + i + " ");
				for (int j = 0; j < t.ny; j++)
				{
					t.putSpectrum(i, j, tempSpec);
					tempans = ArrayOps.subtractPolynomialFit(t.v, tempSpec, degree);
					for (int k = 0; k < t.nlayers; k++)
					{
						t.data[k][i][j] = tempans[0][k];
						fits[k][i][j] = tempans[1][k];
					}
				}
			}
			boolean doEachSpectrum = JOptionPane.showConfirmDialog(null, "To fourier-filter each spectrum, click Yes.\r\n" +
					"To Fourier-filter the map as a whole (i.e. vertical pillars of the 3D FFT) click No, or Cancel.") == JOptionPane.YES_OPTION;
			if (doEachSpectrum){
				
				if (JOptionPane.showConfirmDialog(null, "Save the fits?") == JOptionPane.YES_OPTION)
					Topomap.writeBIN(Topomap.newTopomap(t, fits), fc);
					
				String input = JOptionPane.showInputDialog("Enter upper left corner of the box, spearated by commas.");
				int imin = Integer.parseInt(input.split(",")[0].trim());
				int jmin = Integer.parseInt(input.split(",")[1].trim());
	//			input = JOptionPane.showInputDialog("Enter width and height of the box, separated by commas.");
	//			int width = Integer.parseInt(input.split(",")[0].trim());
	//			int height = Integer.parseInt(input.split(",")[1].trim());
	//			ArrayList<TopomapUtil.FilterBox> box = new ArrayList<TopomapUtil.FilterBox>();
	//			box.add(new TopomapUtil.FilterBox(imin, jmin, width, height));
	//			Topomap ans = TopomapUtil.fourierFilterTheMap(t, box);
	//			for (int i = 0; i < t.nx; i++)
	//				for (int j = 0; j < t.ny; j++)
	//					for (int k = 0; k < t.nlayers; k++)
	//						t.data[k][i][j] = ans.data[k][i][j] + fits[k][i][j];
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++){
						tempans[0] = FFTOps.getFourierFiltered(t.getSpectrum(i, j), new int[] {imin}, new int[] {jmin});
						for (int k = 0; k < t.nlayers; k++)
							t.data[k][i][j] = tempans[0][k] + fits[k][i][j];
					}
			}
			else
			{
				int nranges = Integer.parseInt(JOptionPane.showInputDialog("How many ranges shall we filter?"));
				boolean getFromFile = JOptionPane.showConfirmDialog(null, "Get ranges from a file?") == JOptionPane.YES_OPTION;
				ArrayList<TopomapUtil.FilterBox> boxes = new ArrayList<TopomapUtil.FilterBox>();
				int[] minsX = new int [nranges];
				int[] minsY = new int [nranges];
				int[] widths = new int [nranges];
				int[] heights = new int [nranges];
				if (!getFromFile){
					for (int i = 0; i < nranges; i++)
					{
						String imput2 = JOptionPane.showInputDialog("Enter coordinates of the upper-left corner of the box, comma separated.");
						minsX[i] = Integer.parseInt(imput2.split(",")[0].trim());
						minsY[i] = Integer.parseInt(imput2.split(",")[1].trim());
						imput2 = JOptionPane.showInputDialog("Enter the width and height of the box, comma separated.");
						widths[i] = Integer.parseInt(imput2.split(",")[0].trim());
						heights[i] = Integer.parseInt(imput2.split(",")[1].trim());
						boxes.add(new TopomapUtil.FilterBox(minsX[i], minsY[i], widths[i], heights[i]));
					}
				
					String[] lines = new String[nranges+2];
					lines[0] = "The filter ranges were ";
					lines[1] = "MinX\tMinY\tWidth\tHeight\t";
					for (int i = 0; i < nranges; i++)
					{
						lines[i+2] = minsX[i] + "\t" + minsY[i] + "\t" + widths[i] + "\t" + heights[i];
					}
					JOptionPane.showMessageDialog(null, "Save the record.");
					FileOps.writeLines(fc, lines);
				}
				else
				{
					double[][] table = ColumnIO.readNColumns(FileOps.selectOpen(fc), 4, 2);
					for (int i = 0; i < nranges; i++)
					{
						minsX[i] = (int)table[0][i];
						minsY[i] = (int)table[1][i];
						widths[i] = (int)table[2][i];
						heights[i] = (int)table[3][i];
						boxes.add(new TopomapUtil.FilterBox(minsX[i], minsY[i], widths[i], heights[i]));
						
					}
				}
				
				Topomap temp = TopomapUtil.fourierFilterTheMap(t, boxes);
				FieldOps.copy(temp.data, t.data);
				FieldOps.plusEquals(t.data, fits);
			}

			break;
		}
		case 17:
			TopomapUtil.replaceWithFFTOfSpectra(t);
			break;
		case 18:
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					tempSpec = ArrayOps.getDerivative(t.getSpectrum(i, j));
					for (int k = 0; k < t.nlayers; k++)
						t.data[k][i][j] = tempSpec[k];
				}
			break;
		case 19:
			String input1 = JOptionPane.showInputDialog("Indices to include (exclude upper one)?", "" + 0 + "," + t.nlayers);
			int[] is = new int [] {Integer.parseInt(input1.split(",")[0].trim()),Integer.parseInt(input1.split(",")[1].trim())};
			int nbins = Integer.parseInt(JOptionPane.showInputDialog("How many bins? (Default = Sqrt(number of points)", "" + (int)Math.sqrt((is[1]-is[0])*t.nx*t.ny)));
			double[] bins = ArrayOps.generateArrayInclBoth(ArrayOps.min(t.data), ArrayOps.max(t.data), nbins);
			int[] values = ArrayOps.getHistogram(FieldOps.getArray(t.data, is[0], is[1]), bins);
			GraphDrawerCart.plotGraph(bins, values);
			break;
		case 20:
			Layer gapmap = TopomapUtil.getSimplestGapMap(t, Double.parseDouble(JOptionPane.showInputDialog("Enter the cutoff.")));
			LayerViewer.show(gapmap, 1024, false);
			//plot now the histogram of the gap map properly.
			double dv = (t.v[1]-t.v[0]);
			double[] gapbins = new double [(int)((ArrayOps.max(gapmap.data) - ArrayOps.min(gapmap.data))/dv)];
			for (int i = 0; i < gapbins.length; i++)
				gapbins[i] = -dv/2+i*dv;
			int[] gaphist = ArrayOps.getHistogram(FieldOps.getArray(gapmap.data), gapbins);
			double[] gaphistRelative = new double [gaphist.length];
			
			for (int i = 0; i < gapbins.length; i++){
				gaphistRelative[i] = ((double)gaphist[i]/(double)(t.nx*t.ny));
				System.out.println("" + i + "\t" + (i*dv) + "\t" + gaphistRelative[i]);
			}
			GraphDrawerCart.plotGraph(gapbins, gaphistRelative);
			break;
		case 21:
			int option2 = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice:\r\n"+
					"1 - Smooth the map\r\n"+
					"2 - Subtract a smoothing\r\n" + 
					"3 - Attempt to fix bad pixels\r\n" + 
					"4 - Attempt to Fourier-filter out noise peaks\r\n"));
			switch(option2)
			{
			case 1:
			{
				double li = Double.parseDouble(JOptionPane.showInputDialog("Enter the real-space smoothing length in pixels."));
				double lk = Double.parseDouble(JOptionPane.showInputDialog("Enter the energy smoothing length in pixels."));
				boolean savenow = JOptionPane.showConfirmDialog(null, "Save immediately?") == JOptionPane.YES_OPTION;
				File save = null;
				if (savenow) save = FileOps.selectSave(fc);
				
				FieldOps.copy(FieldOps.getGaussianSmoothing3D(t.data, lk, li, li), t.data);
				if (savenow) Topomap.writeBIN(t, save.toString());
				break;
			}
			case 2:
			{
				boolean openFromFile = JOptionPane.showConfirmDialog(null, "Get sheet from a file?") == JOptionPane.YES_OPTION;
				double li = 0, lk = 0;
				if (!openFromFile){
					li = Double.parseDouble(JOptionPane.showInputDialog("Enter the real-space smoothing length in pixels."));
					lk = Double.parseDouble(JOptionPane.showInputDialog("Enter the energy smoothing length in pixels."));
				}
				FieldOps.minusEquals(t.data, openFromFile ? Topomap.open(fc).data : FieldOps.getGaussianSmoothing3D(t.data, lk, li, li));
				break;
			}
			case 3:
			{
				boolean openFromFile = JOptionPane.showConfirmDialog(null, "Get sheet from a file?") == JOptionPane.YES_OPTION;
				double[][][] smooth = null;
				double li = 0, lk = 0;
				if (openFromFile) smooth = Topomap.open(fc).data;
				else {
					li = Double.parseDouble(JOptionPane.showInputDialog("Enter the real-space smoothing length in pixels."));
					lk = Double.parseDouble(JOptionPane.showInputDialog("Enter the energy smoothing length in pixels."));
				}
				boolean killThePixels = JOptionPane.showConfirmDialog(null, "To KILL the bad pixels click yes. To reduce them, click anything else.") == JOptionPane.YES_OPTION;
				File sheet = null, differences = null;
				boolean saveSheet = openFromFile ? false : JOptionPane.showConfirmDialog(null, "Save the smoothing function?") == JOptionPane.YES_OPTION;
				if (saveSheet) sheet = FileOps.selectSave(fc);
				boolean saveDiff = JOptionPane.showConfirmDialog(null, "Save the differences?") == JOptionPane.YES_OPTION;
				if (saveDiff) differences = FileOps.selectSave(fc);
				double minperc = Double.parseDouble(JOptionPane.showInputDialog("Minimum distrubition cutoff (0 to 1)?"));
				double maxperc = Double.parseDouble(JOptionPane.showInputDialog("Maximum distrubition cutoff (0 to 1)?"));

				if (!openFromFile) smooth = FieldOps.getGaussianSmoothing3D(t.data, lk, li, li);
				FieldOps.minusEquals(t.data, smooth);
				double[][][] diff = FieldOps.cutOffExtremes3D(t.data, minperc, maxperc, killThePixels, saveDiff);
				FieldOps.plusEquals(t.data, smooth);
				
				if (saveSheet)
					Topomap.writeBIN(Topomap.newTopomap(t, smooth), sheet.toString());
				if (saveDiff)
					Topomap.writeBIN(Topomap.newTopomap(t, diff), differences.toString());
				
				JOptionPane.showMessageDialog(null, "Done");
				break;
			}
			case 4:
			{
				boolean openFromFile = JOptionPane.showConfirmDialog(null, "Get sheet from a file?") == JOptionPane.YES_OPTION;
				double[][][] smooth = null;
				double li = 0, lk = 0;
				if (openFromFile) smooth = Topomap.open(fc).data;
				else {
					li = Double.parseDouble(JOptionPane.showInputDialog("Enter the real-space smoothing length in pixels."));
					lk = Double.parseDouble(JOptionPane.showInputDialog("Enter the energy smoothing length in pixels."));
				}
				double[][][] diff = null;
				File sheet = null, differences = null, result = null;
				boolean saveSheet = openFromFile ? false : JOptionPane.showConfirmDialog(null, "Save the smoothing function?") == JOptionPane.YES_OPTION;
				if (saveSheet) sheet = FileOps.selectSave(fc);
				boolean saveDiff = JOptionPane.showConfirmDialog(null, "Save the differences?") == JOptionPane.YES_OPTION;
				if (saveDiff){ differences = FileOps.selectSave(fc);
					diff = new double [t.nlayers][t.nx][t.ny];
				}
				boolean saveNow = JOptionPane.showConfirmDialog(null, "Save the result immediately?") == JOptionPane.YES_OPTION;
				if (saveNow) result = FileOps.selectSave(fc);
				
				boolean doEachSpectrum = JOptionPane.showConfirmDialog(null, "To fourier-filter each spectrum, click Yes.\r\n" +
						"To Fourier-filter the map as a whole (i.e. vertical pillars of the 3D FFT) click No, or Cancel.") == JOptionPane.YES_OPTION;
				
				if (doEachSpectrum){
					int nranges = Integer.parseInt(JOptionPane.showInputDialog("How many ranges shall we filter?"));
					int[] mins = new int [nranges];
					int[] maxs = new int [nranges];
					for (int i = 0; i < nranges; i++)
					{
						String imput2 = JOptionPane.showInputDialog("Enter the minimum and maximum separated by commas. Note maximum will be INCLUDED.");
						mins[i] = Integer.parseInt(imput2.split(",")[0].trim());
						maxs[i] = Integer.parseInt(imput2.split(",")[1].trim());
					}
					if (!openFromFile) smooth = FieldOps.getGaussianSmoothing3D(t.data, lk, li, li);
					FieldOps.minusEquals(t.data, smooth);
					double[] tempSpec2;
					for (int i = 0; i < t.nx; i++)
						for (int j = 0; j < t.ny; j++){
							tempSpec2 = FFTOps.getFourierFiltered(t.getSpectrum(i, j), mins, maxs);
							for (int k = 0; k < t.nlayers; k++)
							{
								if (saveDiff)
									diff[k][i][j] = t.data[k][i][j] - tempSpec2[k];
								t.data[k][i][j] = tempSpec2[k];
							}
						}
				}
				else
				{
					int nranges = Integer.parseInt(JOptionPane.showInputDialog("How many ranges shall we filter?"));
					boolean getFromFile = JOptionPane.showConfirmDialog(null, "Get ranges from a file?") == JOptionPane.YES_OPTION;
					ArrayList<TopomapUtil.FilterBox> boxes = new ArrayList<TopomapUtil.FilterBox>();
					int[] minsX = new int [nranges];
					int[] minsY = new int [nranges];
					int[] widths = new int [nranges];
					int[] heights = new int [nranges];
					if (!getFromFile){
						for (int i = 0; i < nranges; i++)
						{
							String imput2 = JOptionPane.showInputDialog("Enter coordinates of the upper-left corner of the box, comma separated.");
							minsX[i] = Integer.parseInt(imput2.split(",")[0].trim());
							minsY[i] = Integer.parseInt(imput2.split(",")[1].trim());
							imput2 = JOptionPane.showInputDialog("Enter the width and height of the box, comma separated.");
							widths[i] = Integer.parseInt(imput2.split(",")[0].trim());
							heights[i] = Integer.parseInt(imput2.split(",")[1].trim());
							boxes.add(new TopomapUtil.FilterBox(minsX[i], minsY[i], widths[i], heights[i]));
						}
					
						String[] lines = new String[nranges+2];
						lines[0] = "The filter ranges were ";
						lines[1] = "MinX\tMinY\tWidth\tHeight\t";
						for (int i = 0; i < nranges; i++)
						{
							lines[i+2] = minsX[i] + "\t" + minsY[i] + "\t" + widths[i] + "\t" + heights[i];
						}
						JOptionPane.showMessageDialog(null, "Save the record.");
						FileOps.writeLines(fc, lines);
					}
					else
					{
						double[][] table = ColumnIO.readNColumns(FileOps.selectOpen(fc), 4, 2);
						for (int i = 0; i < nranges; i++)
						{
							minsX[i] = (int)table[0][i];
							minsY[i] = (int)table[1][i];
							widths[i] = (int)table[2][i];
							heights[i] = (int)table[3][i];
							boxes.add(new TopomapUtil.FilterBox(minsX[i], minsY[i], widths[i], heights[i]));
							
						}
					}
					if (!openFromFile) smooth = FieldOps.getGaussianSmoothing3D(t.data, lk, li, li);
					FieldOps.minusEquals(t.data, smooth);
					
					Topomap temp = TopomapUtil.fourierFilterTheMap(t, boxes);
					FieldOps.copy(temp.data, t.data);
				}
				
				FieldOps.plusEquals(t.data, smooth);
				
				if (saveSheet)
					Topomap.writeBIN(Topomap.newTopomap(t, smooth), sheet.toString());
				if (saveDiff)
					Topomap.writeBIN(Topomap.newTopomap(t, diff), differences.toString());
				if (saveNow)
					Topomap.writeBIN(t, result.toString());
				
				JOptionPane.showMessageDialog(null, "Done");
				break;
			}
			}
			case 22:
				System.out.println("22");
				double[][] data = new double [t.nx][t.ny];
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
						data[i][j] = ArrayOps.sum(t.getSpectrum(i, j));
				Layer sumlayer = t.getLayer(0);
				sumlayer.data = data;
				LayerViewer.show(sumlayer, 1024, false);
				break;
			
			
//			int nf = Integer.parseInt(JOptionPane.showInputDialog("How many frames?"));
//			String input2 = JOptionPane.showInputDialog("Enter minimum and maximum values of the gap parameter, spearated by commas.");
//			double pppmin = Double.parseDouble(input2.split(",")[0].trim());
//			double pppmax = Double.parseDouble(input2.split(",")[1].trim());
//			double[] ppp = ArrayOps.generateArrayInclBoth(pppmin, pppmax, nf);
//			BufferedImage[] pppstack = new BufferedImage[nf];
//			for (int i = 0; i < nf; i++){
//				Layer xx = TopomapUtil.getSimplestGapMap(t, ppp[i]);
//				pppstack[i] = ImageEditing.getBufferedImage(xx.data, csh.getCurrentCS());
//				Graphics g = pppstack[i].getGraphics();
//				g.setFont(new Font("Arial", Font.PLAIN, 10));
//				g.setColor(csh.getUnusedColor());
//				g.drawString(String.format("%.2f", ppp[i]) + "   " + NumFormat.voltage(ArrayOps.max(xx.data)), 5, t.ny-2);
//				SRAW.writeImage(MovieMaker.avidir + "pppmovie" + MovieMaker.fromInt(i), pppstack[i]);
//			}
//			System.out.println(MovieMaker.BMPtoAVICommand("pppmovie", 0, pppstack.length-1, 5));
//
////			LayerViewer.show(TopomapUtil.getSimplestGapMap(t, Double.parseDouble(JOptionPane.showInputDialog("Enter the cutoff."))), 1024);
//			break;
//			Layer mindidv = TopomapUtil.getSimplestGapMap(t);
			case 23:
			{
				int nbins1 = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
				int[][] bins1 = FieldOps.getPercentileBinsForField(Layer.open(fc).data, nbins1);
				double[][][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins1, Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels (to soften edges)")));
				double[][][] smoothParts = new double [smoothWeight.length][t.nx][t.ny];
				for (int i = 0; i < nbins1; i++)
				{
					smoothParts[i] = FieldOps.spatialFilter(t.data[para], smoothWeight[i]);
				}
				double[] binsArray = ArrayOps.generateArrayNotInclUpper(0, nbins1, nbins1);
				Topomap result = new Topomap(smoothParts, binsArray, t.x, t.y, null);
				Topomap weights = new Topomap(smoothWeight, binsArray, t.x, t.y, null);
				
				Topomap.writeBIN(result, fc);
				Topomap.writeBIN(weights, fc);
				break;
			}
			case 24:
			{
				int nbins1 = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins"));
				Topomap binsource = Topomap.open(fc);
				int[][][] bins1 = new int[t.nlayers][][];
				double[] binsArray = ArrayOps.generateArrayNotInclUpper(0, nbins1, nbins1);
				double L1 = Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels (to soften edges)"));
				int savebin = JOptionPane.showConfirmDialog(null,  "Save the bins?");
				boolean writeEachLayerVsBin = JOptionPane.showConfirmDialog(null, "Would you also like to write a map for each layer, as a function of bin?") == JOptionPane.YES_OPTION;
				File fs = FileOps.selectSave(fc);
				double[][][] smoothWeight = new double [nbins1][t.nx][t.ny];
				double[][][][] smoothParts = new double [nbins1][t.nlayers][t.nx][t.ny];
				for (int j = 0; j < t.nlayers; j++){
					bins1[j] = FieldOps.getPercentileBinsForField(binsource.data[j], nbins1);
					smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins1[j], L1);
					for (int i = 0; i < nbins1; i++)
					{
						smoothParts[i][j] = FieldOps.spatialFilter(t.data[j], smoothWeight[i]);
					}
				}
				if (savebin == JOptionPane.YES_OPTION)
				{
					double[][][] bins2 = new double [t.nlayers][t.nx][t.ny];
					for (int i = 0; i < t.nlayers; i++) bins2[i] = ArrayOps.toDouble(bins1[i]);
					
					Topomap.writeBIN(Topomap.newTopomap(t,bins2), fs.toString() + "bins_unsmoothed.bin");
				}
				if (writeEachLayerVsBin)
				{
					for (int i = 0; i < t.nlayers; i++)
					{
						double[][][] layeri = new double [nbins1][t.nx][t.ny];
						for (int j = 0; j < nbins1; j++)
							layeri[j] = smoothParts[j][i];
						Topomap.writeBIN(new Topomap(layeri, binsArray, t.x, t.y, null), fs.toString() + "layer_" + i + ".bin");					
					}
				}
				
				
				for (int i = 0; i < nbins1; i++)
					Topomap.writeBIN(Topomap.newTopomap(t, smoothParts[i]), fs.toString() + "_" + i + ".bin");
				JOptionPane.showMessageDialog(null, "Done.");
				break;
			}
			case 25:
			{
				double[][][] lap = FieldOps.laplace(t.data);
				for (int i = 0; i < lap.length; i++)
					t.data[i] = lap[i];
				break;
			}
			case 26:
			{
				Topomap m = Topomap.open(fc);
				FieldOps.minusEquals(t.data, m.data);
				break;
			}
			case 27:
			{
				TopomapUtil.saveCropped(t);
				break;
			}
			case 28:
			{
				Printer.copyToClipboard(Printer.arrayVertical(ArrayOps.toArrayLineByLine(t.data[para])));
				break;
			}
			case 29:
			{
				PointSpectra ps = TopomapUtil.getTimeOrderedPointSpectra(t, 0);
				PointSpectra.writeBIN(ps, fc);
				break;
			}
			case 30:
			{
				//export all layers
				String base = FileOps.selectSave(fc).toString();
				String insert = JOptionPane.showInputDialog("Enter what you with to name each layer.", "layer");
				for (int i = 0; i < t.nlayers; i++)
					Layer.writeBIN(t.getLayer(i), base + "_"+insert+"_" + i + ".bin");
				break;
			}
			case 31:
			{
				double initCent = Double.parseDouble(JOptionPane.showInputDialog("Enter your initial guess for the spectrum at the center pixel, in electron volts."));
				double hWidth = Double.parseDouble(JOptionPane.showInputDialog("Enter the half-width of the fitting range in electron volts."));
				Layer[] ans = TopomapUtil.adaptivelyFitAllSpectraToAParabola(t, initCent, hWidth);
				File f = FileOps.selectSave(fc);
				Layer.writeBIN(ans[0], f.toString() + "fitCenter.bin");
				Layer.writeBIN(ans[1], f.toString() + "fitCurvature.bin");
				Layer.writeBIN(ans[2], f.toString() + "rawData.bin");
				Layer.writeBIN(ans[3], f.toString() + "fitOffset.bin");
				break;
			}
			case 32:
			{
				String[] token = JOptionPane.showInputDialog("Enter the Q-vector peak w.r.t. center, in pixels, comma separated.").split(",");
				double[] q = new double [] {Double.parseDouble(token[0]), Double.parseDouble(token[1])};
				double length = Double.parseDouble(JOptionPane.showInputDialog("Enter the length scale in pixles"));
				double[][][] re = new double [t.nlayers][t.nx][t.ny];
				double[][][] im = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++){
					double[][][] temp = FieldOps.getSignalZ_2NM(t.data[i], q, length);
					re[i] = temp[0];
					im[i] = temp[1];
					System.out.print(" " + i);
				}
				new TopomapViewer_complex2(Topomap.newTopomap(t, re), Topomap.newTopomap(t, im),dir, 512);
				System.out.println();
				break;
			}
			case 33:{
				Topomap spf = Topomap.open(fc);
				for (int i = 0; i < t.nlayers; i++)
					t.data[i] = FieldOps.spatialFilter(t.data[i], spf.data[i]);
				break;
			}
			case 34:{
				Layer scGap = TopomapUtil.superconductingGapMap(t);
				LayerViewer.show(scGap, 512, false);
				break;
			}
			case 35:{
				TopomapUtil.getAve(t);
			
				break;
			}

		}
		refreshFFT = true;
		recalculateBounds();
	}
	
	/**
	 * Returns the largest power of two for which p*nsquare is <= 1024
	 * @param nsquare
	 * @return
	 */
	public static int getGoodSize(int nsquare)
	{
		int attempt = 1024/nsquare;
		if (attempt == 1) return nsquare;
		
		double log2 = Math.log(attempt)/Math.log(2);
		return (int)Math.pow(2, (int)log2)*nsquare;
	}
	public static class SliderPanel extends JPanel implements ChangeListener
	{
		TopomapViewer parent;
		JFrame frame;
		public JSlider s;
		public JSlider min, max;
		int parentmember;
		int oldvalue = 0;
		int oldminv = 0, oldmaxv = 999;
		static int npts = 1001;
		
		public SliderPanel(TopomapViewer parent, JFrame frame)
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
			frame.setSize(3*npts+40, 80);
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
				parent.resetCSHFromSlider(((double)min.getValue()/1000), ((double)max.getValue()/1000));
				parent.formImage();
				parent.refresh = true;
				parent.repaint();
				oldminv = min.getValue();
				oldmaxv = max.getValue();
				return;
			}
//			else if (s.getValue() == oldvalue) return;
//			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
			else
			{
				parent.para = s.getValue();
				parent.twoV[0] = parent.t.v[parent.para];
				parent.twoV[1] = parent.t.v[parent.para];
				parent.spec.repaint();
				
				parent.resetGraphics(parent.changeScaleWithLayer, ((double)min.getValue()/1000), ((double)max.getValue()/1000));
				frame.setTitle("" + parent.t.v[parent.para] + "   (Layer " + parent.para + ")");
				System.out.println(s.getValue());
				oldvalue = s.getValue();
			}
		}
		
	}
	
	public BufferedImage[] getImageArrayHighSymPt(){
		if (latt == null) latt = AtomicCoordinatesSet.open(fc);
		int symPt = Integer.parseInt(JOptionPane.showInputDialog("Choose the high-symmetry point:\r\n"
				+ "1 - 1st Bragg peak\r\n"
				+ "2 - 2nd Bragg peak\r\n"
				+ "3 - 1st Rt2 point\r\n"
				+ "4 - 2nd Rt2 point\r\n"
				+ "5 - The origin\r\n"));
		
		int[] rawCent = null;
		if (symPt == 1) rawCent = AtomicCoordinatesSet.generateBragg(latt, t.nx)[0];
		if (symPt == 2) rawCent = AtomicCoordinatesSet.generateBragg(latt, t.nx)[1];
		if (symPt == 3) rawCent = TopomapUtil.MapCuttingMethods.getRt2Point(latt, 0, t);
		if (symPt == 4) rawCent = TopomapUtil.MapCuttingMethods.getRt2Point(latt, 1, t);
		if (symPt == 5) rawCent = new int[] {0, 0};
		
		rawCent[0] += t.nx/2; rawCent[1] += t.ny/2;
		System.out.println("" + rawCent[0] + "\t" + rawCent[1]);
		int size = Integer.parseInt(JOptionPane.showInputDialog("How many pixels do you want in each frame?"));
		
		boolean useCurrentScale = JOptionPane.showConfirmDialog(null, "Lock the color scale?") == JOptionPane.YES_OPTION; 
		String[] tok;
		ColorScale1D scale = null;
		if (useCurrentScale) {
			tok = JOptionPane.showInputDialog("Enter the limits of the color scale, comma separated.", "" + csh.getScale().getMin() + "," + csh.getScale().getMax()).split(",");
			scale = ColorScales.getNew(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), csh.getCurrentCS());
		}
		
		double minperc = useCurrentScale ? 0 : Double.parseDouble(JOptionPane.showInputDialog("Minimum distribution cutoff?"));
		double maxperc = useCurrentScale ? 1 : Double.parseDouble(JOptionPane.showInputDialog("Maximum distribution cutoff?"));
		int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
		int imin, imax;
		String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
		imin = Integer.parseInt(input.split(",")[0].trim());
		imax = Integer.parseInt(input.split(",")[1].trim());
		
		int label = Integer.parseInt(JOptionPane.showInputDialog("Enter your labeling choice\r\n0 - Normal bias voltage \r\n1 - Index label\r\n2 - No label"));
		BufferedImage[] stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, label, minperc, maxperc, 1, useCurrentScale ? scale : null);
		BufferedImage[] substack = new BufferedImage[imax-imin];
		
		for (int i = imin; i < imax; i++)
			substack[i-imin] = stack[i];
		for (int i = 0; i < substack.length; i++)
		{
			substack[i] = ImageEditing.getSubset(substack[i], rawCent[0]-size/2, size, rawCent[1]-size/2, size);
			substack[i] = ImageEditing.getEnlarged(substack[i], resizeFactor);
			ImageEditing.writeBiasVoltage(substack[i], csh.getCurrentCS(), label, t.v[i+imin]);
		}
		return substack;
	}

	public static void main(String[] args)
	{
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		Topomap t = Topomap.open(fc);
		//Changed to make the window display as large as possible while fitting in your display
		//If you want to make it smaller, increase the value of makeWindowSmaller from 0
		java.awt.Dimension screenSize= Toolkit.getDefaultToolkit().getScreenSize();
		int makeSmaller=0;
		int largerDim= Math.max(t.nx, t.ny);
		int scalefactor=(int)(screenSize.getHeight()/largerDim)- makeSmaller;
		int sizechoice = largerDim*scalefactor;
		new TopomapViewer(t, fc.getCurrentDirectory().toString() + "\\", sizechoice);
		System.out.println(t.nx);
	}

}
package util.fileops;

import image.ImageEditing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.border.LineBorder;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;

import drawing.GraphDrawerCart;
import drawing.LayerViewer;
import drawing.TopomapViewer;
import main.SRAW;
import util.ArrayOps;
import util.color.ColorScale;
import util.color.ColorScales;

import util.robot.Robo;


public class MultiFileViewer {

	
	TimeStampedLayer[] everything;
	PointSpectra[] spectra;
	File[] spectraFiles;
	File[] files=null;
	int nf;
	int sf = 0; // selected file
	
	BufferedImage baseImage, drawImage;
	int imageSize;
	Layer coordinateSystem;
	ImageDrawer drawer;
	MiniImageDrawer miniDrawer;
	ExplorerBox exp;
	GraphDrawerCart g = null;
	int miniSize = 384;
	LayerViewer lv; TopomapViewer mv;
	JFileChooser fco;
	JFileChooser fci;
	
	Robo rob = new Robo();
	
	//Menu bar for ImageDrawer
	static JMenuBar menuBar = new JMenuBar();
	
	static int colorScaleIndex = 17;
	
	public MultiFileViewer(int size) {
		this.imageSize = size;
		Topomap.setStdDir();
		fco = new JFileChooser(Topomap.stddir);
		fci = new JFileChooser(NanonisFileOps.nanonisDir);
		baseImage = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
		drawImage = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
	}

	public void populateFiles(String dir)
	{
		sf = 0;
		System.out.println(dir);
		files = new File(dir).listFiles();
		
		if(files == null)
		{
			JOptionPane.showMessageDialog(null, "Choose folder with Nanonis files");
			fco.showSaveDialog(null);
			files = new File(fco.getCurrentDirectory().toString()).listFiles();
		}
		
		int n = 0;
		int m = 0;
		for (int i = 0; i < files.length; i++)
		{
			System.out.println(files[i]);
			if (files[i].toString().endsWith(".sxm") || files[i].toString().endsWith(".3ds"))
				n++;
			if (files[i].toString().endsWith(".dat"))
				m++;
		}
		
		everything = new TimeStampedLayer[n];
		spectra = new PointSpectra[m];
		spectraFiles = new File[m];
		nf = n;
		n = 0;
		m = 0;
		Layer layer;
		for (int i = 0; i < files.length; i++)
		{
			if (files[i].toString().endsWith(".sxm")){
				layer = NanonisFileOps.getLayerZFromScan(files[i]);
				everything[n] = new TimeStampedLayer(files[i].lastModified(), layer, 0, files[i], layer.fileIsComplete || layer.nx*layer.nx < NanonisFileOps.nanThreshold);
				everything[n].lineTime = NanonisFileOps.getLineTimeFromScan(files[i].toString())[0];
			//	everything[n].angle = 0;//NanonisFileOps.getAngleFromFile(files[i].toString());
				n++;
			}
			else if (files[i].toString().endsWith(".3ds")){
				layer = NanonisFileOps.getLayerZFrom3ds(files[i]);
				everything[n] = new TimeStampedLayer(files[i].lastModified(), layer, 0, files[i], layer.fileIsComplete  || layer.nx*layer.nx < NanonisFileOps.nanThreshold);
				//everything[n].angle = 0;//NanonisFileOps.getAngleFromFile(files[i].toString());
				n++;
			}
			if (files[i].toString().endsWith(".dat")){
				spectra[m] = NanonisFileOps.getFromDat(files[i].toString());
				spectraFiles[m++] = files[i];
			}
			
		}
		
		ArrayOps.quicksort(everything);
		SimpleDateFormat form = new SimpleDateFormat();

		for (int i = 0; i < everything.length; i++)
			System.out.println(form.format(new Date(everything[i].timestamp)));
	}
	
	public void populateImage()
	{
		double minX = Double.MAX_VALUE, maxX = -Double.MAX_VALUE;
		double minY = Double.MAX_VALUE, maxY = -Double.MAX_VALUE;
		
		for (int i = 0; i < nf; i++)
		{
			minX = Math.min(minX, ArrayOps.min(everything[i].t.x));
			maxX = Math.max(maxX, ArrayOps.max(everything[i].t.x));
			minY = Math.min(minY, ArrayOps.min(everything[i].t.y));
			maxY = Math.max(maxY, ArrayOps.max(everything[i].t.y));
		}
		
		double dx = maxX - minX, dy = maxY - minY;
		double center;
		if (dx > dy)
		{
			center = (maxY + minY)/2;
			maxY = center + dx/2;
			minY = center - dx/2;
		}
		else if (dy > dx)
		{
			center = (maxX + minX)/2;
			maxX = center + dy/2;
			minX = center - dy/2;
		}
		double[] x = ArrayOps.generateArrayInclBoth(minX, maxX, imageSize);
		double[] y = ArrayOps.generateArrayInclBoth(maxY, minY, imageSize);
		
		double[][] data = new double [imageSize][imageSize];
		coordinateSystem = new Layer(data, x, y, 0, 0);
		
		for (int i = 0; i < imageSize; i++)
			for (int j = 0; j < imageSize; j++)
			{
				baseImage.setRGB(i, j, Color.BLACK.getRGB());
				for (int k = 0; k < nf; k++)
				{
//					if (i % 100 == 0 && k == 0) System.out.println(i);
					if (everything[k].t.containsMetric(x[i], y[j]) && everything[k].showThis)
						baseImage.setRGB(i, j, everything[k].cs.of(everything[k].t.evaluateAtMetricInt(x[i], y[j])).getRGB());
				}
			}
		if (exp != null) exp.dispose();
		exp = null;
		if (drawer != null) drawer.dispose();
		drawer = null;
		if (miniDrawer != null) miniDrawer.dispose();
		miniDrawer = null;

		drawSelectionBox();
		if (drawer == null)
		{
			drawer = new ImageDrawer();
			drawer.setSize(imageSize + 20, imageSize + 80);
			drawer.image = drawImage;
			drawer.parent = this;
			
		    // File menu
		    JMenu fileMenu = new JMenu("File");
		    menuBar.add(fileMenu);
		    
	        //Actions->export as bin
		    JMenuItem exportMI = new JMenuItem("Export selection as a .bin");
		    fileMenu.add(exportMI);
	        class exportAL implements ActionListener{
	        	Robo rob;
	        	exportAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar('x');}
	        }
	        exportMI.addActionListener(new exportAL(rob));
	        
	        //Attempt to convert to topomap
		    JMenuItem topoConvertMI = new JMenuItem("Convert directory to topomap");
		    fileMenu.add(topoConvertMI);
	        class topoConvertAL implements ActionListener{
	        	Robo rob;
	        	topoConvertAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar('m');}
	        }
	        topoConvertMI.addActionListener(new topoConvertAL(rob));
	        
	        //Actions->Spawn a viewer for selected map
		    JMenuItem spawnViewMI = new JMenuItem("Spawn viewer for selected map");
		    fileMenu.add(spawnViewMI);
	        class spawnViewAL implements ActionListener{
	        	Robo rob;
	        	spawnViewAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar('o');}
	        }
	        spawnViewMI.addActionListener(new spawnViewAL(rob));
		    
	        //Selection Menu
	        JMenu selectMenu = new JMenu("Selection");
	        menuBar.add(selectMenu);
	        
		    //Actions->increment selection (note - doing things the cheesy way and typing a key)
		    JMenuItem incSelectMI = new JMenuItem("Select next map (hotkey a)");
		    selectMenu.add(incSelectMI);
	        class incSelectAL implements ActionListener{
	        	Robo rob;
	        	incSelectAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar('a');}
	        }
	        incSelectMI.addActionListener(new incSelectAL(rob));
	        
	        //Actions->decrement selection
		    JMenuItem decSelectMI = new JMenuItem("Select previous map (hotkey d)");
		    selectMenu.add(decSelectMI);
	        class decSelectAL implements ActionListener{
	        	Robo rob;
	        	decSelectAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar('d');}
	        }
	        decSelectMI.addActionListener(new decSelectAL(rob));
	        
	        //Other menu
	        JMenu otherMenu = new JMenu("Other");
	        menuBar.add(otherMenu);
	        
	        //Other-> refresh display
		    JMenuItem refreshMI = new JMenuItem("Refresh display (hotkey spacebar)");
		    otherMenu.add(refreshMI);
	        class refreshAL implements ActionListener{
	        	Robo rob;
	        	refreshAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar(' ');}
	        }
	        refreshMI.addActionListener(new refreshAL(rob));
	        
	        //Other-> do time signal stuff
		    JMenuItem timeSigMI = new JMenuItem("Do time signal stuff");
		    otherMenu.add(timeSigMI);
	        class timeSigAL implements ActionListener{
	        	Robo rob;
	        	timeSigAL(Robo robPrime){rob=robPrime;}
	        	public void actionPerformed(ActionEvent e){rob.typeChar('t');}
	        }
	        timeSigMI.addActionListener(new timeSigAL(rob));
	        
			
			drawer.setJMenuBar(menuBar);
			drawer.setVisible(true);
			drawer.repaint();
			drawer.addKeyListener(drawer);
			drawer.addMouseListener(drawer);
			drawer.addMouseMotionListener(drawer);
			drawer.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		}
		if (miniDrawer == null)
		{
			miniDrawer = new MiniImageDrawer();
			miniDrawer.setSize(miniSize + 20, miniSize + 50);
			miniDrawer.setLocation(imageSize + 20, 0);
			miniDrawer.image = ImageEditing.getBufferedImage_FixedSizeInt(everything[sf].t.data, everything[sf].cs, miniSize, miniSize);
			miniDrawer.setTitle("Selected Item");
			miniDrawer.setVisible(true);
			miniDrawer.repaint();
			miniDrawer.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		}
		if (exp == null)
		{
			exp = new ExplorerBox();
			exp.setSize(1920 - (imageSize + miniSize + 40), imageSize + 50);
			exp.setLocation(imageSize + miniSize + 40, 0);
			exp.parent = this;
			exp.setTitle("File List");
//			exp.setLayout(new BoxLayout(exp, BoxLayout.Y_AXIS));
			exp.initPeerDirList();
			exp.initFileList();
			exp.initSubDirList();
			exp.finalizeInit();
			exp.setVisible(true);
			exp.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		}
	}

	public void exportFile() {
		
		String dest;
		if (everything[sf].f.toString().endsWith(".3ds"))
		{
			dest = FileOps.selectSave(fco).toString();
			boolean cropIt = JOptionPane.showConfirmDialog(null, "Crop any pixels?") == JOptionPane.YES_OPTION;

			if (!cropIt)
				NanonisFileOps.convert3dsToBins(everything[sf].f.toString(), dest);
			else
				NanonisFileOps.convert3dsToBinsCropped(everything[sf].f.toString(), dest);
		}
		else if (everything[sf].f.toString().endsWith(".sxm"))
		{
			dest = FileOps.selectSave(fco).toString();
			NanonisFileOps.convertScanToBins(everything[sf].f.toString(), dest, false);
		}
	}

	public void attemptTopomapConversion() {
		String dir = everything[sf].f.getParent() + "\\";
		if (JOptionPane.showConfirmDialog(null, "Convert the whole directory to topomaps?") == JOptionPane.YES_OPTION)
		{
			NanonisFileOps.convertDirToTopomap2(dir);
		}
		JOptionPane.showMessageDialog(null, "Done");
	}

	public void drawSelectionBox()
	{
		ImageEditing.copy(baseImage, drawImage);
		double[] r0, r1;
		r0 = coordinateSystem.getPixelCoords(ArrayOps.min(everything[sf].t.x), ArrayOps.max(everything[sf].t.y));
		r1 = coordinateSystem.getPixelCoords(ArrayOps.max(everything[sf].t.x), ArrayOps.min(everything[sf].t.y));
		Graphics g = drawImage.getGraphics();
		g.setColor(ColorScales.getUnusedColor(colorScaleIndex));
		g.drawRect((int)r0[0], (int)r0[1], (int)(r1[0]-r0[0]+1), (int)(r1[1]-r0[1]+1));
		String description;
		description = "" + everything[sf].t.nx + " x " + everything[sf].t.ny+ "; ";
		description += String.format("%.1f", (everything[sf].t.xLength*1e9)) + "n X " + String.format("%.1f", (everything[sf].t.xLength*1e9)) + "n; ";
		description += "Bias = " + String.format("%.1f", (everything[sf].t.v*1000)) + "m";
		description += " Setpoint = " + String.format("%.1f", (everything[sf].t.current*1e12)) + "p";
		description += " Line time " + String.format("%.4f", everything[sf].lineTime) + " s"; 
		description += " Angle " + (int)everything[sf].angle;
		if (drawer != null){
			
			drawer.setTitle(description);
			drawer.repaint();
		}
		if (miniDrawer != null) {
			miniDrawer.image = ImageEditing.getBufferedImage_FixedSizeInt(everything[sf].t.data, everything[sf].cs, miniSize, miniSize);
			miniDrawer.repaint();
			miniDrawer.setTitle(new File(everything[sf].f.getParent()).getName() + "\\" + everything[sf].f.getName());
		}
	}
	
	public void spawnViewer() {
		// TODO Auto-generated method stub
		String name = everything[sf].f.getName();
		String message = "What data do you want to view?\r\n";
		if (name.endsWith(".3ds"))
		{
			String[] theNames = NanonisFileOps.getLayerNamesFrom3ds(everything[sf].f.toString());
			message += "0 - The Z Layer\r\n";
			for (int i = 0; i < theNames.length; i++)
				message += "" + (i+1) + " - " + theNames[i] + "\r\n";
			int o = Integer.parseInt(JOptionPane.showInputDialog(message));
			if (o == 0)
			{
				Layer t = NanonisFileOps.getLayerZFrom3ds(everything[sf].f);
				lv = new LayerViewer(t, Topomap.stddir, 512);
				lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			}
			else{
				Topomap[] m = NanonisFileOps.loadTopomapFrom3ds(everything[sf].f.toString());
				mv = new TopomapViewer(m[o-1], Topomap.stddir, 512);
				mv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
				if (theNames[o-1].contains("X"))
					mv.setName("lockin");
				else
					mv.setName(theNames[o-1]);
			}
		}
		if (name.endsWith(".sxm"))
		{
			String[] theNames = NanonisFileOps.loadLayerNamesFromScan(everything[sf].f.toString());
			for (int i = 0; i < theNames.length; i++)
				message += "" + i + " - " + theNames[i] + "\r\n";
			int o = Integer.parseInt(JOptionPane.showInputDialog(message).trim());
			Layer[] m = NanonisFileOps.loadLayerFromScan(everything[sf].f.toString(), true);
			lv = new LayerViewer(m[o], Topomap.stddir, 1024);
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		}
	}
	
	public static class ExplorerBox extends JFrame implements ActionListener
	{
		public MultiFileViewer parent;
		
		public JPanel fileList;
		public JRadioButton[] radio;
		
		public JPanel subDirList;
		public JButton[] sub;
		public JButton up;
		
		public JPanel peerDirList;
		public JRadioButton[] peerRadio;
		int previousPeer;
		public File[] peers;
		
		public JPanel all = null; 
		JScrollPane scroller;
		
		public void initFileList()
		{
			if (all == null)
			{
				all = new JPanel();
				all.setLayout(new BoxLayout(all, BoxLayout.X_AXIS));
			}
			fileList = new JPanel();
			fileList.setLayout(new BoxLayout(fileList, BoxLayout.Y_AXIS));
			fileList.setBorder(new LineBorder(Color.GRAY));
			fileList.add(new JLabel("Files:"));
			JPanel[] lines = new JPanel[parent.nf];
			radio = new JRadioButton[parent.nf];
			ButtonGroup group = new ButtonGroup();
			for (int i = 0; i < parent.nf; i++)
			{
				lines[i] = new JPanel();
				lines[i].setLayout(new BoxLayout(lines[i], BoxLayout.X_AXIS));
				radio[i] = new JRadioButton();
				radio[i].setActionCommand("r" + i);
				radio[i].addActionListener(this);
				group.add(radio[i]);
				
				lines[i].add(radio[i]);
				lines[i].add(new JLabel(parent.everything[i].f.getName()));
				fileList.add(lines[i]);
			}
			radio[0].setSelected(true);
			all.add(fileList);
		}
		
		public void initSubDirList()
		{
			if (all == null)
			{
				all = new JPanel();
				all.setLayout(new BoxLayout(all, BoxLayout.X_AXIS));
			}
			subDirList = new JPanel();
			subDirList.setLayout(new BoxLayout(subDirList, BoxLayout.Y_AXIS));
			
			up = new JButton("Go up one");
			up.setActionCommand("up");
			up.addActionListener(this);
			subDirList.add(up);
			
			int nsd = 0;
			ArrayList<Integer> subDirIs = new ArrayList<Integer>();
			for (int i = 0; i < parent.files.length; i++)
			{
				if (parent.files[i].isDirectory())
				{
					nsd++;
					subDirIs.add(i);
				}
			}
		
			if (nsd == 0){
				all.add(subDirList);
				return;
				
			}
			
			subDirList.add(new JLabel("Subdirectories:"));
			sub = new JButton[nsd];
			int nfs  = 0;
			for (int i = 0; i < nsd; i++)
			{
				sub[i] = new JButton();
				nfs = new File(parent.files[subDirIs.get(i)].toString()).listFiles().length;
				sub[i].setText(parent.files[subDirIs.get(i)].getName() + " - " + nfs + " files");
				
				sub[i].setActionCommand("s" + subDirIs.get(i));
				sub[i].addActionListener(this);
				subDirList.add(sub[i]);
			}
			
			all.add(subDirList);
		}
		
		public void initPeerDirList()
		{
			if (all == null)
			{
				all = new JPanel();
//				all.setLayout(new BoxLayout(all, BoxLayout.X_AXIS));
			}
			
			String upDir = new File(parent.files[0].getParent()).getParent();
			File[] p = new File (upDir).listFiles();
			ArrayList<Integer> pDirIs = new ArrayList<Integer>();
			
			peerDirList = new JPanel();
			peerDirList.setLayout(new BoxLayout(peerDirList, BoxLayout.Y_AXIS));
			peerDirList.add(new JLabel("Peer Directories:"));
			
			int npd = 0;
			for (int i = 0; i < p.length; i++)
			{
				if (p[i].isDirectory())
				{
					npd++;
					pDirIs.add(i);
				}
			}
			
			peers = new File[npd];
			peerRadio = new JRadioButton[npd];
			ButtonGroup group = new ButtonGroup();
			int nfilesP;
			for (int i = 0; i < peers.length; i++)
			{
				peers[i] = p[pDirIs.get(i)];
				peerRadio[i] = new JRadioButton();
				nfilesP = peers[i].listFiles().length;
				peerRadio[i].setText(peers[i].getName() + " - " + nfilesP + " files");
				peerRadio[i].setActionCommand("p" + i);
				peerRadio[i].addActionListener(this);
				if (peers[i].toString().equals(parent.files[0].getParent()))
				{
					peerRadio[i].setSelected(true);
					previousPeer = i;
				}
				group.add(peerRadio[i]);
				peerDirList.add(peerRadio[i]);
				
			}
			
			all.add(peerDirList);
			
		}
		
		public void finalizeInit()
		{
			scroller = new JScrollPane(all);  
			this.add(scroller);
		}
		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			
			if (arg0.getActionCommand().equals("up"))
			{
				String upDir = new File(parent.files[0].getParent()).getParent();
				System.out.println(upDir);
				parent.populateFiles(upDir);
				parent.populateImage();
			}
			
			for (int i = 0; i < parent.nf; i++)
			{
				if (arg0.getActionCommand().equals("r" + i)){
					radio[i].setSelected(true);
					parent.sf = i;
					parent.drawSelectionBox();
				}
			}
			for (int i = 0; i < parent.files.length; i++)
			{
				if (arg0.getActionCommand().equals("s" + i))
				{
					parent.populateFiles(parent.files[i].toString());
					parent.populateImage();
				}
			}
			for (int i = 0; i < peers.length; i++)
			{
				if (arg0.getActionCommand().equals("p" + i))
				{
					if (previousPeer != i)
					{
						previousPeer = i;
						parent.populateFiles(peers[i].toString());
						parent.populateImage();
					}
				}
			}
				
		}
	}
	
	public static class ImageDrawer extends JFrame implements KeyListener, MouseListener, MouseMotionListener
	{
		public BufferedImage image;
		public MultiFileViewer parent;
		Point mouseDownCompCoords = null;
		
		public boolean drawingSpecSites = false;
		public int selectedSpectra = 0;
		public static Color basic = Color.red;
		public static Color special = Color.GREEN;
		
		public void paint(Graphics g)
		{
			g.clearRect(0, 0, this.getWidth(), this.getHeight());
			g.drawImage(image, 10, 55, null);
			
			if (drawingSpecSites)
			{
				for (int i = 0; i < parent.spectra.length; i++){
					double[] specPt = parent.coordinateSystem.getPixelCoords(parent.spectra[i].x[0], parent.spectra[i].y[0]);
					g.setColor(i == selectedSpectra ? special : basic);
					drawPlus(g, (int)Math.round(specPt[0]), (int)Math.round(specPt[1]), 2);
				}
			}
			menuBar.paintImmediately(0,0,this.getWidth(),50);
		}
		public static void drawPlus(Graphics g, int x, int y, int halflength)
		{
			g.drawLine(x-halflength, y, x+halflength, y);
			g.drawLine(x, y-halflength, x, y+halflength);
		}

		@Override
		public void keyTyped(KeyEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void keyReleased(KeyEvent arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public void keyPressed(KeyEvent arg0) {
			// TODO Auto-generated method stub
			switch (arg0.getKeyChar())
			{
			case 'a':
				parent.sf--;
				parent.sf += parent.nf;
				parent.sf %= parent.nf;
				parent.drawSelectionBox();
				parent.exp.radio[parent.sf].setSelected(true);
				break;
			case 'd':
				parent.sf++;
				parent.sf += parent.nf;
				parent.sf %= parent.nf;
				parent.exp.radio[parent.sf].setSelected(true);
				parent.drawSelectionBox();
				break;
			case 'q': if (parent.spectra.length > 0){
				selectedSpectra--;
				selectedSpectra += parent.spectra.length;
				selectedSpectra %= parent.spectra.length;
				repaint();
				if (drawingSpecSites)
				{
					parent.g.setXY(new GraphDrawerCart.GraphObject(parent.spectra[selectedSpectra].v, parent.spectra[selectedSpectra].average), true, true, 0);
					parent.g.repaint();
					parent.g.setTitle(parent.spectraFiles[selectedSpectra].getName());
				}
				}break;
			case 'e':if (parent.spectra.length > 0){
				selectedSpectra++;
				selectedSpectra += parent.spectra.length;
				selectedSpectra %= parent.spectra.length;
				repaint();
				parent.exp.radio[parent.sf].setSelected(true);
				if (drawingSpecSites)
				{
					parent.g.setXY(new GraphDrawerCart.GraphObject(parent.spectra[selectedSpectra].v, parent.spectra[selectedSpectra].average), true, true, 0);
					parent.g.repaint();
					parent.g.setTitle(parent.spectraFiles[selectedSpectra].getName());
				}
				}break;
			case 'x':
				parent.exportFile();
				break;
			case 'o':
				parent.spawnViewer();
				break;
			case ' ':
				repaint();
				break;
			case 'm':
				parent.attemptTopomapConversion();
				break;
			case 'p':
				if (!drawingSpecSites) SRAW.writeImage(FileOps.selectSave(parent.fci).toString(), image);
				else {
					for (int i = 0; i < parent.spectra[selectedSpectra].v.length; i++)
						System.out.println("" + parent.spectra[selectedSpectra].v[i] + "\t" + parent.spectra[selectedSpectra].average[i]);
				}
				break;
			case 'r':
				drawingSpecSites = !drawingSpecSites;
				repaint();
				if (drawingSpecSites){
					parent.g = GraphDrawerCart.getNew(parent.spectra[selectedSpectra].v, parent.spectra[selectedSpectra].average);
					parent.g.setLocation(parent.imageSize + 20, parent.miniSize+50);
					parent.g.showWindow();
					parent.g.setSize(parent.miniSize + 20, parent.miniSize + 50);
					parent.g.setTitle(parent.spectraFiles[selectedSpectra].getName());
				}
				else
				{
					parent.g.dispose();
					parent.g = null;
				}
				break;
			case 't':
				NanonisFileOps.doTimeSignalStuff(parent.everything[parent.sf].f.toString());
				break;
			case 'l'://Treat it as a line cut
				NanonisFileOps.convert3dsToBinsLineCut_PointSpec(parent.everything[parent.sf].f.toString(), FileOps.selectSave(parent.fco).toString());
				break;
			}
		}

		@Override
		public void mouseClicked(MouseEvent e) {
		}
		@Override
		public void mouseEntered(MouseEvent e) {
		}
		@Override
		public void mouseExited(MouseEvent e) {
		}
		@Override
		public void mousePressed(MouseEvent e) {
			mouseDownCompCoords = e.getPoint();
		}
		@Override
		public void mouseReleased(MouseEvent e) {
			mouseDownCompCoords = null;
		}
		@Override
		public void mouseDragged(MouseEvent e) {
			this.setLocation(e.getLocationOnScreen().x - mouseDownCompCoords.x, e.getLocationOnScreen().y - mouseDownCompCoords.y);
		}
		@Override
		public void mouseMoved(MouseEvent e) {
		}
	}
	public static class MiniImageDrawer extends JFrame
	{
		public BufferedImage image;
		public MultiFileViewer parent;
		
		public void paint(Graphics g)
		{
			g.clearRect(0, 0, this.getWidth(), this.getHeight());
			g.drawImage(image, 10, 40, null);
		}

	}
	
	public static class TimeStampedLayer implements Comparable
	{
		public double angle;
		public long timestamp;
		public Layer t;
		public boolean showThis;
		public File f;
		ColorScale cs;
		
		public double lineTime;
		
		public double xShift = 0, yShift = 0;
		
		int type; // 0 for Layer, 1 for map
		
		public TimeStampedLayer(long timestamp, Layer o, int type, File f, boolean showThis) {
			super();
			this.timestamp = timestamp;
			this.t = o;
			this.type = type;
			this. f = f;
			this.showThis = showThis;
//			/*if (t.nx == t.ny)*/ RHKFileOps.doFitting(t, 11);
			cs = ColorScales.getNew(t.data, colorScaleIndex);
			String nan = "NaN";
			this.showThis = this.showThis && !(nan.equals("" + cs.getMax()) || nan.equals("" + cs.getMin()));
		}
		public TimeStampedLayer(File scan) {
			super();
			this. f = scan;
			this.timestamp = f.lastModified();
			this.t = NanonisFileOps.getLayerZFromScan(f);
			this.type = 0;

			this.showThis = true;
//			/*if (t.nx == t.ny)*/ RHKFileOps.doFitting(t, 11);
			cs = ColorScales.getNew(t.data, colorScaleIndex);
			String nan = "NaN";
			this.showThis = this.showThis && !(nan.equals("" + cs.getMax()) || nan.equals("" + cs.getMin()));
		}


		@Override
		/**
		 * Calling arrayOps.quickSort results in an array sorted from earliest to latest.
		 * 
		 */
		public int compareTo(Object arg0) {
			// TODO Auto-generated method stub
			return (int) (timestamp - ((TimeStampedLayer)arg0).timestamp);
		}
	}
	
	public static void main(String[] args)
	{
		MultiFileViewer mfv = new MultiFileViewer(1024);
//		String rawDir = ColumnIO.readString("rawDir.txt");
		mfv.populateFiles(NanonisFileOps.nanonisDir);
//		System.out.println(rawDir);
//		mfv.populateFiles(rawDir);
		mfv.populateImage();
//		System.out.println(rawDir);
//		SRAW.writeImage(FileOps.selectSave(null).toString(), mfv.baseImage);
	}

	public void setEverything(TimeStampedLayer[] o) {
		everything = o;
		nf = o.length;
		files = new File[nf];
		for (int i = 0; i < nf; i++)
			files[i] = o[i].f;
	}

}


package util.robot;

import java.awt.Color;
import java.awt.MouseInfo;
import java.awt.PointerInfo;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import java.sql.Date;
import java.text.SimpleDateFormat;


public class NanonisTopomapTaker implements KeyListener{


	public final Color TIPCOLOR = new Color(220, 0, 44);
	
	public JFrame killFrame;
	public NanonisNumber[] bias;
	public NanonisNumber[] current = null;
	
	public int firstScan = 0; // = 0 for down scan, 1 for up scan
	
	public int[] downScan, upScan, biasBox, currentBox;
	
	public int safetyMargin;
	
	public final String WARNING = "WARNING! This program uses a robot to operate the computer's human interface devices. \r\n"+
			"Do NOT move any window after the mouse locations have been set.\r\n" +
			"Failure to obey all instructions PERFECTLY may result in consequences including but not limited to: \r\n"+
			"           Non-completion or invalidity of the dataset\r\n"+
			"           Tip crash\r\n"+
			"           Shutting down of the computer\r\n" +
			"           Corporal punishment\r\n";
	
	Robo r;
	
	int stage = 0;
	int layerTime, extraTime, initDelay, settleTime;
	public NanonisTopomapTaker()
	{
//		System.out.println(NanonisNumber.parseNum("500m"));
		
		JOptionPane.showMessageDialog(null, WARNING, "Warning", JOptionPane.WARNING_MESSAGE);
	
		
		int nlayers = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of layers the topomap shall have."));
		
		bias = new NanonisNumber[nlayers];
		
		int equallySpaced;
		if (nlayers > 2) equallySpaced = JOptionPane.showConfirmDialog(null, "Shall the layers of the topomap be equally spaced in energy?", "Layer spacing", JOptionPane.YES_NO_OPTION);
		else equallySpaced = JOptionPane.YES_OPTION;
		String bLimits;
		String[] tokens;
		try {
			if (equallySpaced == JOptionPane.YES_OPTION)
			{
				bLimits = JOptionPane.showInputDialog("Enter, in Nanonis format (e.g. 500m) the limits of the voltage range of the topomap.\r\n" +
						"The two limits must be in chronological order, separated by commas.");
				tokens = bLimits.split(",");
				for (int i = 0; i < tokens.length; i++)
					tokens[i] = tokens[i].trim();
				double first = NanonisNumber.parseNum(tokens[0]).toDouble();
				double second = NanonisNumber.parseNum(tokens[1]).toDouble();
				
				double[] numbers = generateArrayInclBoth(first, second, nlayers);
				for (int i = 0; i < numbers.length; i++)
					bias[i] = NanonisNumber.parseDouble(numbers[i]);
			}
			else
			{
				bLimits = JOptionPane.showInputDialog("Enter in Nanonis format (e.g. 500m) the voltage of each layer, in chronological order, separated by commas.");
				tokens = bLimits.split(",");
				for (int i = 0; i < tokens.length; i++)
					tokens[i] = tokens[i].trim();
				if (tokens.length != nlayers)
				{
					JOptionPane.showMessageDialog(null, "You entered the wrong number of voltages. The program will exit.\r\n"+
							"(Expected " + nlayers + "; got " + tokens.length +")");
					System.exit(0);
				}
				for (int i = 0; i < tokens.length; i++)
					bias[i] = NanonisNumber.parseNum(tokens[i]);
			}
		}
		catch(Exception e)
		{
			JOptionPane.showMessageDialog(null, "Something went wrong. The program will exit. (Corporal punishment?)");
			System.exit(0);
		}
		
		String biases = "";
		for (int i = 0; i < nlayers; i++)
			biases += bias[i].toString() + "\r\n";
		
		JOptionPane.showMessageDialog(null, "The following numbers will be typed in the Bias box: \r\n" + biases);
			
		
		int currentOption;
		if (nlayers > 2) currentOption = JOptionPane.showConfirmDialog(null, "Do you want to use a different setpoint for each layer?", "vary the current?", JOptionPane.YES_NO_OPTION);
		else currentOption = JOptionPane.YES_OPTION;
		if (currentOption == JOptionPane.YES_OPTION)
		{
			current = new NanonisNumber[bias.length];
			equallySpaced = JOptionPane.showConfirmDialog(null, "Shall setpoint values be equally spaced?", "Setpoint spacing", JOptionPane.YES_NO_OPTION);
			
			try {
				if (equallySpaced == JOptionPane.YES_OPTION)
				{
					bLimits = JOptionPane.showInputDialog("Enter, in Nanonis format (e.g. 50p) the limits of the setpoint range.\r\n" +
							"The two limits must be in chronological order, separated by commas.");
					tokens = bLimits.split(",");
					for (int i = 0; i < tokens.length; i++)
						tokens[i] = tokens[i].trim();
					double first = NanonisNumber.parseNum(tokens[0]).toDouble();
					double second = NanonisNumber.parseNum(tokens[1]).toDouble();
					
					double[] numbers = generateArrayInclBoth(first, second, nlayers);
					for (int i = 0; i < numbers.length; i++)
						current[i] = NanonisNumber.parseDouble(numbers[i]);
				}
				else
				{
					bLimits = JOptionPane.showInputDialog("Enter in Nanonis format (e.g. 50p) the setpoint for each layer, in chronological order, separated by commas.");
					tokens = bLimits.split(",");
					for (int i = 0; i < tokens.length; i++)
						tokens[i] = tokens[i].trim();
					if (tokens.length != nlayers)
					{
						JOptionPane.showMessageDialog(null, "You entered the wrong number of voltages. The program will exit.\r\n"+
								"(Expected " + nlayers + "; got " + tokens.length +")");
						System.exit(0);
					}
					for (int i = 0; i < tokens.length; i++)
						current[i] = NanonisNumber.parseNum(tokens[i]);
				}
			}
			catch(Exception e)
			{
				JOptionPane.showMessageDialog(null, "Something went wrong. The program will exit. (Corporal punishment?)");
				e.printStackTrace();
				System.exit(0);
			}
			String currents = "";
			for (int i = 0; i < nlayers; i++)
				currents += current[i].toString() + "\r\n";
			
			JOptionPane.showMessageDialog(null, "The following numbers will be typed in the Setpoint box: \r\n" + currents);
			
		}
	}
	public static double[] generateArrayInclBoth(double min, double max, int npts)
	{
		if (npts == 1) return new double[] {min};
		double dm = (max-min)/(npts-1);
		double [] a = new double [npts];
		for (int i = 0; i < a.length; i++)
			a[i] = min + i*dm;
		return a;
	}

	private void getCoordinatesOfPoints()
	{
		killFrame = new JFrame();
		killFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		killFrame.addKeyListener(this);
		killFrame.show();
		
		String intro = "In this part you will mouse over various objects and press spacebar to register their location.\r\n" +
				"It is absolutely essential that you get this right since otherwise the mouse will click on the wrong things.\r\n" + 
				"You can abort at any time BEFORE starting the map by closing the Java window which has just appeared." ;
		
		JOptionPane.showMessageDialog(null, intro);
		
		guideUser();
	}
	public void guideUser()
	{
		switch(stage)
		{
		case 0:
			JOptionPane.showMessageDialog(null, "(After closing this box) Move the mouse over the DOWN SCAN button and press Spacebar");
			break;
		case 1:
			JOptionPane.showMessageDialog(null, "(After closing this box) Move the mouse over the UP SCAN button and press Spacebar");
			break;
		case 2:
			JOptionPane.showMessageDialog(null, "(After closing this box) Move the mouse over the Bias input box and press Spacebar");
			break;
		case 3:
			if (current != null) 
				JOptionPane.showMessageDialog(null, "(After closing this box) Move the mouse over the Setpoint input box and press Spacebar");
			else
				takeTime();
			break;
		case 4:
			takeTime();
			break;
//		case 4:
//			JOptionPane.showMessageDialog(this, "Now move the mouse over the Save button and press a key.");
//			break;
		}
	}

	private void takeTime()
	{
		try{
			String t = JOptionPane.showInputDialog(null, "Enter the time which NANONIS estimates a single scan will take, in the form mm:ss");
			String[] tokens = t.split(":");
			layerTime = Integer.parseInt(tokens[0])*60 + Integer.parseInt(tokens[1]);
			
			t = JOptionPane.showInputDialog(null, "Enter, in seconds, the amount of ADDITIONAL time you want to wait after each scan (to be safe).");
			extraTime = (int)Double.parseDouble(t);
			
			t = JOptionPane.showInputDialog(null, "It's a good idea to let the tip stabilze for a bit after we change the bias" + (current != null ? " (and setpoint).\r\n" : ".\r\n")+
			"Enter the time you want to wait in this manner, in seconds.", "Stabalization time", JOptionPane.INFORMATION_MESSAGE);
			
			settleTime = (int)Double.parseDouble(t);
			
			int o = JOptionPane.showConfirmDialog(null, "To start with a down scan, click Yes.\r\n" +
					"To start with an up scan, click No.", "Scan Direction", JOptionPane.YES_NO_OPTION);
			if (o == JOptionPane.NO_OPTION)
				firstScan = 1;
		
			
			long now = System.currentTimeMillis();
			long dn = now + bias.length*(settleTime + layerTime + extraTime)*1000;
			Date done = new Date(dn);
			SimpleDateFormat form = new SimpleDateFormat();
			SimpleDateFormat day = new SimpleDateFormat("EE");
			String updown = firstScan == 1 ? "an up scan" : "a down scan";
			t = JOptionPane.showInputDialog(null, "The map will begin shortly, with " + updown + ".\r\n" +
			"We estimate that it will finish on " + day.format(done) + " " + form.format(done) + ".\r\n" +
			"Enter in the box below, in seconds, the delay between when you press OK and when the first scan will start.\r\n"+
					"(To abort, press Ctrl+Shift+Esc to open the Task Manager. In the process tab, find javaw.exe and End Process.)");
			
			
			initDelay = (int)Double.parseDouble(t);
			
		}
		catch (Exception e)
		{
			JOptionPane.showMessageDialog(null, "You made a mistake in input. The program will exit.");
			System.exit(0);
		}
		takeTheMap();
	}
	private void takeTheMap()
	{
		r = new Robo();
		Robo.wait(initDelay*1000);
		
		for (int i = 0; i < bias.length; i++)
		{
			killFrame.setTitle("Layer " + (i+1) + " of " + bias.length);
			takeOneLayer(i);
			Robo.wait((layerTime+extraTime)*1000);
		}
		
		JOptionPane.showMessageDialog(null, "The map is done.");
		System.exit(0);
	}
	
	private void takeOneLayer(int i)
	{
		r.goTo(biasBox);
		r.wait(50);
		r.doubleclick();
		r.wait(50);
//		r.typeChar('a', true);
		r.typeString(bias[i].toString());
		r.pressEnter();
		
		if (current != null)
		{
			r.goTo(currentBox);
			r.wait(50);
			r.doubleclick();
			r.wait(50);
//			r.typeChar('a', true);
			r.typeString(current[i].toString());
			r.pressEnter();
		}
		
		r.wait(settleTime*1000);
		if (i %2 == firstScan)
			r.goTo(downScan);
		else
			r.goTo(upScan);
		r.wait(50);
		r.click();
	}
	
	
	public static void main(String[] args)
	{
		NanonisTopomapTaker nt = new NanonisTopomapTaker();
		nt.getCoordinatesOfPoints();
	}

	public static class NanonisNumber
	{
		double num;
		String sfx;
		
		public NanonisNumber(double num, String sfx)
		{
			this.num = num;
			this.sfx = sfx;
		}
		
		public static NanonisNumber parseNum(String input)
		{
			if (input.endsWith("f") || input.endsWith("p") || input.endsWith("n") || input.endsWith("u") || input.endsWith("m"))
				return new NanonisNumber(Double.parseDouble(input.substring(0, input.length()-1)), input.substring(input.length()-1));
			else return new NanonisNumber(Double.parseDouble(input), "");
		}
		public static NanonisNumber parseDouble(double d)
		{
			if (Math.abs(d) < 1e-12)
				return new NanonisNumber(d*1e15, "f");
			if (Math.abs(d) < 1e-9)
				return new NanonisNumber(d*1e12, "p");
			if (Math.abs(d) < 1e-6)
				return new NanonisNumber(d*1e9, "n");
			if (Math.abs(d) < 1e-3)
				return new NanonisNumber(d*1e6, "u");
			if (Math.abs(d) < 1)
				return new NanonisNumber(d*1000, "m");
			else return new NanonisNumber(d, "");

		}
		
		public String toString()
		{
			return "" + String.format("%.3f", num) + sfx;
		}
		
		public double toDouble()
		{
			if (sfx.equals("m"))
				return num/1000.0;
			if (sfx.equals("u"))
				return num/1e6;
			if (sfx.equals("n"))
				return num/1e9;
			if (sfx.equals("p"))
				return num/1e12;
			if (sfx.equals("f"))
				return num/1e15;
			return num;
		}
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
		PointerInfo p = MouseInfo.getPointerInfo();
		int x = (int)p.getLocation().x;
		int y = (int)p.getLocation().y;
		switch(stage)
		{
		case 0:
			downScan = new int [] {x, y};
			break;
		case 1:
			upScan = new int [] {x, y};
			break;
		case 2:
			biasBox = new int [] {x, y};
			break;
		case 3:
			currentBox = new int [] {x, y};
			break;
		}
		stage++;
		guideUser();

	}
}

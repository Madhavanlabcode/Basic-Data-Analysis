package util.fileops;

import java.awt.image.BufferedImage;
import java.io.File;

import javax.swing.JFileChooser;

import main.SRAW;

public class FileOps {

	public static String openText(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			return ColumnIO.readString(fc.getSelectedFile().toString());
		else return null;
	}
	public static String[] openLines(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			return ColumnIO.readLines(fc.getSelectedFile());
		else return null;
	}
	public static String openText()
	{
		return openText(new JFileChooser(Topomap.stddir));
	}
	
	public static String getDir(String fullpath)
	{
		return fullpath.substring(0, fullpath.lastIndexOf("\\")+1);
	}

	public static double[][] openBin(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			return ColumnIO.readSquareTable(fc.getSelectedFile().toString());
		else return null;
	}
	public static double[][] openBin(String dir)
	{
		if (dir.equals(""))
			return openBin(new JFileChooser(Topomap.stddir));
		else return openBin(new JFileChooser(dir));
	}
	
	public static double[][] openTable(String dir)
	{
		JFileChooser fc;
		if (dir.equals(""))
			fc = new JFileChooser(Topomap.stddir);
		else
			fc = new JFileChooser(dir);
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			if (fc.getSelectedFile().toString().endsWith(".txt"))
				return ColumnIO.readAllColumns(fc.getSelectedFile(), null);
			else
				return ColumnIO.readSquareTable(fc.getSelectedFile().toString());
		else return null;
	}
	public static double[][] openTable(JFileChooser fc)
	{
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			if (fc.getSelectedFile().toString().endsWith(".txt"))
				return ColumnIO.readAllColumns(fc.getSelectedFile(), null);
			else
				return ColumnIO.readSquareTable(fc.getSelectedFile().toString());
		else return null;
	}
	public static double[][][] open2Comp(JFileChooser fc, boolean polar)
	{
		File f1 = null, f2 = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f1 = fc.getSelectedFile();
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f2 = fc.getSelectedFile();
		return ColumnIO.readSquareTables(f1.toString(), f2.toString(), polar);
	}
	public static double[][][] open2Comp(boolean polar)
	{
		JFileChooser fc = new JFileChooser (Topomap.stddir);
		File f1 = null, f2 = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f1 = fc.getSelectedFile();
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f2 = fc.getSelectedFile();
		return ColumnIO.readSquareTables(f1.toString(), f2.toString(), polar);
	}
	

	public static void writeTableASCII(double[][] data)
	{
		writeTableASCII(new JFileChooser(Topomap.stddir), data);
	}
	public static void writeTableBIN(double[][] data)
	{
		writeTableBIN(new JFileChooser(Topomap.stddir), data);
	}
	public static void writeString(String s)
	{
		writeString(new JFileChooser(Topomap.stddir), s);
	}
	public static void writeLines(String[] lines)
	{
		writeLines(new JFileChooser(Topomap.stddir), lines);
	}
	public static void writeTableASCII(JFileChooser fc, double[][] data)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			ColumnIO.writeTable(data, fc.getSelectedFile().toString() + (fc.getSelectedFile().toString().endsWith(".txt") ? "" : ".txt"));		
	}
	public static void writeTableBIN(JFileChooser fc, double[][] data)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			ColumnIO.writeBin(data, fc.getSelectedFile().toString());		
	}
	public static void writeImage(JFileChooser fc, double[][] data)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			SRAW.writeImage(fc.getSelectedFile().toString(), data);		
	}
	public static void writeImage(JFileChooser fc, BufferedImage img)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			SRAW.writeImage(fc.getSelectedFile().toString(), img);		
	}

	public static void writeString(JFileChooser fc, String s)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			ColumnIO.writeString(s, fc.getSelectedFile().toString());		
	}
	public static void writeLines(JFileChooser fc, String[] lines)
	{
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			ColumnIO.writeLines(lines, fc.getSelectedFile().toString().endsWith(".txt") ? fc.getSelectedFile().toString() : fc.getSelectedFile().toString() + ".txt");		
	}
	public static void writeBINandImage(String s, double[][] data)
	{
		ColumnIO.writeBin(data, s + ".dat");
		SRAW.writeImage(s, data);
	}
	public static void writeBINandImage(String s, int[][] data)
	{
		ColumnIO.writeBin(data, s + ".dat");
		SRAW.writeImage(s, data);
	}
	
	public static File selectOpen(JFileChooser fc)
	{
		if (fc == null) fc = new JFileChooser(Topomap.stddir);
		File f = null;
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION)
			f = fc.getSelectedFile();
		return f;
	}
	public static File selectSave(JFileChooser fc)
	{
		if (fc == null) fc = new JFileChooser(Topomap.stddir);
		File f = null;
		if (fc.showSaveDialog(null) == JFileChooser.APPROVE_OPTION)
			f = fc.getSelectedFile();
		return f;
	}
	public static String selectDir(JFileChooser fc)
	{
		if (fc == null) fc = new JFileChooser(Topomap.stddir);
		File f = null;
		fc.showSaveDialog(null);
		f = fc.getCurrentDirectory();
		return f.toString() + "\\";
		
	}
	public static String addExtraBackslashes(File f)
	{
		String s = f.toString();
//		String[] g = s.split("\\");
//		String x = "";
//		int i = 0;
//		for (i = 0; i < g.length-1; i++)
//			x += g[i] + "\\\\";
//		x += g[i];
		String x = s.replace("\\", "\\\\");
		return x;
	}
	public static JFileChooser getFC()
	{
		return new JFileChooser(Topomap.stddir);
	}
	public static String getDoubleName(File f)
	{
		return new File(f.getParent()).getName() + "\\" + f.getName();
	}
	public static void mkdir(String dir) {
		if (!new File(dir).exists()) new File(dir).mkdirs();
	}
}

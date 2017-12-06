package util;

import java.awt.Toolkit;
import java.awt.datatransfer.StringSelection;

import javax.swing.JOptionPane;

public class Printer {
	public static void printlnHorizontal(double[] x)
	{
		for (int i = 0; i < x.length-1; i++)
			System.out.print(x[i] + "\t");
		System.out.println(x[x.length-1]);
	}
	public static void printlnHorizontal(int[] x)
	{
		for (int i = 0; i < x.length-1; i++)
			System.out.print(x[i] + "\t");
		System.out.println(x[x.length-1]);
	}
	public static void printlnVertical(double[] x)
	{
		for (int i = 0; i < x.length; i++)
			System.out.println(x[i]);
	}
	public static void printlnVertical(int[] x)
	{
		for (int i = 0; i < x.length; i++)
			System.out.println(x[i]);
	}
	public static void printOrderedPairs(double[][][] ans)
	{
		
		for (int i = 0; i < ans[0].length; i++)
		{
			for (int j = 0; j < ans[0][i].length; j++)
			{
				System.out.print("" + ans[0][i][j] + "\t" + ans[1][i][j]);
			}
		}
//		return ans;
		
	}
	public static String arrayVertical(double[] x)
	{
		String s = "";
		for (int i = 0; i < x.length-1; i++)
			s += x[i] + "\r\n";
		s += x[x.length-1];
		return s;
	}
	public static String arrayLnHorizontal(double[] x)
	{
		String s = "";
		for (int i = 0; i < x.length-1; i++)
			s += x[i] + "\t";
		s += x[x.length-1] + "\r\n";
		return s;
	}
	public static String vectorP(double[] x)
	{
		String s = "(";
		for (int i = 0; i < x.length-1; i++)
			s += x[i] + ", ";
			s += x[x.length-1]+")";
		return s;
	}
	public static String vectorPFormat(double[] x)
	{
		String s = "(";
		for (int i = 0; i < x.length-1; i++)
			s += NumFormat.scientific(x[i],3) + ", ";
			s += NumFormat.scientific(x[x.length-1],3)+")";
		return s;
	}
	public static String vectorP(double x, double y)
	{
		return "(" + x + ", " + y + ")";
	}
	public static String vectorP(int[] x)
	{
		return "(" + x[0] + ", " + x[1] + ")";
	}
	public static String[] getHistLines(double[] bins, int[] hist) {
		String[] lines = new String[bins.length];
		for (int i = 0; i < bins.length; i++)
			lines[i] = bins[i] + "\t" + hist[i];
		return lines;
	}
	
	public static String[] getColumnSeries(double[][] x, double [][] y)
	{
		int maxHeight = -1;
		for (int i = 0; i < x.length; i++)
			maxHeight = Math.max(x[i].length, maxHeight);
		
		String[] ans = new String[maxHeight];
		for (int i = 0; i < maxHeight; i++)
		{
			ans[i] = "";
			for (int j = 0; j < x.length; j++)
			{
				ans[i] += i < x[j].length ? "" + x[j][i] + "\t" + y[j][i] + "\t" : "\t\t";
			}
		}
		return  ans;
	}
	public static void printColumnSeries(double[][] x, double [][] y)
	{
		String[] lines = getColumnSeries(x, y);
		for (int i = 0; i < lines.length; i++)
			System.out.println(lines[i]);
		
		//optional:
		copyToClipboard(lines);
	}
	public static void copyToClipboard(String s)
	{
		Toolkit.getDefaultToolkit ().getSystemClipboard ().setContents(new StringSelection(s), null);
	}
	public static void copyToClipboard(String[] lines){
		String sum = "";
		for (int i = 0; i < lines.length; i++)
			sum += lines[i] + "\r\n";
		
		StringSelection ss = new StringSelection(sum);
		Toolkit.getDefaultToolkit ().getSystemClipboard ().setContents (ss, null);
		
	}
	public static String getTable(double[][] a) {
		String ans = "";
		for (int i = 0; i < a[0].length; i++){
			for (int j = 0; j < a.length; j++)
				ans += a[j][i] + "\t";
			if (i < a[0].length - 1)
				ans +="\r\n" ;
		}
		return ans;
	}
	public static int[] getTwoInts(){
		String s = JOptionPane.showInputDialog("Enter two integers separated by commas.");
		String[] tok = s.split(",");
		return new int[] {Integer.parseInt(tok[0]), Integer.parseInt(tok[1])};
	}
	public static int[] getTwoInts(String title){
		String s = JOptionPane.showInputDialog(null, title, "Enter the integers separated by commas.", JOptionPane.QUESTION_MESSAGE);
		String[] tok = s.split(",");
		return new int[] {Integer.parseInt(tok[0]), Integer.parseInt(tok[1])};
	}
	public static int getAnInt(String title){
		String s = JOptionPane.showInputDialog(null, title, "Enter an integer.", JOptionPane.QUESTION_MESSAGE);
		return Integer.parseInt(s);
	}
	public static int getAnInt(String title, int initChoice){
		String s = JOptionPane.showInputDialog(title, "" + initChoice);
		return Integer.parseInt(s);
	}
	public static double getADouble(String title){
		String s = JOptionPane.showInputDialog(null, title, "Enter an integer.", JOptionPane.QUESTION_MESSAGE);
		return Double.parseDouble(s);
	}
	public static double[] getTwoDoubles(String title){
		String s = JOptionPane.showInputDialog(null, title, "Enter the numbers, separated by commas.", JOptionPane.QUESTION_MESSAGE);
		String[] tok = s.split(",");
		return new double[] {Double.parseDouble(tok[0]), Double.parseDouble(tok[1])};
	}
}

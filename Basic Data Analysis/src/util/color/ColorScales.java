package util.color;

import java.awt.Color;

import util.ArrayOps;
import util.Complex;

public class ColorScales {
	//Color scale 1D: table
	// 0 - BRYW
	// 1 - Grey
	// 2 - BlueWhiteRed
	// 3 - BlackBlueWhite
	//Basic color scales: Black to Primary Color to White. (512 points.)
	
	public static Color bryw[] = null;
	public static Color myc2d[][] = null;
	public static Color rgb2d[][] = null; //Red = 0, g = 2pi/3, b = 4pi/3.
	public static Color blueWhiteRed[] = null;
	public static Color grey[] = null;
	public static Color[] blackBlueWhite = null;
	public static Color[][] blackPrimaryWhite = new Color[][] {null, null, null};
	public static Color[][] blackPrimarySecondaryWhite = new Color[][] {null, null, null, null, null, null};
	public static Color[] greyBounded = null;
	//This one is based on 512 values. 1st color goes from 0 to 255 in the first 256. 3rd goes in last 256. Middle color grows half as fast through all 512.
	public static Color[][] nanonisBlueType = new Color[][] {null, null, null, null, null, null};
	public static Color[] rainbow = null;
	public static Color[] rhkDefault = null;
	public static Color[] purpleWhiteOrange = null;
	
	public static final int NSCALES = 21;
	static {
		initBRYW();
		initMYC2d();
	}
	
	public static void initBRYW()
	{
		bryw = new Color[768];
		for (int i = 0; i < 768; i++)
		{
			if (i < 256)
				bryw[i] = new Color(i, 0, 0); 
			else if (i < 512)
				bryw[i] = new Color(255, i-256, 0);
			else
				bryw[i] = new Color(255, 255, i-512);
		}
	}
	public static void initMYC2d(){
		myc2d = new Color[256][768];
		for (int j = 0; j < 256; j++)
			for (int i = 0; i < 768; i++)
			{
				if (i < 256)
					myc2d[j][i] = new Color(j, (j*i)/255, (j*(255 - i))/255); 
				else if (i < 512)
					myc2d[j][i] = new Color((j*(511 - i))/255, j, (j*(i-256))/255);
				else
					myc2d[j][i] = new Color((j*(i-512))/255, (j*(767 - i))/255, j);
			}
	}
	public static void initRGB2d(){
		rgb2d = new Color[256][768];
		for (int j = 0; j < 256; j++)
			for (int i = 0; i < 768; i++)
			{
				if (i < 256)
					myc2d[j][i] = new Color((j*(255 - i))/255, (j*i)/255, 0); 
				else if (i < 512)
					myc2d[j][i] = new Color(0, (j*(511 - i))/255, (j*(i-256))/255);
				else
					myc2d[j][i] = new Color((j*(i-512))/255, 0, (j*(767 - i))/255);
			}
	}
	public static void initBlueWhiteRed(){
		blueWhiteRed = new Color[512];
		for (int i = 0; i < 512; i++)
		{
			if (i < 256)
				blueWhiteRed[i] = new Color(i, i, 255); 
			else
				blueWhiteRed[i] = new Color(255, 511-i, 511-i);
		}
	}
	public static void initPurpleWhiteOrange(){
		int nc = 512;
		purpleWhiteOrange = new Color[nc];
		double[] r0 = ArrayOps.generateArrayNotInclUpper(96, 255, nc/2);
		double[] g0 = ArrayOps.generateArrayNotInclUpper(0, 255, nc/2);
		double[] b0 = ArrayOps.generateArrayNotInclUpper(192, 255, nc/2);
		double[] r1 = ArrayOps.generateArrayNotInclUpper(255, 192, nc/2);
		double[] g1 = ArrayOps.generateArrayNotInclUpper(255, 96, nc/2);
		double[] b1 = ArrayOps.generateArrayNotInclUpper(255, 0, nc/2);
		for (int i = 0; i < 256; i++)
			purpleWhiteOrange[i] = new Color((int)r0[i], (int)g0[i], (int)b0[i]);
		for (int i = 0; i < 256; i++)
			purpleWhiteOrange[i+256] = new Color((int)r1[i], (int)g1[i], (int)b1[i]);
//			if (i < 256)
//				purpleWhiteOrange[i] = new Color(i/2+128, i, 255); 
//			else
//				purpleWhiteOrange[i] = new Color(255, 383-i/2, 511-i);
//			if (i < 256)
//				purpleWhiteOrange[i] = new Color(128, i/2, 255 - i/2); 
//			else
//				purpleWhiteOrange[i] = new Color(i/2, 128, 255-i/2);
//		}
	}
	
	public static void initGrey()
	{
		grey = new Color [256];
		for (int i = 0; i < grey.length; i++)
			grey[i] = new Color(i, i, i);
	}
	public static void initRHK()
	{
		rhkDefault = getFromHSVValues(2, 59, 0.96769, 0.90137, 0.32778, 1, 256);
	}
	public static void initGreyBounded()
	{
		greyBounded = new Color [256];
		for (int i = 1; i < 255; i++)
			greyBounded[i] = new Color(i, i, i);
		greyBounded[0] = Color.blue;
		greyBounded[255] = Color.red;
	}
	public static void initBlackBlueWhite()
	{
		blackBlueWhite = new Color[512];
		for (int i = 0; i < 768; i++)
		{
			if (i < 256)
				blackBlueWhite[i] = new Color(0, 0, i); 
			else if (i < 512)
				blackBlueWhite[i] = new Color(i-256, i-256, 255);
		}
	}
	public static void initBlackPrimaryWhite()
	{
		blackPrimaryWhite = new Color[3][512];
		for (int i = 0; i < 768; i++)
		{
			if (i < 256)
			{
				blackPrimaryWhite[0][i] = new Color(i, 0, 0); 
				blackPrimaryWhite[1][i] = new Color(0, i, 0); 
				blackPrimaryWhite[2][i] = new Color(0, 0, i); 
			}
			else if (i < 512)
			{
				blackPrimaryWhite[0][i] = new Color(255, i-256, i-256); 
				blackPrimaryWhite[1][i] = new Color(i-256, 255, i-256); 
				blackPrimaryWhite[2][i] = new Color(i-256, i-256, 255); 
			}
		}
	}
	public static void initBlackPrimarySecondaryWhite()
	{
		blackPrimarySecondaryWhite = new Color[6][768];
		for (int i = 0; i < 768; i++)
		{
			if (i < 256)
			{
				blackPrimarySecondaryWhite[0][i] = new Color(i, 0, 0); 
				blackPrimarySecondaryWhite[1][i] = new Color(i, 0, 0); 
				blackPrimarySecondaryWhite[2][i] = new Color(0, i, 0); 
				blackPrimarySecondaryWhite[3][i] = new Color(0, i, 0); 
				blackPrimarySecondaryWhite[4][i] = new Color(0, 0, i); 
				blackPrimarySecondaryWhite[5][i] = new Color(0, 0, i); 
			}
			else if (i < 512)
			{
				blackPrimarySecondaryWhite[0][i] = new Color(255, i-256, 0);
				blackPrimarySecondaryWhite[1][i] = new Color(255, 0, i-256);
				blackPrimarySecondaryWhite[2][i] = new Color(i-256, 255, 0); 
				blackPrimarySecondaryWhite[3][i] = new Color(0, 255, i-256); 
				blackPrimarySecondaryWhite[4][i] = new Color(0, i-256, 255); 
				blackPrimarySecondaryWhite[5][i] = new Color(i-256, 0, 255); 
			}
			else
			{
				blackPrimarySecondaryWhite[0][i] = new Color(255, 255, i-512);
				blackPrimarySecondaryWhite[1][i] = new Color(255, i-512, 255);
				blackPrimarySecondaryWhite[2][i] = new Color(255, 255, i-512); 
				blackPrimarySecondaryWhite[3][i] = new Color(i-512, 255, 255); 
				blackPrimarySecondaryWhite[4][i] = new Color(i-512, 255, 255); 
				blackPrimarySecondaryWhite[5][i] = new Color(255, i-512, 255); 
				
			}
		}
	}
	public static void initNanonisBlueType()
	{
		nanonisBlueType = new Color[6][512];
		for (int i = 0; i < 512; i++)
		{
			if (i < 256)
			{
				nanonisBlueType[0][i] = new Color(i, 0, i/2); 
				nanonisBlueType[1][i] = new Color(i, i/2, 0); 
				nanonisBlueType[2][i] = new Color(0, i, i/2); 
				nanonisBlueType[3][i] = new Color(i/2, i, 0); 
				nanonisBlueType[4][i] = new Color(i/2, 0, i); 
				nanonisBlueType[5][i] = new Color(0, i/2, i); 
			}
			else
			{
				nanonisBlueType[0][i] = new Color(255, i-256, i/2);
				nanonisBlueType[1][i] = new Color(255, i/2, i-256);
				nanonisBlueType[2][i] = new Color(i-256, 255, i/2); 
				nanonisBlueType[3][i] = new Color(i/2, 255, i-256); 
				nanonisBlueType[4][i] = new Color(i/2, i-256, 255); 
				nanonisBlueType[5][i] = new Color(i-256, i/2, 255); 
			}
		}
	}
	public static void initRainbow()
	{
		rainbow = new Color[336];
		int[] r = new int [336];
		int[] g = new int [336];
		int[] b = new int [336];
		for (int i = 0; i < 80; i++)
			r[i] = (int)(128 - i*128.0/80);
		for (int i = 80; i < 160; i++)
			r[i] = 0;
		for (int i = 160; i < 256; i++)
			r[i] = (int)((i-160)*(255.0/96));
		for (int i = 256; i < 336; i++)
			r[i] = 255;

		for (int i = 80; i < 144; i++)
			g[i] = (int)((i-80)*255.0/64);
		for (int i = 144; i < 256; i++)
			g[i] = 255;
		for (int i = 256; i < 336; i++)
			g[i] = (int)(255 - (i-256)*255.0/80);
		
		for (int i = 0; i < 48; i++)
			b[i] = 127 + (int)(i*128.0/48);
		for (int i = 48; i < 144; i++)
			b[i] = 255;
		for (int i = 144; i < 168; i++)
			b[i] = (int)(255 - (i-144)*255.0/24);
		
		for (int i = 0; i < 336; i++){
			rainbow[i] = new Color(r[i], g[i], b[i]);
//			System.out.println(rainbow[i]);
		}
		
	
	}
	
	public static void initColors(int i)
	{
		switch (i){
		case 0: initBRYW();
			break;
		case 1: initGrey();
			break;
		case 2: initBlueWhiteRed();
			break;
		case 3:
		case 4:
		case 5:
			initBlackPrimaryWhite();
			break;
		case 6:
		case 7:
		case 8:
		case 9:
		case 10:
			initBlackPrimarySecondaryWhite();
			break;
		case 11:
			initGreyBounded();
			break;
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
			initNanonisBlueType();
			break;
		case 18:
			initRainbow();
			break;
		case 19:
			initRHK();
			break;
		case 20: 
			initPurpleWhiteOrange();
			break;
//		case 4: initBlackBlueCyanWhite();
//			break;
//		case 5: initBlackGreenYellowWhite();
//			break;
		}
	}
	public static Color[] getColors(int i)
	{
		switch (i){
		case 0:	return bryw;
		case 1: return grey;
		case 2: return blueWhiteRed;
		case 3:
		case 4:
		case 5:
			return blackPrimaryWhite[i-3];
		case 6:
		case 7:
		case 8:
		case 9:
		case 10:
			return blackPrimarySecondaryWhite[i-5];
		case 11: return greyBounded;
		case 12:
		case 13:
		case 14:
		case 15:
		case 16:
		case 17:
			return nanonisBlueType[i-12];
		case 18: return rainbow;
		case 19: return rhkDefault;
		case 20: return purpleWhiteOrange;
		//		case 4: return blackBlueCyanWhite;
//		case 5: return blackGreenYellowWhite;
		}
		return null;
	}
	public static class Greyscale extends ColorScale1D
	{
		public Greyscale(double max, double min) {
			super(grey, min, max);
		}
	}
	
	public static class LinearBRYW extends ColorScale1D
	{
		public LinearBRYW(double max, double min) {
			super(bryw, min, max);
		}
	}
	public static class BlueWhiteRed extends ColorScale1D
	{
		public BlueWhiteRed(double max, double min) {
			super(blueWhiteRed, min, max);
		}
	}
	public static class PurpleWhiteOrange extends ColorScale1D
	{
		public PurpleWhiteOrange(double max, double min) {
			super(purpleWhiteOrange, min, max);
		}
	}
	public static class BlackBlueWhite extends ColorScale1D
	{
		public BlackBlueWhite(double max, double min) {
			super(blackPrimaryWhite[2], min, max);
		}
	}
	public static class BlackRedWhite extends ColorScale1D
	{
		public BlackRedWhite(double max, double min) {
			super(blackPrimaryWhite[0], min, max);
		}
	}
	public static class BlackGreenWhite extends ColorScale1D
	{
		public BlackGreenWhite(double max, double min) {
			super(blackPrimaryWhite[1], min, max);
		}
	}
	public static class BlackRedMagentaWhite extends ColorScale1D
	{
		public BlackRedMagentaWhite(double max, double min) {
			super(blackPrimarySecondaryWhite[1], min, max);
		}
	}
	public static class BlackGreenYellowWhite extends ColorScale1D
	{
		public BlackGreenYellowWhite(double max, double min) {
			super(blackPrimarySecondaryWhite[2], min, max);
		}
	}
	public static class BlackGreenCyanWhite extends ColorScale1D
	{
		public BlackGreenCyanWhite(double max, double min) {
			super(blackPrimarySecondaryWhite[3], min, max);
		}
	}
	public static class BlackBlueCyanWhite extends ColorScale1D
	{
		public BlackBlueCyanWhite(double max, double min) {
			super(blackPrimarySecondaryWhite[4], min, max);
		}
	}
	public static class BlackBlueMagentaWhite extends ColorScale1D
	{
		public BlackBlueMagentaWhite(double max, double min) {
			super(blackPrimarySecondaryWhite[5], min, max);
		}
	}
	
	public static class GreyBounded extends ColorScale1D
	{
		public GreyBounded(double max, double min) {
			super(greyBounded, min, max);
		}
	}
	
	public static class Rainbow extends ColorScale1D{
		public Rainbow(double max, double min){
			super(rainbow, min, max);
		}
	}
	
	//The nanonisBlue color scale is the actual nanonis blue color scale.
	//Nanonis is in no way responsible for the other color scales on the same 
	//pattern. They are labeled Nanonis only for convenience.
	public static class NanonisBlue extends ColorScale1D
	{
		public NanonisBlue(double max, double min) {
			super(nanonisBlueType[5], min, max);
		}
	}
	public static class NanonisBlueAlt extends ColorScale1D
	{
		public NanonisBlueAlt(double max, double min) {
			super(nanonisBlueType[4], min, max);
		}
	}
	public static class NanonisGreen extends ColorScale1D
	{
		public NanonisGreen(double max, double min) {
			super(nanonisBlueType[3], min, max);
		}
	}
	public static class NanonisGreenAlt extends ColorScale1D
	{
		public NanonisGreenAlt(double max, double min) {
			super(nanonisBlueType[2], min, max);
		}
	}
	public static class NanonisRed extends ColorScale1D
	{
		public NanonisRed(double max, double min) {
			super(nanonisBlueType[1], min, max);
		}
	}
	public static class NanonisRedAlt extends ColorScale1D
	{
		public NanonisRedAlt(double max, double min) {
			super(nanonisBlueType[0], min, max);
		}
	}
	public static class RHKDefault extends ColorScale1D
	{
		public RHKDefault(double max, double min) {
			super(rhkDefault, min, max);
		}
	}

	
	
	//We ASSUME that all values are positive
	public static class LogBRYW implements ColorScale
	{
		double max, min;
		double delta;	
		
		public LogBRYW(double max, double min) {
			if (bryw == null) initBRYW();
			this.max = Math.log(max);
			this.min = Math.log(min);
			delta = this.max - this.min;
		}

		public Color of(double x) {
			int c = (int)((Math.log(x) - min)*768/delta);
			c = Math.min(767, c);
			c = Math.max(0, c);
			return bryw[c];
		}
		public void renormalize(double max, double min)
		{
			this.max = Math.log(max);
			this.min = Math.log(min);
			delta = max - min;
		}
		public double getMin()
		{
			return min;
		}
		public double getMax()
		{
			return max;
		}

	}

	public static class MYC2d implements ColorScale2d
	{

		double period, min, max, delta;
		public MYC2d(double max, double min, double period) {
			this.max = max;
			this.min = min;
			this.period = period;
			delta = max - min;
		}
		public void renormalize(double max, double min) {
			this.max = max;
			this.min = min;
			delta = max - min;
		}

		//a is the amplitude, b is the phase
		public Color of(double a, double b) {
			int amp = (int)((a - min)*256/delta);
			amp = Math.min(255, amp);
			amp = Math.max(0, amp);
			
			
			int phase = (int)((768*b)/period);
			while (phase<0) phase += 768;
			
			phase %= 768;
			
			// TODO Auto-generated method stub
			return myc2d[amp][phase];
		}
		public Color of(double[] z) {
			double a = Complex.mag(z);
			double b = Complex.phase(z);
			int amp = (int)((a - min)*256/delta);
			amp = Math.min(255, amp);
			amp = Math.max(0, amp);
			
			
			int phase = (int)((768*b)/period);
			while (phase<0) phase += 768;
			
			phase %= 768;
			
			// TODO Auto-generated method stub
			return myc2d[amp][phase];
		}
		
	}

	public static ColorScale1D getNew(double min, double max, int type)
	{
		if (getColors(type) == null) initColors(type);
		switch(type)
		{
		case 0: return new LinearBRYW(max, min);
		case 1: return new Greyscale(max, min);
		case 2: return new BlueWhiteRed(max, min);
		case 3: return new BlackRedWhite(max, min);
		case 4: return new BlackGreenWhite(max, min);
		case 5: return new BlackBlueWhite(max, min);
		case 6: return new BlackRedMagentaWhite(max, min);
		case 7: return new BlackGreenYellowWhite(max, min);
		case 8: return new BlackGreenCyanWhite(max, min);
		case 9: return new BlackBlueCyanWhite(max, min);
		case 10: return new BlackBlueMagentaWhite(max, min);
		case 11: return new GreyBounded(max, min);
		case 12: return new NanonisRedAlt(max, min);
		case 13: return new NanonisRed(max, min);
		case 14: return new NanonisGreenAlt(max, min);
		case 15: return new NanonisGreen(max, min);
		case 16: return new NanonisBlueAlt(max, min);
		case 17: return new NanonisBlue(max, min);
		case 18: return new Rainbow(max, min);
		case 19: return new RHKDefault(max, min);
		case 20: return new PurpleWhiteOrange(max, min);
		
//		case 5:
		//		case 4: return new BlackBlueCyanWhite(max, min);
//		case 5: return new BlackGreenYellowWhite(max, min);
		}
		return null;
	}
	public static ColorScale1D getNew(double[][] data, int type)
	{
		double max = ArrayOps.max(data);
		double min = ArrayOps.min(data);
		return getNew(min, max, type);
	}
	
	public static ColorScale1D getNew(double[][] data, int type, boolean invert)
	{
		double max = ArrayOps.max(data);
		double min = ArrayOps.min(data);
		if (invert)
			return getNew(max, min, type);
		else
			return getNew(min, max, type);
	}
	//This is a color which can be used for writing
	public static Color getUnusedColor(int i) {
		switch(i)
		{
		case 0: return Color.blue;
		case 1: return Color.green;
		case 2: return Color.black;
		case 3: return Color.green;
		case 4: return Color.red;
		case 5: return Color.red;
		case 6: return Color.green;
		case 7: return Color.blue;
		case 8: return Color.red;
		case 9: return Color.red;
		case 10: return Color.green;
		case 11: return Color.magenta;
		case 12: return Color.green;
		case 13: return Color.blue;
		case 14: return Color.red;
		case 15: return Color.blue;
		case 16: return Color.green;
		case 17: return Color.red;
		case 18: return Color.black;
		case 19: return Color.blue;
		case 20: return Color.cyan;
		}
		
		return null;
	}
	public static Color getUnusedColor2(int i) {
		switch(i)
		{
		case 0: return Color.green;
		case 1: return Color.red;
		case 2: return Color.green;
		case 3: return Color.cyan;
		case 4: return Color.magenta;
		case 5: return Color.yellow;
		case 6: return Color.yellow;
		case 7: return Color.magenta;
		case 8: return Color.magenta;
		case 9: return Color.magenta;
		case 10: return Color.yellow;
		case 11: return Color.yellow;
		case 12: return Color.cyan;
		case 13: return Color.cyan;
		case 14: return Color.magenta;
		case 15: return Color.magenta;
		case 16: return Color.yellow;
		case 17: return Color.yellow;
		case 18: return Color.MAGENTA;
		case 19: return Color.GREEN;
		case 20: return Color.GREEN;
		}
		
		return null;
	}
	
	public static int getIndex(ColorScale scale)
	{
	 if (scale instanceof LinearBRYW) return 0;
	 if (scale instanceof Greyscale) return 1;
	 if (scale instanceof BlueWhiteRed) return 2;
	 if (scale instanceof BlackRedWhite) return 3;
	 if (scale instanceof BlackGreenWhite) return 4;
	 if (scale instanceof BlackBlueWhite) return 5;
	 if (scale instanceof BlackRedMagentaWhite) return 6;
	 if (scale instanceof BlackGreenYellowWhite) return 7;
	 if (scale instanceof BlackGreenCyanWhite) return 8;
	 if (scale instanceof BlackBlueCyanWhite) return 9;
	 if (scale instanceof BlackBlueMagentaWhite) return 10;
	 if (scale instanceof GreyBounded) return 11;
	 if (scale instanceof NanonisRedAlt) return 12;
	 if (scale instanceof NanonisRed) return 13;
	 if (scale instanceof NanonisGreenAlt) return 14;
	 if (scale instanceof NanonisGreen) return 15;
	 if (scale instanceof NanonisBlueAlt) return 16;
	 if (scale instanceof NanonisBlue) return 17;
	 if (scale instanceof Rainbow) return 18;
	 if (scale instanceof RHKDefault) return 19;
	 if (scale instanceof PurpleWhiteOrange) return 20;
	 return -1;
	}
	
	public static Color getUnusedColor(ColorScale c) {
		return getUnusedColor(getIndex(c));
	}
	
	/**
	 * h goes from 0 to 360
	 * c from 0 to 1
	 * v from 0 to 1
	 * output goes from 0 to 255 each.
	 * @param h
	 * @param c
	 * @param v
	 * @return
	 */
	public static int[] getRGBFromHSV(double h, double s, double v)
	{
		double C = v*s;
		double hPrime = h/60;
		double X = C*(1 - Math.abs((hPrime % 2) - 1));
		
		double r1, g1, b1;
		int intHP = (int)hPrime;
		switch(intHP)
		{
			case 0:
				r1 = C; g1 = X; b1 = 0; break;
			case 1:
				r1 = X; g1 = C; b1 = 0; break;
			case 2:
				r1 = 0; g1 = C; b1 = X; break;
			case 3:
				r1 = 0; g1 = X; b1 = C; break;
			case 4:
				r1 = X; g1 = 0; b1 = C; break;
			case 5:
				r1 = C; g1 = 0; b1 = X; break;
			default:
				r1 = 0; g1 = 0; b1 = 0; break;
		}
		
		double m = v - C;
		double r, g, b;
		r = r1 + m; g = g1 + m; b = b1 + m;
		return new int [] {(int)(255*r), (int)(255*g), (int)(255*b)};
	}
	public static Color[] getFromHSVValues(double h1, double h2, double s1, double s2, double v1, double v2, int ncolors)
	{
		double[] h = ArrayOps.generateArrayInclBoth(h1, h2, ncolors);
		double[] s = ArrayOps.generateArrayInclBoth(s1, s2, ncolors);
		double[] v = ArrayOps.generateArrayInclBoth(v1, v2, ncolors);
		
		Color[] ans = new Color[ncolors];
		for (int i = 0; i < ncolors; i++)
		{
			int[] rgb = getRGBFromHSV(h[i], s[i], v[i]);
			ans[i] = new Color(rgb[0], rgb[1], rgb[2]);
		}
		return ans;
	}

}



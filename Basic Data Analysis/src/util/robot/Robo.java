package util.robot;

import java.awt.AWTException;
import java.awt.Robot;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
public class Robo {

	
	public Robot r = null;
	
	public static int[] cs = new int [128];
	static boolean[] shift = new boolean[128];
	//the last index indicates if a shift is necessary.
	
	int dt = 10;
	
	static{
		cs[32] = KeyEvent.VK_SPACE; shift[32] = false;
		cs[33] = KeyEvent.VK_1; shift[33] = true;
		cs[34] = KeyEvent.VK_QUOTE; shift[34] = true;
		cs[35] = KeyEvent.VK_3; shift[35] = true;
		cs[36] = KeyEvent.VK_4; shift[36] = true;
		cs[37] = KeyEvent.VK_5; shift[37] = true;
		cs[38] = KeyEvent.VK_7; shift[38] = true;
		cs[39] = KeyEvent.VK_QUOTE; shift[39] = false;
		cs[40] = KeyEvent.VK_9; shift[40] = true;
		cs[41] = KeyEvent.VK_0; shift[41] = true;
		cs[42] = KeyEvent.VK_8; shift[42] = true;
		cs[43] = KeyEvent.VK_EQUALS; shift[43] = true;
		cs[44] = KeyEvent.VK_COMMA; shift[44] = false;
		cs[45] = KeyEvent.VK_MINUS; shift[45] = false;
		cs[46] = KeyEvent.VK_PERIOD; shift[46] = false;
		cs[47] = KeyEvent.VK_SLASH; shift[47] = false;
		cs[48] = KeyEvent.VK_0; shift[48] = false;
		cs[49] = KeyEvent.VK_1; shift[49] = false;
		cs[50] = KeyEvent.VK_2; shift[50] = false;
		cs[51] = KeyEvent.VK_3; shift[51] = false;
		cs[52] = KeyEvent.VK_4; shift[52] = false;
		cs[53] = KeyEvent.VK_5; shift[53] = false;
		cs[54] = KeyEvent.VK_6; shift[54] = false;
		cs[55] = KeyEvent.VK_7; shift[55] = false;
		cs[56] = KeyEvent.VK_8; shift[56] = false;
		cs[57] = KeyEvent.VK_9; shift[57] = false;
		cs[58] = KeyEvent.VK_SEMICOLON; shift[58] = true;
		cs[59] = KeyEvent.VK_SEMICOLON; shift[59] = false;
		cs[60] = KeyEvent.VK_COMMA; shift[60] = true;
		cs[61] = KeyEvent.VK_EQUALS; shift[61] = false;
		cs[62] = KeyEvent.VK_PERIOD; shift[62] = true;
		cs[63] = KeyEvent.VK_SLASH; shift[63] = true;
		cs[64] = KeyEvent.VK_2; shift[64] = true;
		cs[65] = KeyEvent.VK_A; shift[65] = true;
		cs[66] = KeyEvent.VK_B; shift[66] = true;
		cs[67] = KeyEvent.VK_C; shift[67] = true;
		cs[68] = KeyEvent.VK_D; shift[68] = true;
		cs[69] = KeyEvent.VK_E; shift[69] = true;
		cs[70] = KeyEvent.VK_F; shift[70] = true;
		cs[71] = KeyEvent.VK_G; shift[71] = true;
		cs[72] = KeyEvent.VK_H; shift[72] = true;
		cs[73] = KeyEvent.VK_I; shift[73] = true;
		cs[74] = KeyEvent.VK_J; shift[74] = true;
		cs[75] = KeyEvent.VK_K; shift[75] = true;
		cs[76] = KeyEvent.VK_L; shift[76] = true;
		cs[77] = KeyEvent.VK_M; shift[77] = true;
		cs[78] = KeyEvent.VK_N; shift[78] = true;
		cs[79] = KeyEvent.VK_O; shift[79] = true;
		cs[80] = KeyEvent.VK_P; shift[80] = true;
		cs[81] = KeyEvent.VK_Q; shift[81] = true;
		cs[82] = KeyEvent.VK_R; shift[82] = true;
		cs[83] = KeyEvent.VK_S; shift[83] = true;
		cs[84] = KeyEvent.VK_T; shift[84] = true;
		cs[85] = KeyEvent.VK_U; shift[85] = true;
		cs[86] = KeyEvent.VK_V; shift[86] = true;
		cs[87] = KeyEvent.VK_W; shift[87] = true;
		cs[88] = KeyEvent.VK_X; shift[88] = true;
		cs[89] = KeyEvent.VK_Y; shift[89] = true;
		cs[90] = KeyEvent.VK_Z; shift[90] = true; //Now I know my ABC's next time won't you sing with me.
		cs[91] = KeyEvent.VK_OPEN_BRACKET; shift[91] = false;
		cs[92] = KeyEvent.VK_BACK_SLASH; shift[92] = false;
		cs[93] = KeyEvent.VK_CLOSE_BRACKET; shift[93] = false;
		cs[94] = KeyEvent.VK_6; shift[94] = true;
		cs[95] = KeyEvent.VK_MINUS; shift[95] = true;
		cs[96] = KeyEvent.VK_QUOTE; shift[96] = false; //unknown
		for (int i = 97; i < 97 + 26; i++)
			{cs[i] = cs[65 + (i-97)]; shift[i] = false;}
		cs[123] = KeyEvent.VK_OPEN_BRACKET; shift[123] = true;
		cs[124] = KeyEvent.VK_BACK_SLASH; shift[124] = true;
		cs[125] = KeyEvent.VK_CLOSE_BRACKET; shift[125] = true;
		cs[126] = cs[96]; shift[126] = true;
		cs[127] = KeyEvent.VK_SPACE; shift[127] = false;
		
		
	}
	
	public Robo()
	{
		try {
			r = new Robot();
		} catch (AWTException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public void typeChar(char c)
	{
		int i = (int)c;
//		System.out.println(i);
		if (shift[i])
		{
			wait(dt);
			r.keyPress(KeyEvent.VK_SHIFT);
			wait(dt);
			r.keyPress(cs[i]);
			wait(dt);
			r.keyRelease(cs[i]);
			wait(dt);
			r.keyRelease(KeyEvent.VK_SHIFT);
			wait(dt);
		}
		else
		{
			wait(dt);
			r.keyPress(cs[i]);
			wait(dt);
			r.keyRelease(cs[i]);
			wait(dt);
		}

	}
	public void typeChar(char c, boolean ctrl)
	{
		int i = (int)c;
//		System.out.println(i);
		if (ctrl)
		{
			wait(dt);
			r.keyPress(KeyEvent.VK_CONTROL);
		}
		if (shift[i])
		{
			wait(dt);
			r.keyPress(KeyEvent.VK_SHIFT);
		}
		wait(dt);
		r.keyPress(cs[i]);
		wait(dt);
		r.keyRelease(cs[i]);
		wait(dt);
		if (shift[i])
		{
			r.keyRelease(KeyEvent.VK_SHIFT);
			wait(dt);
		}
		if (ctrl)
		{
			r.keyRelease(KeyEvent.VK_CONTROL);
			wait(dt);
		}
	}
	
	public void typeString(String s)
	{
		int l = s.length();
		for (int i = 0; i < l; i++)
		{
			typeChar(s.charAt(i));
		}
	}
	public void typeString(String s, int extrawait)
	{
		int l = s.length();
		for (int i = 0; i < l; i++)
		{
			typeChar(s.charAt(i));
			wait(extrawait);
		}
	}
	public void keyType(int key)
	{
		wait(dt);
		r.keyPress(key);
		wait(dt);
		r.keyRelease(key);
		wait(dt);
	}
	public void pressEnter()
	{
		wait(dt);	
		r.keyPress(KeyEvent.VK_ENTER);
		wait(dt);
		r.keyRelease(KeyEvent.VK_ENTER);
		wait(dt);

	}
	
	static int indexOf(int digit)
	{
		return digit + 48;
	}
	
	public void goTo(int[] spot)
	{
		wait(dt);
		r.mouseMove(spot[0], spot[1]);
		wait(dt);
	}
	public void flickMouse()
	{
		wait(dt);
		r.mouseMove(0, 0);
		wait(dt);
		r.mouseMove(100, 0);
		wait(dt);
		r.mouseMove(100, 100);
		wait(dt);
		r.mouseMove(0, 100);
		wait(dt);
	}
	
	public void click()
	{
		wait(dt);
		r.mousePress(InputEvent.BUTTON1_MASK);
		wait(dt);
		r.mouseRelease(InputEvent.BUTTON1_MASK);
		wait(dt);
	}
	
	public void doubleclick()
	{
		wait(dt);
		click();
		click();
		wait(dt);
	}
	
	public void keyTypeCtrl(int key)
	{
		wait(dt);
		r.keyPress(KeyEvent.VK_CONTROL);
		wait(dt);
		r.keyPress(key);
		wait(dt);
		r.keyRelease(key);
		wait(dt);
		r.keyRelease(KeyEvent.VK_CONTROL);
		wait(dt);
	}
	public void keyMaskType(int outer, int inner)
	{
		wait(dt);
		r.keyPress(outer);
		wait(dt);
		r.keyPress(inner);
		wait(dt);
		r.keyRelease(inner);
		wait(dt);
		r.keyRelease(outer);
		wait(dt);
	}
	
	public static void wait(int milli)
	{
		try {
			Thread.sleep(milli);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args)
	{
		char[] ch = new char[1000];
		for (int i = 0; i < 200; i++)
		{
			ch[i] = (char)i;
			System.out.println(i + "\t" + ch[i]);
		}
		wait(3000);
		Robo r = new Robo();
		r.typeString("EasyBMPtoAVI -filebase psi -start 0 -end 100 -framerate 20 -output out.avi");
		r.pressEnter();
		
	}
	
}

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "calplot.h"

float border;
float title;
int Button_Color;
int Button_Color_Light;
int Button_Color_Dark;
int Button_Color_Fore;
int Button_Color_Back;
int black;
int mycursor = 0;

#define	CURSOR_ARROW	0
#define	CURSOR_XORARROW	1
#define	CURSOR_CROSS	2
#define	CURSOR_PLUS	3
#define	CURSOR_BOX	4
#define	CURSOR_RUBBER	5
#define	CURSOR_OFF	6
void setcursor( int);
void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
		float *xup, float *yup);
void dotitle(void);
void dobox(void);
void dorect(void);
void dotri(void);
void doline(void);
void doclear(void);
static int inside(float xv, float yv, 
	float xlb, float ylb, float xhb, float yhb);

struct menu { float xl; float yl; float xh 
	; float yh ; char *str ; int action;  } ;
static struct menu m[] = { 
	{  -1.0, -1.0, -1.0, -1.0, " Line \0" , 1},
	{  -1.0, -1.0, -1.0, -1.0, " Tri  \0" , 2},
	{  -1.0, -1.0, -1.0, -1.0, " Rect \0" , 3},
	{  -1.0, -1.0, -1.0, -1.0, " Back \0" , 4},
	{  -1.0, -1.0, -1.0, -1.0, " Fore \0" , 5},
	{  -1.0, -1.0, -1.0, -1.0, " Clear\0" , 6},
	{  -1.0, -1.0, -1.0, -1.0, " Red  \0" , 7},
	{  -1.0, -1.0, -1.0, -1.0, "Orange\0" , 8},
	{  -1.0, -1.0, -1.0, -1.0, " Green\0" , 9},
	{  -1.0, -1.0, -1.0, -1.0, " Cyan \0" ,10},
	{  -1.0, -1.0, -1.0, -1.0, " Blue \0" ,11},
	{  -1.0, -1.0, -1.0, -1.0, " XOR  \0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " Arrow\0" ,13},
	{  -1.0, -1.0, -1.0, -1.0, "XorCro\0" ,14},
	{  -1.0, -1.0, -1.0, -1.0, "XorArr\0" ,15},
	{  -1.0, -1.0, -1.0, -1.0, "XorPls\0" ,16},
	{  -1.0, -1.0, -1.0, -1.0, " Gread\0" ,17},
	{  -1.0, -1.0, -1.0, -1.0, " Exit \0" ,18}
} ;

/*
	{  -1.0, -1.0, -1.0, -1.0, "SYMBOL\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, "NUMBER\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, "FACTOR\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, "CIRCLE\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " SFILL\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, "CENTER\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " LEFT \0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " RIGHT\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " WIDTH\0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " FONT \0" ,12},
	{  -1.0, -1.0, -1.0, -1.0, " ARROW\0" ,12},
*/
#define	ON	1
#define OFF	0
#define XL	0.25
#define	YL	1.0
#define	XH	8.25
#define	YH	7.4
#define        MIN(a,b) ( (b) > (a) ? (a):(b) )
#define MAX(a,b) ( (a) > (b) ? (a):(b) )

/* globals */
int kolor = 1;
static int Xoring = OFF;
static int butrev = OFF ;



int main(void)
{
	float xl,yl,xh,yh;
	int nm, i;
	int ii;
	int cmd;
	float xv, yv;
	char c[2];

	ginitf("INTER","CDRAW");
	border = 8.0/60.0;
	kolor = 1;
	black = 1;
	

	nm = sizeof(m)/sizeof(struct menu);
	for(i=0 ; i < nm ; i++){
		draw_button(8.5,7.0 -i*0.4,m[i].str,&xl,&yl,&xh,&yh);
		m[i].xl = xl;
		m[i].yl = yl;
		m[i].xh = xh;
		m[i].yh = yh;
	}
	/* loop for values */
	cmd = -1;
	dobox();
	gcursor("arrow");
	for(; ;){
		curaxy(&xv, &yv, c);
		cmd = -1;
		for(i=0 ; i < nm ; i++){
			if(inside(xv,yv,m[i].xl,m[i].yl,m[i].xh,m[i].yh))
			{
				cmd = m[i].action;
				gmesg(m[i].str);
				ii = i;
				break;
			}
		}
		if(cmd > 0){
		switch(cmd){
			case 1:
				butrev = ON;
		draw_button(8.5,7.0 -ii*0.4,m[ii].str,&xl,&yl,&xh,&yh);
				doline();
				break;
			case 2:
				butrev = ON;
		draw_button(8.5,7.0 -ii*0.4,m[ii].str,&xl,&yl,&xh,&yh);
				dotri();
				break;
			case 3:
				butrev = ON;
		draw_button(8.5,7.0 -ii*0.4,m[ii].str,&xl,&yl,&xh,&yh);
				dorect();
				break;
			case 4:
				kolor = 0;
				break;
			case 5:
				kolor = 1;
				break;
			case 6:
				doclear();
				newpen(1);
				kolor = 1;
				break;
			case 7:
				kolor = 1000;
				break;
			case 8:
				kolor = 1025;
				break;
			case 9:
				kolor = 1050;
				break;
			case 10:
				kolor = 1075;
				break;
			case 11:
				kolor = 1100;
				break;
			case 12:
				if(Xoring){
					Xoring = OFF;
				} else {
					Xoring = ON;
				}
				break;
			case 13:
				mycursor = CURSOR_ARROW;
				setcursor(mycursor);
				break;
			case 14:
				mycursor = CURSOR_CROSS;
				setcursor(mycursor);
				break;
			case 15:
				mycursor = CURSOR_XORARROW;
				setcursor(mycursor);
				break;
			case 16:
				mycursor = CURSOR_PLUS;
				setcursor(mycursor);
				break;
			case 17:
				newpen(1);
				gread("PLTTST.PLT",XL,YL,XL-XL,YL-YL,
				XH-XL,YH-YL,2,1.0,1.0);
				newpen(kolor);
				break;
			case 18:
				goto end;
			default:
				break;
		}
		butrev = OFF ;
		draw_button(8.5,7.0 -ii*0.4,m[ii].str,&xl,&yl,&xh,&yh);
		dotitle();
		}
	}
end:
	gmesg("Closing Session - Any Key to End");
	pend();
	return 0;
}

static int inside(float xv, float yv, float xlb, 
	float ylb, float xhb, float yhb)
{
	if(xv >= xlb && xv <= xhb && yv >= ylb && yv <= yhb)
		return 1;
	else
		return 0;
}

void draw_button(float xl, float yl, char *str, float *xlw, float *ylw,
		float *xup, float *yup)
{
	int lstr;
	int nc;
	float bdr = 0.04;
	float title = 20.0/60.0;
	float x[6], y[6];
	lstr = (int)strlen(str);
	/* add a space */
	nc = lstr*8;
	gclip("off", XL, YL, XH, YH);
		if(butrev == OFF){
			Button_Color_Dark  = 1100;
			Button_Color       = 1081;
			Button_Color_Light = 1071;
		} else {
			Button_Color_Dark  = 1071;
			Button_Color       = 1081;
			Button_Color_Light = 1101;
		}
		/* color */
		if(black){
			Button_Color_Fore  =    0;
			Button_Color_Back  =    1;
		} else {
			Button_Color_Fore  =    1;
			Button_Color_Back  =    0;
		}

	newpen(Button_Color_Back);
	x[0] = xl ;
	y[0] = yl ;
	x[1] = xl + lstr*0.125 ;
	y[1] = yl ;
	x[2] = xl + lstr*0.125 ;
	y[2] = yl + title ;
	x[3] = xl ;
	y[3] = yl + title ;
	shadep(4,x,y);

	newpen(Button_Color_Dark);
	x[0] = xl;
	y[0] = yl;
	x[1] = xl + lstr*0.125;
	y[1] = yl;
	x[2] = xl + lstr*0.125;
	y[2] = yl + title ;
	x[3] = xl + lstr*0.125 - bdr;
	y[3] = yl + title - bdr ;
	x[4] = xl + lstr*0.125 - bdr;
	y[4] = yl + bdr ;
	x[5] = xl + bdr;
	y[5] = yl + bdr ;
	shadep(6,x,y);

	newpen(Button_Color_Light);
	x[0] = xl;
	y[0] = yl;
	x[1] = xl + bdr;
	y[1] = yl + bdr ;
	x[2] = xl + bdr;
	y[2] = yl + title - bdr ;
	x[3] = xl + lstr*0.125 - bdr;
	y[3] = yl + title - bdr ;
	x[4] = xl + lstr*0.125 ;
	y[4] = yl + title  ;
	x[5] = xl  ;
	y[5] = yl + title  ;
	shadep(6,x,y);

	newpen(Button_Color);
	x[0] = xl + bdr;
	y[0] = yl + bdr;
	x[1] = xl + lstr*0.125 - bdr;
	y[1] = yl + bdr;
	x[2] = xl + lstr*0.125 - bdr;
	y[2] = yl + title - bdr;
	x[3] = xl + bdr;
	y[3] = yl + title - bdr;
	shadep(4,x,y);
	newpen(Button_Color_Fore);
	gwrtxt(xl+0.0700,yl+0.090,str,0);
	*xlw = xl;
	*ylw = yl;
	*xup = xl + lstr*0.125;
	*yup = yl + title;
	gclip("on", XL, YL, XH, YH);
}


void doclear(void)
{
	newpen(0);
	shader(XL,YL,XH,YH, 0, 0, 0.01, 0.01);
	dobox();
}

void doline(void)
{
	char c[2];
	float xv1, xv2, xv3, yv1, yv2, yv3;
	curaxy(&xv1, &yv1, c);
	gcursor("Rubber");
	curaxy(&xv2, &yv2, c);
	setcursor(mycursor);
	if( inside(xv1,yv1,XL,YL,XH,YH) && inside(xv2,yv2,XL,YL,XH,YH)){
		/* this sequence is required to make the XOR work
			on the PC */
		if(Xoring==ON){
			newpen(3001);
		}
		newpen(kolor);
		plot(xv1, yv1, 3);
		plot(xv2, yv2, 2);
		if(Xoring==ON){
			newpen(3000);
		}
	}
}

void dotri(void)
{
	char c[2];
	float xv1, xv2, xv3, yv1, yv2, yv3;
	curaxy(&xv1, &yv1, c);
	gcursor("Rubber");
	curaxy(&xv2, &yv2, c);
	curaxy(&xv3, &yv3, c);
	setcursor(mycursor);
	if( inside(xv1,yv1,XL,YL,XH,YH) 
			&& inside(xv3,yv3,XL,YL,XH,YH)
			&& inside(xv2,yv2,XL,YL,XH,YH)){
		if(Xoring==ON){
			newpen(3001);
		}
		newpen(kolor);
		shadet(xv1,yv1,xv2,yv2,xv3,yv3, 0, 0, 0.01, 0.01);
		if(Xoring ==ON){
			newpen(3000);
		}
	}
}

void dorect(void)
{
	char c[2];
	float xv1, xv2, xv3, yv1, yv2, yv3;
	curaxy(&xv1, &yv1, c);
	gcursor("Box");
	curaxy(&xv2, &yv2, c);
	setcursor(mycursor);
	if( inside(xv1,yv1,XL,YL,XH,YH) && inside(xv1,yv1,XL,YL,XH,YH)){
		if(Xoring ==ON){
			newpen(3001);
		}
		newpen(kolor);
		shader(xv1,yv1,xv2,yv2, 0, 0, 0.01, 0.01);
		if(Xoring ==ON){
			newpen(3000);
		}
	}
}

void dobox(void)
{
	newpen(1);
	dotitle();
	plot(XL,YL,3);
	plot(XH,YL,2);
	plot(XH,YH,2);
	plot(XL,YH,2);
	plot(XL,YL,2);
	plot(XL,YL,3);
	gclip("on", XL, YL, XH, YH);
}

void dotitle(void)
{
	char strout[11];
	gclip("off", XL, YL, XH, YH);
	newpen(2);
	sprintf(strout,"COLOR=%4d",kolor);
	gwrtxt(1.0,7.5,strout,2);
	sprintf(strout,"XOR  =%4d",Xoring);
	gwrtxt(3.0,7.5,strout,2);
	gclip("on", XL, YL, XH, YH);
	newpen(kolor);
}

void setcursor( int curs)
{
	if(curs <= CURSOR_ARROW)
		gcursor("arrow"); 
	else if(curs == CURSOR_XORARROW)
		gcursor("xorarrow"); 
	else if(curs == CURSOR_CROSS)
		gcursor("cross"); 
	else if(curs == CURSOR_PLUS)
		gcursor("plus"); 
	else if(curs == CURSOR_BOX )
		gcursor("box"); 
	else if(curs == CURSOR_RUBBER)
		gcursor("rubber"); 
	else if(curs == CURSOR_OFF)
		gcursor("off"); 
	else
		gcursor("arrow");
}


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <windows.h>
#include <glut.h>

#include "header.hpp"
#include "bitmap_image.hpp"

using namespace std;

#define pi (2*acos(0.0))

extern vector<Object*> objects;
extern vector<Point3> lights;


Point3 pos(0,-100,10);
Point3 l(0,1,0);
Point3 r(1,0,0);
Point3 u(0,0,1);

double rx;
double ry;
double rz;
double n1;
double n2;

double c = cos(.01);
double s = sin(.01);

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes = 1;
double angle;

int windowWidth = 500;
int windowHeight = 500;

int imageHeight, imageWidth;
double fovy = 80;

double thickness = 0.5;

void drawPoint(Point3 point){

    glColor3f(1.0,1.0,0);

    glBegin(GL_QUADS);
    {
        glVertex3f(point.x + thickness, point.y, point.z + thickness);
        glVertex3f(point.x - thickness, point.y, point.z + thickness);
        glVertex3f(point.x - thickness, point.y, point.z - thickness);
        glVertex3f(point.x + thickness, point.y, point.z - thickness);
    }
    glEnd();

}

void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}



void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}

void Capture(){

    ///debug
    for(int i=0; i<objects.size(); i++){
        cout<<objects[i]->color[0]<<"  "<<objects[i]->color[1]<<"  "<<objects[i]->color[2]<<endl;
	}

    bitmap_image image(imageWidth,imageWidth);

    for(int i=0;i<imageWidth;i++){
        for(int j=0;j<imageWidth;j++){
            image.set_pixel(i,j,0,0,0);
        }
    }

    double plane_distance = (windowHeight/2)/tan(fovy*pi/360);

    Point3 *topleft = new Point3();
    topleft->x = pos.x + l.x * plane_distance - r.x*windowWidth/2 + u.x*windowHeight/2;
    topleft->y = pos.y + l.y * plane_distance - r.y*windowWidth/2 + u.y*windowHeight/2;
    topleft->z = pos.z + l.z * plane_distance - r.z*windowWidth/2 + u.z*windowHeight/2;

    double du = (windowWidth*1.0) / imageWidth;
    double dv = (windowHeight*1.0) / imageHeight;

    for (int i=0; i<imageWidth; i++){
        //cout<<i<<endl;
        for(int j=0;j<imageWidth;j++){
            //cout<<"i j "<<i<<" "<<j<<endl;
            int nearest = -1;
            double t_min = 9999999;

            Point3 *corner = new Point3();
            corner->x = topleft->x + r.x*j*du - u.x*i*dv;
            corner->y = topleft->y + r.y*j*du - u.y*i*dv;
            corner->z = topleft->z + r.z*j*du - u.z*i*dv;

            Point3 dir(corner->x-pos.x, corner->y-pos.y, corner->z-pos.z);
            Ray ray(pos, dir);

            double color[3];

            for(int k=0; k<objects.size();k++){
                double t = objects[k]->intersect(&ray,color,0);
                if(t<=0)continue;
                else if(t<t_min){
                    t_min = t;
                    nearest = k;
                }
            }
            if(nearest!=-1){
                double t = objects[nearest]->intersect(&ray,color,1);
                //cout<<objects[k]->color[0]<<"  "<<objects[k]->color[1]<<"  "<<objects[k]->color[2]<<endl;
                image.set_pixel(j,i,color[0],color[1],color[2]);
                //image.set_pixel(j,i,255,0,0);
            }
        }
    }
    image.save_image("output.bmp");
    exit(0);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
            rx = r.x*c + l.x*s;
            ry = r.y*c + l.y*s;
            rz = r.z*c + l.z*s;
            n1 = sqrt(rx*rx + ry*ry + rz*rz);

		    l.x = l.x*c - rx*s;
		    l.y = l.y*c - ry*s;
		    l.z = l.z*c - rz*s;
            n2 = sqrt(l.x*l.x + l.y*l.y + l.z*l.z);

            r.x = rx/n1;
		    r.y = ry/n1;
		    r.z = rz/n1;

            l.x = l.x/n2;
		    l.y = l.y/n2;
		    l.z = l.z/n2;
			break;

		case '2':
            rx = r.x*c - l.x*s;
            ry = r.y*c - l.y*s;
            rz = r.z*c - l.z*s;
            n1 = sqrt(rx*rx + ry*ry + rz*rz);

		    l.x = l.x*c + rx*s;
		    l.y = l.y*c +  ry*s;
		    l.z = l.z*c + rz*s;
            n2 = sqrt(l.x*l.x + l.y*l.y + l.z*l.z);

            r.x = rx/n1;
		    r.y = ry/n1;
		    r.z = rz/n1;

            l.x = l.x/n2;
		    l.y = l.y/n2;
		    l.z = l.z/n2;
			break;

		case '3':

            rx = u.x*c - l.x*s;
		    ry = u.y*c - l.y*s;
		    rz = u.z*c - l.z*s;
            n1 = sqrt(rx*rx + ry*ry + rz*rz);

		    l.x = l.x*c + rx*s;
		    l.y = l.y*c + ry*s;
		    l.z = l.z*c + rz*s;
            n2 = sqrt(l.x*l.x + l.y*l.y + l.z*l.z);

            u.x = rx/n1;
		    u.y = ry/n1;
		    u.z = rz/n1;

            l.x = l.x/n2;
		    l.y = l.y/n2;
		    l.z = l.z/n2;
			break;

		case '4':

            rx = u.x*c + l.x*s;
		    ry = u.y*c + l.y*s;
		    rz = u.z*c + l.z*s;
            n1 = sqrt(rx*rx + ry*ry + rz*rz);


		    l.x = l.x*c - rx*s;
		    l.y = l.y*c - ry*s;
		    l.z = l.z*c - rz*s;
            n2 = sqrt(l.x*l.x + l.y*l.y + l.z*l.z);

            u.x = rx/n1;
		    u.y = ry/n1;
		    u.z = rz/n1;

            l.x = l.x/n2;
		    l.y = l.y/n2;
		    l.z = l.z/n2;
			break;

		case '5':
            rx = u.x*c + r.x*s;
		    ry = u.y*c + r.y*s;
		    rz = u.z*c + r.z*s;
            n1 = sqrt(rx*rx + ry*ry + rz*rz);


		    r.x = r.x*c - rx*s;
		    r.y = r.y*c - ry*s;
		    r.z = r.z*c - rz*s;
            n2 = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);

            u.x = rx/n1;
		    u.y = ry/n1;
		    u.z = rz/n1;

            r.x = r.x/n2;
		    r.y = r.y/n2;
		    r.z = r.z/n2;
			break;

		case '6':
            rx = u.x*c - r.x*s;
		    ry = u.y*c - r.y*s;
		    rz = u.z*c - r.z*s;
            n1 = sqrt(rx*rx + ry*ry + rz*rz);

		    r.x = r.x*c + rx*s;
		    r.y = r.y*c + ry*s;
		    r.z = r.z*c + rz*s;
            n2 = sqrt(r.x*r.x + r.y*r.y + r.z*r.z);

            u.x = rx/n1;
		    u.y = ry/n1;
		    u.z = rz/n1;

            r.x = r.x/n2;
		    r.y = r.y/n2;
		    r.z = r.z/n2;
			break;

        case '0':
            Capture();
            break;

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_UP:		//up arrow key
			pos.x += 2*l.x;
			pos.y += 2*l.y;
			pos.z += 2*l.z;
			break;
		case GLUT_KEY_DOWN:		// down arrow key
			pos.x -= 2*l.x;
			pos.y -= 2*l.y;
			pos.z -= 2*l.z;
			break;

		case GLUT_KEY_RIGHT:
            pos.x += 2*r.x;
			pos.y += 2*r.y;
			pos.z += 2*r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x -= 2*r.x;
			pos.y -= 2*r.y;
			pos.z -= 2*r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += 2*u.x;
			pos.y += 2*u.y;
			pos.z += 2*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x -= 2*u.x;
			pos.y -= 2*u.y;
			pos.z -= 2*u.z;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			break;
		case GLUT_KEY_END:
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
			}
			break;

		case GLUT_RIGHT_BUTTON:
			//........
			break;

		case GLUT_MIDDLE_BUTTON:
			//........
			break;

		default:
			break;
	}
}





void display(){

	//clear the display
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(0,0,0,0);	//color black
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	/********************
	/ set-up camera here
	********************/
	//load the correct matrix -- MODEL-VIEW matrix
	glMatrixMode(GL_MODELVIEW);

	//initialize the matrix
	glLoadIdentity();

    //gluLookAt(0,0,200,	0,0,0,	0,1,0);
    gluLookAt(pos.x, pos.y, pos.z, pos.x + l.x, pos.y + l.y, pos.z + l.z, u.x, u.y, u.z);

    //drawSquare(50);

	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	//drawAxes();

	for(int i=0; i<objects.size(); i++){
        objects[i]->draw();
	}
    for(int i=0; i<lights.size(); i++){
        drawPoint(lights[i]);
	}


    //drawPoint(Point3(20,10,30),.5);

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}


void init(){
	//codes for initialization
	//clear the screen
	glClearColor(0,0,0,0);
	/************************
	/ set-up projection here
	************************/
	//load the PROJECTION matrix
	glMatrixMode(GL_PROJECTION);
	//initialize the matrix
	glLoadIdentity();
	//give PERSPECTIVE parameters
	gluPerspective(fovy,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

void loadTestData(){

    Object *temp;

    Point3 center(0,0,10);
    temp = new Sphere(center,10);
    temp->setColor(255,0,0);
    temp->setCoefficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);

    objects.push_back(temp);

    Point3 center2(15,20,20);
    temp = new Sphere(center2,10);
    temp->setColor(0,255,0);
    temp->setCoefficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);

    objects.push_back(temp);

    Point3 center3(20,0,100);
    temp = new Sphere(center3,10);
    temp->setColor(0,255,255);
    temp->setCoefficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);

    objects.push_back(temp);

    Point3 center4(25,50,10);
    temp = new Sphere(center4,10);
    temp->setColor(255,255,0);
    temp->setCoefficients(0.4,0.2,0.2,0.2);
    temp->setShine(1);

    objects.push_back(temp);

    Point3 light1(-50,50,50);
    lights.push_back(light1);

//    temp=new Floor(1000, 20);
//    temp->setCoefficients(0.4,0.2,0.2,0.2);
//    temp->setShine(1);
//    objects.push_back(temp);

    imageWidth = imageHeight = 768;

}



void test(){
    Point3 *topleft = new Point3();
    topleft->x = 1;
    topleft->y = 2;
    topleft->z = 3;

    cout<<topleft->x;
}

int main(int argc, char **argv){

    test();
    loadTestData();
    cout<<imageWidth<<endl;

	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}

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
extern vector<Point> lights;


Point pos(0,-100,10);
Point l(0,1,0);
Point r(1,0,0);
Point u(0,0,1);

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

int recursion_level;

int imageHeight, imageWidth;
double fovy = 80;

double thickness = 0.5;

void drawPoint(Point point){

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

    Point *topleft = new Point();
    topleft->x = pos.x + l.x * plane_distance - r.x*windowWidth/2 + u.x*windowHeight/2;
    topleft->y = pos.y + l.y * plane_distance - r.y*windowWidth/2 + u.y*windowHeight/2;
    topleft->z = pos.z + l.z * plane_distance - r.z*windowWidth/2 + u.z*windowHeight/2;

    double du = (windowWidth*1.0) / imageWidth;
    double dv = (windowHeight*1.0) / imageHeight;

    for (int i=0; i<imageWidth; i++){
        for(int j=0;j<imageWidth;j++){
            int nearest = -1;
            double t_min = 9999999;

            Point *corner = new Point();
            corner->x = topleft->x + r.x*j*du - u.x*i*dv;
            corner->y = topleft->y + r.y*j*du - u.y*i*dv;
            corner->z = topleft->z + r.z*j*du - u.z*i*dv;

            double dirx = corner->x - pos.x;
            double diry = corner->y - pos.y;
            double dirz = corner->z - pos.z;
            Point dir(dirx, diry, dirz);

            Ray ray(pos, dir);

            double dummyColorAt[3];

            for(int k=0; k<objects.size();k++){

                double t = objects[k]->intersect(&ray,dummyColorAt,0);

                if(t<=0)continue;

                if(t<t_min){
                    t_min = t;
                    nearest = k;
                }
            }

            //double colorAt[3];

            if(nearest!=-1){
                objects[nearest]->intersect(&ray,dummyColorAt,1);
                image.set_pixel(j,i,dummyColorAt[0],dummyColorAt[1],dummyColorAt[2]);
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
			pos.x += 10*l.x;
			pos.y += 10*l.y;
			pos.z += 10*l.z;
			break;
		case GLUT_KEY_DOWN:		// down arrow key
			pos.x -= 10*l.x;
			pos.y -= 10*l.y;
			pos.z -= 10*l.z;
			break;

		case GLUT_KEY_RIGHT:
            pos.x += 10*r.x;
			pos.y += 10*r.y;
			pos.z += 10*r.z;
			break;
		case GLUT_KEY_LEFT:
			pos.x -= 10*r.x;
			pos.y -= 10*r.y;
			pos.z -= 10*r.z;
			break;

		case GLUT_KEY_PAGE_UP:
		    pos.x += 10*u.x;
			pos.y += 10*u.y;
			pos.z += 10*u.z;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos.x -= 10*u.x;
			pos.y -= 10*u.y;
			pos.z -= 10*u.z;
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
	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
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

void animate() {
    glutPostRedisplay();
}

void loadTestData(){
    imageWidth = 768;
    imageHeight = 768;
    recursion_level = 3;

    Object *temp;

//    Point center(40,0,10);
//    temp = new Sphere(center,10);
//    temp->setColor(1,0,0);
//    temp->setCoefficients(0.4,0.2,0.2,0.2);
//    temp->setShine(10);
//    objects.push_back(temp);
//
//    Point center2(-30,60,20);
//    temp = new Sphere(center2,20);
//    temp->setColor(0,0,1);
//    temp->setCoefficients(0.2,0.2,0.4,0.2);
//    temp->setShine(15);
//    objects.push_back(temp);

//    Point center3(-15.0, 15.0, 45.0);
//    temp = new Sphere(center3,15.0);
//    temp->setColor(1, 1, 0);
//    temp->setCoefficients(0.4,0.3,0.1,0.2);
//    temp->setShine(5);
//    objects.push_back(temp);

    Point light1(70,-10,70);
    lights.push_back(light1);

    Point light2(-70,-10,70);
    lights.push_back(light2);

    temp=new Floor(1000, 20);
    temp->setCoefficients(0.4,0.2,0.2,0.2);
    temp->setShine(5);
    objects.push_back(temp);
Point q[3];
    ///triangles
//    Point a1(-60, 0, 0);
//    Point b1(-30, 0, 0);
//    Point c1(-45, 15, 15);
//    q[0] = a1;q[1] = b1;q[2] = c1;
//    temp = new Triangle(q);
//    temp->setColor(1, 0, 0);
//    temp->setCoefficients(0.4, 0.2, 0.1, 0.3);
//    temp->setShine(5);
//    objects.push_back(temp);
//
//    Point a2(0, 0, 0);
//    Point b2(30, 0, 0);
//    Point c2(15, -15, 15);
//    q[0] = a2;q[1] = b2;q[2] = c2;
//    temp = new Triangle(q);
//    temp->setColor(1, 0, 0);
//    temp->setCoefficients(0.4, 0.2, 0.1, 0.3);
//    temp->setShine(5);
//    objects.push_back(temp);
//
//    Point a3(30, 60, 0);
//    Point b3(50, 30, 0);
//    Point c3(50, 45, 50);
//    q[0] = a3;q[1] = b3;q[2] = c3;
//    temp = new Triangle(q);
//    temp->setColor(1, 0, 0);
//    temp->setCoefficients(0.4, 0.2, 0.1, 0.3);
//    temp->setShine(5);
//    objects.push_back(temp);
//
//    temp = new genQuad(0, 0, 0 , 0, 0, 20);
//    temp->setVar(1, 1, 1, 0, 0, 0, 0, 0, 0, 100);
//    temp->setColor(0, 1, 0);
//    temp->setCoefficients(0.4, 0.2, 0.1, 0.3);
//    temp->setShine(10);
//    objects.push_back(temp);
//
//    temp = new genQuad(0, 0, 0 , 0, 0, 15);
//    temp->setVar(0.0625, 0.04, 0.04, 0, 0, 0, 0, 0, 0, -36);
//    temp->setColor(1, 0, 0);
//    temp->setCoefficients(0.4, 0.2, 0.1, 0.3);
//    temp->setShine(15);
//    objects.push_back(temp);
}

void LoadActualData(){
    ifstream file;
    file.open("scene.txt");
    Object *temp;
    string type;

    temp=new Floor(1000, 20);
    temp->setCoefficients(0.4,0.2,0.2,0.2);
    temp->setShine(5);
    objects.push_back(temp);

    file>>recursion_level;
    cout<<recursion_level<<endl;
    file>>imageWidth;
    imageHeight = imageWidth;
    cout<<imageWidth<<endl;
    int numberOfObjects, numberOfLights;
    file>>numberOfObjects;
    cout<<numberOfObjects<<endl;

    for(int i=0; i<numberOfObjects; i++){
        file>>type;
        cout<<i<<". "<<type<<endl;

        if(!type.compare("sphere")){
            double x,y,z,a,b,c,d,radius,c1,c2,c3,shine;

            file>>x>>y>>z;
            file>>radius;
            file>>c1>>c2>>c3;
            file>>a>>b>>c>>d;
            file>>shine;

            cout<<" "<<x<<" "<<y<<" "<<z<<" "<<endl;
            cout<<" "<<radius<<" "<<endl;
            cout<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<endl;
            cout<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<endl;
            cout<<" "<<shine<<" "<<endl;

            Point center(x,y,z);
            temp = new Sphere(center,radius);
            temp->setColor(c1,c2,c3);
            temp->setCoefficients(a,b,c,d);
            temp->setShine(shine);
            objects.push_back(temp);
        }

        else if(!type.compare("triangle")){
            double x1,y1,z1,x2,y2,z2,x3,y3,z3,c1,c2,c3,a,b,c,d,shine;

            file>>x1>>y1>>z1;
            file>>x2>>y2>>z2;
            file>>x3>>y3>>z3;
            file>>c1>>c2>>c3;
            file>>a>>b>>c>>d;
            file>>shine;

            cout<<" "<<x1<<" "<<y1<<" "<<z1<<" "<<endl;
            cout<<" "<<x2<<" "<<y2<<" "<<z2<<" "<<endl;
            cout<<" "<<x3<<" "<<y3<<" "<<z3<<" "<<endl;
            cout<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<endl;
            cout<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<endl;
            cout<<" "<<shine<<" "<<endl;

            Point a2(x1,y1,z1);
            Point b2(x2,y2,z2);
            Point cc2(x3,y3,z3);
            Point q[3];q[0] = a2;q[1] = b2;q[2] = cc2;
            temp = new Triangle(q);
            temp->setColor(c1,c2,c3);
            temp->setCoefficients(a,b,c,d);
            temp->setShine(shine);
            objects.push_back(temp);
        }

        else if(!type.compare("general")){
            double A,B,C,D,E,F,G,H,I,J,x,y,z,a,b,c,d,c1,c2,c3,shine,l,w,h;

            file>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            file>>x>>y>>z>>l>>w>>h;
            file>>c1>>c2>>c3;
            file>>a>>b>>c>>d;
            file>>shine;

            cout<<" "<<A<<" "<<B<<" "<<C<<" "<<D<<" "<<E<<" "<<F<<" "<<G<<" "<<H<<" "<<I<<" "<<J<<endl;
            cout<<" "<<x<<" "<<y<<" "<<z<<" "<<l<<" "<<w<<" "<<h<<endl;
            cout<<" "<<c1<<" "<<c2<<" "<<c3<<" "<<endl;
            cout<<" "<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
            cout<<" "<<shine<<endl;

            temp = new genQuad(x,y,z,l,w,h);
            temp->setVar(A,B,C,D,E,F,G,H,I,J);
            temp->setColor(c1,c2,c3);
            temp->setCoefficients(a,b,c,d);
            temp->setShine(shine);
            objects.push_back(temp);
        }

    }
    file>>numberOfLights;
    cout<<numberOfLights<<endl;
    for(int i=0; i<numberOfLights; i++){
        cout<<i<<". light"<<endl;
        double a,b,c;
        file>>a>>b>>c;
        cout<<a<<" "<<b<<" "<<c<<endl;
        Point light(a,b,c);
        lights.push_back(light);
    }
    file.close();
}

void freeMemory(){
    vector<Point>().swap(lights);
    vector<Object*>().swap(objects);
}

int main(int argc, char **argv){
    LoadActualData();
//    loadTestData();
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
    freeMemory();
	return 0;
}

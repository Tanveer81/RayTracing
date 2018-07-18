#include <bits/stdc++.h>
using namespace std;
class Point3{
public:
    double x,y,z;
    Point3(){}
    Point3(double x, double y, double z){
        this->x = x;
        this->y = y;
        this->z = z;
    }
};

Point3 normalize(Point3 b){
    Point3 a;
    double d = sqrt(b.x*b.x+ b.y*b.y+ b.z*b.z);
    a.x = b.x / d;
    a.y = b.y / d;
    a.z = b.z / d;
    return a;
}

class Ray{
public:
    Point3 start, dir;
    Ray(Point3 a, Point3 b){
        start = a;
        dir = normalize(b);
    }
};



double dot(Point3 a, Point3 b){
    double dot = a.x*b.x + a.y*b.y + a.z*b.z;
    return dot;
}

void cross(double a[], double b[], double c[]){
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=-a[0]*b[2]+a[2]*b[0];///
    c[2]=a[0]*b[1]-a[1]*b[0];
}

class Object{
public:
    Point3 reference_point;
    double height, width, length;
    int shine;
    double color[3];
    double co_efficients[4];
    double source_factor = 1.0;
    double refractive_index = 1.5;
    Object(){}

    virtual void draw() = 0;

    virtual double calculateT(Ray *r) = 0;

    virtual Point3 getNormal(Point3 intersectionPoint) = 0;

    virtual double intersect(Ray *r, double current_color[3], int level) = 0;

    void setColor(double a, double b, double c){
        this->color[0] = a*255;
        this->color[1] = b*255;
        this->color[2] = c*255;
    }

    void setShine(int sh){
        shine = sh;
    }

    void setCoefficients(double ambient, double diffuse, double specular, double reflection){
        co_efficients[0] = ambient;
        co_efficients[1] = diffuse;
        co_efficients[2] = specular;
        co_efficients[3] = reflection;
    }

    Point3 getReflection(Ray* r, Point3 normal){
        double d = dot(r->dir,normal);
        Point3 reflection(r->dir.x-normal.x*2.0*d, r->dir.y-normal.y*2.0*d, r->dir.z-normal.z*2.0*d);
        reflection = normalize(reflection);
        return reflection;
    }
};


vector<Object*> objects;
vector<Point3> lights;


class Sphere : public Object{
public:
    Sphere(Point3 center, double radius){
        reference_point = center;
        length = radius;
    }
    void draw(){
        glPushMatrix();{
        glColor3f(color[0],color[1],color[2]);
        GLUquadric *quad;
        quad = gluNewQuadric();
        glTranslatef(reference_point.x,reference_point.y,reference_point.z);
        gluSphere(quad,length,500,20);//radius,slice,stack
        }glPopMatrix();
    }


    double calculateT(Ray *r){
        double t;

        Point3 Ro(r->start.x - reference_point.x, r->start.y - reference_point.y, r->start.z - reference_point.z);
        double a = dot(r->dir, r->dir);
        double b = 2*dot(r->dir,Ro);
        double c = dot(Ro,Ro) - length*length;
        double d = b*b - 4*a*c;
        if(d<0)return -1;
        d = sqrt(d);
        double t1 = (-b - d) / (2.0 * a);
        double t2 = (-b + d) / (2.0 * a);
        if(t2 > t1){
            double temp = t1;
            t1 = t2;
            t2 = temp;
        }
        return t2;
    }


    Point3 getNormal(Point3 p){
        Point3 normal(p.x-reference_point.x, p.y-reference_point.y, p.z-reference_point.z);
        normal = normalize(normal);
        return normal;
    }


    double intersect(Ray *r, double current_color[3], int level){
        double t = calculateT(r);

        if(t<=0)return -1;

        if(level == 0)return t;

        for(int i=0; i<3;i++)
            current_color[i] = color[i] * co_efficients[0];

        Point3 intersec(r->start.x+r->dir.x*t, r->start.y+r->dir.y*t, r->start.z+r->dir.z*t);
        Point3 normal = getNormal(intersec);
        Point3 reflection = getReflection(r , normal);

        for(int i=0; i<lights.size();i++){
            Point3 direction(lights[i].x-intersec.x, lights[i].y-intersec.y, lights[i].z-intersec.z);
            direction = normalize(direction);
            Point3 start(intersec.x+direction.x*1.0, intersec.y+direction.y*1.0, intersec.z+direction.z*1.0);
            Ray L(start, direction);
            bool obstacleFlag = false;
            double len = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);

            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);

                if(t<0 || t>len )continue;
                obstacleFlag = true;
                break;
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = dot(reflection,r->dir);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int k=0; k<3; k++){
                    current_color[k] += source_factor*lambert*co_efficients[1]*color[k];
                    current_color[k] += source_factor*phong*co_efficients[2]*color[k];
                }
            }
        }

        return t;
    }

};

void drawTile(double a)
{
	glBegin(GL_QUADS);{
		glVertex3f(0,0,0);
		glVertex3f(a,0,0);
		glVertex3f(a,a,0);
		glVertex3f(0,a,0);
	}glEnd();
}

class Floor : public Object{
public:
    Floor(double FloorWidth, double TileWidth){
        reference_point.x = -FloorWidth/2;
        reference_point.y = -FloorWidth/2;
        reference_point.z = 0.0;
        length=TileWidth;
        width = FloorWidth;
    }
    void draw(){
        int tiles = abs(width / length);
        for(int i=0;i<tiles;i++){
            for(int j=0; j<tiles;j++){
                if((i+j) % 2 == 0)glColor3f(0,0,0);
                else glColor3f(1,1,1);
                glPushMatrix();{
                    glTranslatef(reference_point.x+length*i,reference_point.y+length*j,reference_point.z);
                    drawTile(length);
                }glPopMatrix();
            }
        }
    }

    double calculateT(Ray *r){
        return 0;
    }

    Point3 getNormal(Point3 intersectionPoint){
        Point3 normal(2,5,4);
        return normal;
    }
    double intersect(Ray *r, double current_color[3], int level){
        return 0;
    }

};



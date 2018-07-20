#include <bits/stdc++.h>
using namespace std;
extern int recursion_level;
#include "bitmap_image.hpp"
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

Point3 cross(Point3 a, Point3 b){
    Point3 c(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
    return c;
}

class Object{
public:
    Point3 reference_point;
    double height, width, length;
    int shine;
    double color[3];
    double co_efficients[4];
    double source_factor = 1.0;
    double refractive_index = 2.0;
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

    Point3 getRefraction(Ray *r, Point3 normal){
        double N_dot_I = dot(normal,r->dir);
        double k = 1.0 - refractive_index * refractive_index * (1.0 - N_dot_I * N_dot_I);
        Point3 refraction(0, 0, 0);
        if(k>=0){
            refraction.x = refractive_index * r->dir.x - normal.x*(refractive_index * N_dot_I + sqrt(k));
            refraction.y = refractive_index * r->dir.y - normal.y*(refractive_index * N_dot_I + sqrt(k));
            refraction.z = refractive_index * r->dir.z - normal.z*(refractive_index * N_dot_I + sqrt(k));
            refraction = normalize(refraction);
        }
        return refraction;
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
        if(t1 < t2)t = t1;
        else t = t2;
        return t;
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

        for(int i=0; i<3; i++)
            current_color[i] = color[i] * co_efficients[0];

        Point3 intersec(r->start.x+r->dir.x*t, r->start.y+r->dir.y*t, r->start.z+r->dir.z*t);
        Point3 normal = getNormal(intersec);
        Point3 reflection = getReflection(r , normal);
        Point3 refraction = getRefraction(r , normal);

        for(int i=0; i<lights.size();i++){
            Point3 direction(lights[i].x-intersec.x, lights[i].y-intersec.y, lights[i].z-intersec.z);
            direction = normalize(direction);
            Point3 start(intersec.x+direction.x*1.0, intersec.y+direction.y*1.0, intersec.z+direction.z*1.0);
            Ray L(start, direction);
            bool obstacleFlag = false;
            double len = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);

            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);

                if( abs(t)>len || t<0 )continue;
                obstacleFlag = true;
                break;
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = pow(dot(reflection,r->dir),shine);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int k=0; k<3; k++){
                    current_color[k] += source_factor*lambert*co_efficients[1]*color[k];
                    current_color[k] += source_factor*phong*co_efficients[2]*color[k];
                }
            }

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double t_min = 9999999;
                double reflected_color[3];

                for(int k=0; k<objects.size();k++){
                    double t = objects[k]->calculateT(&reflectionRay);
                    if(t<=0)continue;
                    else if(t<t_min){
                        t_min = t;
                        nearest = k;
                    }
                }
                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int k=0; k<3; k++){
                        current_color[k] += reflected_color[k] * co_efficients[3];
                    }
                }

                start.x = intersec.x + refraction.x * 1.0;
                start.y = intersec.y + refraction.y * 1.0;
                start.z = intersec.z + refraction.z * 1.0;

                Ray refractionRay(start, refraction);
                nearest = -1;
                t_min = 9999999;
                double refracted_color[3];

                for(int k=0; k<objects.size();k++){

                    double t = objects[k]->calculateT(&refractionRay);

                    if(t<=0)continue;

                    else if(t<t_min){
                        t_min = t;
                        nearest = k;
                    }
                }

                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&refractionRay,refracted_color,level+1);

                    for (int k=0; k<3; k++){
                        current_color[k] += refracted_color[k] * refractive_index;
                    }
                }
            }

            for(int c=0; c<3; c++){
                if(current_color[c] < 0)current_color[c] = 0;
                else if(current_color[c] > 255)current_color[c] = 255;
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
        reference_point.z = 0;
        length=TileWidth;
        width = FloorWidth;
    }
    void draw(){
        //int tiles = abs(width / length);
        int tiles = abs(reference_point.x*2/length);
        for(int i=0;i<tiles;i++){
            for(int j=0; j<tiles;j++){
                if((i+j) % 2)glColor3f(0,0,0);
                else glColor3f(1,1,1);
                glPushMatrix();{
                    glTranslatef(reference_point.x+length*i,reference_point.y+length*j,reference_point.z);
                    drawTile(length);
                }glPopMatrix();
            }
        }
    }

    double calculateT(Ray *r){
        Point3 normal = getNormal(reference_point);
        double nRo = dot(normal,r->start);
        double nRd = dot(normal, r->dir);
        double t = ( (-1) * nRo ) / nRd;
        return t;
    }

    Point3 getNormal(Point3 intersectionPoint){
        Point3 normal(0,0,1);
        normalize(normal);
        return normal;
    }

    double intersect(Ray *r, double current_color[3], int level){

        double t = calculateT(r);
        Point3 intersec(r->start.x+r->dir.x*t, r->start.y+r->dir.y*t, r->start.z+r->dir.z*t);
        if(intersec.x<reference_point.x || intersec.y<reference_point.y || intersec.x>((-1) * reference_point.x) || intersec.y>((-1) * reference_point.y)){
            return -1;
        }

        int x = (intersec.x-reference_point.x)/length;
        int y = (intersec.y-reference_point.y)/length;

        if((x + y)%2){
            for(int i=0; i<3; i++)color[i] = 0;
        }
        else {
            for(int i=0; i<3; i++)color[i] = 255;
        }
        for(int i=0; i<3; i++){
            current_color[i] = color[i] * co_efficients[0];
        }

        Point3 normal = getNormal(intersec);
        Point3 reflection = getReflection(r , normal);
        Point3 refraction = getRefraction(r , normal);

        for(int i=0; i<lights.size();i++){
            Point3 direction(lights[i].x-intersec.x, lights[i].y-intersec.y, lights[i].z-intersec.z);
            direction = normalize(direction);
            Point3 start(intersec.x+direction.x*1.0, intersec.y+direction.y*1.0, intersec.z+direction.z*1.0);
            Ray L(start, direction);
            bool obstacleFlag = false;
            double len = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);

            for(int j=0; j<objects.size(); j++){

                double t = objects[j]->calculateT(&L);

                if(abs(t)>len || t<0)continue;

                obstacleFlag = true;

                break;
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = pow(dot(reflection,r->dir),shine);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int k=0; k<3; k++){
                    current_color[k] += source_factor*lambert*co_efficients[1]*color[k];
                    current_color[k] += source_factor*phong*co_efficients[2]*color[k];
                }
            }

            if(level == 0) return t ;

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double t_min = 9999999;
                double reflected_color[3];

                for(int k=0; k<objects.size();k++){

                    double t = objects[k]->calculateT(&reflectionRay);

                    if(t<=0)continue;

                     if(t<t_min){
                        t_min = t;

                        nearest = k;

                    }
                }
                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int k=0; k<3; k++){
                        current_color[k] += reflected_color[k] * co_efficients[3];
                    }
                }
            }

            for(int c=0; c<3; c++){
                if(current_color[c] < 0)current_color[c] = 0;
                else if(current_color[c] > 255)current_color[c] = 255;
            }
        }
        return t;
    }

};


class Triangle : public Object{
public:
    Point3 a,b,c;
    Point3 p[3];
    Triangle(Point3 a[]){
        for(int i=0; i<3; i++){
            p[i]=a[i];
            cout<<p[i].x<<" "<<p[i].y<<" "<<p[i].z<<endl;
        }
    }

    void draw(){
        glColor3f(color[0],color[1],color[2]);
        glBegin(GL_TRIANGLES);
        {
            for(int i=0; i<3; i++){
                glVertex3f(p[i].x,p[i].y,p[i].z);
            }
        }
        glEnd();
    }

    Point3 getNormal(Point3 intersectionPoint){

        Point3 a(p[1].x-p[0].x,p[1].y-p[0].y,p[1].z-p[0].z);

        Point3 b(p[2].x-p[0].x,p[2].y-p[0].y,p[2].z-p[0].z);

        Point3 normal = cross(a,b);

        normalize(normal);

        return normal;
    }


    double calculateT(Ray *r){
        double t;
        double epsilon = 0.0000001;
        double a,f,u,v;

        Point3 edge1(p[1].x-p[0].x, p[1].y-p[0].y, p[1].z-p[0].z);
        Point3 edge2(p[2].x-p[0].x, p[2].y-p[0].y, p[2].z-p[0].z);

        Point3 h = cross(r->dir,edge2);
        a = dot(edge1,h);
        if (a > -epsilon && a < epsilon)return -1;
        f = 1.0 / a;

        Point3 s(r->start.x-p[0].x,r->start.y-p[0].y,r->start.z-p[0].z);
        u = f * dot(s,h);
        if (u < 0 || u > 1)return -1;

        Point3 q = cross(s,edge1);
        v = f * dot(r->dir, q);
        if (v < 0 || u + v > 1)return -1;

        t = dot(edge2, q) *f;
        if(t > epsilon)return t;

        return -1;
    }


    double intersect(Ray *r, double current_color[3], int level){
        double t = calculateT(r);

        if(t<=0)return -1;

        if(level == 0)return t;

        for(int i=0; i<3; i++)
            current_color[i] = color[i] * co_efficients[0];

        Point3 intersec(r->start.x+r->dir.x*t, r->start.y+r->dir.y*t, r->start.z+r->dir.z*t);
        Point3 normal = getNormal(intersec);
        Point3 reflection = getReflection(r , normal);
        Point3 refraction = getRefraction(r , normal);

        for(int i=0; i<lights.size();i++){
            Point3 direction(lights[i].x-intersec.x, lights[i].y-intersec.y, lights[i].z-intersec.z);
            direction = normalize(direction);
            Point3 start(intersec.x+direction.x*1.0, intersec.y+direction.y*1.0, intersec.z+direction.z*1.0);
            Ray L(start, direction);
            bool obstacleFlag = false;
            double len = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);

            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);

                if(t<0 || abs(t)>len )continue;
                obstacleFlag = true;
                break;
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = pow(dot(reflection,r->dir),shine);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int k=0; k<3; k++){
                    current_color[k] += source_factor*lambert*co_efficients[1]*color[k];
                    current_color[k] += source_factor*phong*co_efficients[2]*color[k];
                }
            }

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double t_min = 9999999;
                double reflected_color[3];

                for(int k=0; k<objects.size();k++){

                    double t = objects[k]->calculateT(&reflectionRay);

                    if(t<=0)continue;

                    if(t<t_min){
                        t_min = t;

                        nearest = k;
                    }
                }

                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int k=0; k<3; k++){
                        current_color[k] += reflected_color[k] * co_efficients[3];
                    }
                }
            }

            for(int c=0; c<3; c++){
                if(current_color[c] < 0)current_color[c] = 0;
                else if(current_color[c] > 255)current_color[c] = 255;
            }
        }

        return t;
    }
};



class genQuad : public Object{
public:
    double A,B,C,D,E,F,G,H,I,J;
    genQuad(double A,double B,double C,double D,double E,double F,double G,
            double H,double I,double J,Point3 r, double length,double width,
            double height)
    {
        this->A =A ;
        this->B =B ;
        this->C =C ;
        this->D =D ;
        this->E =E ;
        this->F =F ;
        this->G =G ;
        this->H =H ;
        this->I =I ;
        this->J =J ;
        this->reference_point = r;
        this->length = length;
        this->width = width;
        this->height = height;
    }

    void draw(){}

    Point3 getNormal(Point3 intersec){
        Point3 normal(2*A*intersec.x + D*intersec.y + F*intersec.z + G,
                      2*B*intersec.y + D*intersec.x + E*intersec.z + H,
                      2*C*intersec.z + E*intersec.y + F*intersec.x + I);
        normalize(normal);
        return normal;
    }


    double calculateT(Ray *r){
        double X0 = r->start.x;
        double Y0 = r->start.y;
        double Z0 = r->start.z;

        double X1 = r->dir.x;
        double Y1 = r->dir.y;
        double Z1 = r->dir.z;

        double a = A*X1*X1 + B*Y1*Y1 + C*Z1*Z1 + D*X1*Y1 + E*Y1*Z1 + F*Z1*X1;
        double b = 2*A*X0*X1 + 2*B*Y0*Y1 + 2*C*Z0*Z1 + D*X0*Y1 + D*X1*Y0 + E*Y0*Z1 + E*Y1*Z0 + F*Z0*X1 + F*Z1*X0 + G*X1 + H*Y1 + I*Z1;
        double c = A*X0*X0 + B*Y0*Y0 + C*Z0*Z0 + D*X0*Y0 + E*Y0*Z0 + F*Z0*X0 + G*X0 + H*Y0 + I*Z0 + J;

        double d = b*b - 4*a*c;

        if(d<0)return -1;

        double t1 = (-b + sqrt(d)) / (2.0*a);
        double t2 = (-b - sqrt(d)) / (2.0*a);

        Point3 p1(r->start.x+r->dir.x*t1, r->start.y+r->dir.y*t1, r->start.z+r->dir.z*t1);
        Point3 p2(r->start.x+r->dir.x*t2, r->start.y+r->dir.y*t2, r->start.z+r->dir.z*t2);

        bool fx1 = length>0 && (reference_point.x > p1.x || reference_point.x+length < p1.x);
        bool fy1 = width>0 && (reference_point.y > p1.y || reference_point.y+width < p1.y);
        bool fz1 = height>0 && (reference_point.y > p1.y || reference_point.y+height < p1.z);

        bool fx2 = length>0 && (reference_point.x > p2.x || reference_point.x+length < p2.x);
        bool fy2 = width>0 && (reference_point.y > p2.y || reference_point.y+width < p2.y);
        bool fz2 = height>0 && (reference_point.y > p2.y || reference_point.y+height < p2.z);

        bool f1 = fx1 || fy1 || fz1;
        bool f2 = fx2 || fy2 || fz2;

        double t;
        if(f1 && f2)t = -1;
        if(f2)t = t1;
        if(f1)t = t2;
        if(t2>t1)t = t1;
        else t = t2;

        return t;
    }


    double intersect(Ray *r, double current_color[3], int level){
        double t = calculateT(r);
        if(t<=0)return -1;
        if(level == 0)return t;
        for(int i=0; i<3; i++)
            current_color[i] = color[i] * co_efficients[0];

        Point3 intersec(r->start.x+r->dir.x*t, r->start.y+r->dir.y*t, r->start.z+r->dir.z*t);
        Point3 normal = getNormal(intersec);
        Point3 reflection = getReflection(r , normal);
        Point3 refraction = getRefraction(r , normal);

        for(int i=0; i<lights.size();i++){
            Point3 direction(lights[i].x-intersec.x, lights[i].y-intersec.y, lights[i].z-intersec.z);
            direction = normalize(direction);
            Point3 start(intersec.x+direction.x*1.0, intersec.y+direction.y*1.0, intersec.z+direction.z*1.0);
            Ray L(start, direction);
            bool obstacleFlag = false;
            double len = sqrt(direction.x*direction.x + direction.y*direction.y + direction.z*direction.z);

            for(int j=0; j<objects.size(); j++){

                double t = objects[j]->calculateT(&L);

                if(abs(t)>len || t<0  )continue;

                obstacleFlag = true;

                break;
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = pow(dot(reflection,r->dir),shine);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int k=0; k<3; k++){
                    current_color[k] += source_factor*lambert*co_efficients[1]*color[k];
                    current_color[k] += source_factor*phong*co_efficients[2]*color[k];
                }
            }

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double t_min = 9999999;
                double reflected_color[3];

                for(int k=0; k<objects.size();k++){
                    double t = objects[k]->calculateT(&reflectionRay);
                    if(t<=0)continue;
                    else if(t<t_min){
                        t_min = t;
                        nearest = k;
                    }
                }
                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int k=0; k<3; k++){
                        current_color[k] += reflected_color[k] * co_efficients[3];
                    }
                }
            }

            for(int c=0; c<3; c++){
                if(current_color[c] < 0)current_color[c] = 0;
                else if(current_color[c] > 255)current_color[c] = 255;
            }
        }

        return t;
    }
};

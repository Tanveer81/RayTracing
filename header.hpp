#include <bits/stdc++.h>
using namespace std;
extern int recursion_level;
#include "bitmap_image.hpp"
class Point{
public:
    double x,y,z;

    Point(double a, double b, double c){
        x = a;y = b;z = c;
    }

    Point(){}
};

Point normalize(Point b){
    Point a;
    double d = sqrt(b.x*b.x+ b.y*b.y+ b.z*b.z);
    a.x = b.x / d;
    a.y = b.y / d;
    a.z = b.z / d;
    return a;
}

class Ray{
public:
    Point start, dir;
    Ray(Point a, Point b){
        start = a;
        dir = normalize(b);
    }
};


double dot(Point a, Point b){
    double dot = a.x*b.x + a.y*b.y + a.z*b.z;
    return dot;
}

Point cross(Point a, Point b){
    double x = a.y*b.z - a.z*b.y;
    double y = a.z*b.x - a.x*b.z;
    double z = a.x*b.y - a.y*b.x;
    Point c(x,y,z);
    return c;
}

class Object{
public:
    Point reference_point;
    double height, width, length;
    int shine;
    double color[3];
    double co_efficients[4];

    virtual void draw() = 0;
    virtual double calculateT(Ray *r) = 0;
    virtual Point getNormal(Point intersectionPoint) = 0;
    virtual double intersect(Ray *r, double current_color[3], int level) = 0;
    virtual void setVar(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J) = 0;

    void setColor(double a, double b, double c){
        this->color[0] = a*255;
        this->color[1] = b*255;
        this->color[2] = c*255;
    }

    void setShine(int sh){
        this->shine = sh;
    }

    Object(){}

    void setCoefficients(double ambient, double diffuse, double specular, double reflection){
        this->co_efficients[0] = ambient;
        this->co_efficients[1] = diffuse;
        this->co_efficients[2] = specular;
        this->co_efficients[3] = reflection;
    }

    double source_factor = 1.0;
    double refractive_index = 2.0;

    Point getReflection(Ray* r, Point normal){
        double d = dot(r->dir,normal);
        double a = r->dir.x-normal.x*2.0*d;
        double b = r->dir.y-normal.y*2.0*d;
        double c = r->dir.z-normal.z*2.0*d;
        Point reflection(a,b,c);
        reflection = normalize(reflection);
        return reflection;
    }

    Point getRefraction(Ray *r, Point normal){
        double N_dot_I = dot(normal,r->dir);
        double k = 1.0 - refractive_index * refractive_index * (1.0 - N_dot_I * N_dot_I);
        Point refraction;
        if(k<0){
            refraction.x = 0;
            refraction.y = 0;
            refraction.z = 0;
        }
        else{
            refraction.x = refractive_index * r->dir.x - normal.x*(refractive_index * N_dot_I + sqrt(k));
            refraction.y = refractive_index * r->dir.y - normal.y*(refractive_index * N_dot_I + sqrt(k));
            refraction.z = refractive_index * r->dir.z - normal.z*(refractive_index * N_dot_I + sqrt(k));
            refraction = normalize(refraction);
        }
        return refraction;
    }
};


vector<Object*> objects;
vector<Point> lights;


class Sphere : public Object{
public:
    Sphere(Point center, double radius){
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

    void setVar(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J){}

    double calculateT(Ray *r){
        double t;
        Point Ro(r->start.x - reference_point.x, r->start.y - reference_point.y, r->start.z - reference_point.z);
        double a = dot(r->dir, r->dir);
        double b = 2*dot(r->dir,Ro);
        double c = dot(Ro,Ro) - length*length;
        double d = b*b - 4*a*c;
        if(d<0)return -1;
        d = sqrt(d);
        double t1 = (-b - d) / (2.0 * a);
        double t2 = (-b + d) / (2.0 * a);
        t = (t1 > t2) ? t2 : t1;
        return t;
    }


    Point getNormal(Point p){
        double a = p.x-reference_point.x;
        double b = p.y-reference_point.y;
        double c = p.z-reference_point.z;
        Point normal(a,b,c);
        normal = normalize(normal);
        return normal;
    }


    double intersect(Ray *r, double current_color[3], int level){

        double t = calculateT(r);

        if(t<=0)return -1;

        if(level == 0)return t;

        double intersecx = r->start.x+r->dir.x*t;
        double intersecy = r->start.y+r->dir.y*t;
        double intersecz = r->start.z+r->dir.z*t;
        Point intersec(intersecx, intersecy, intersecz);

        Point normal = getNormal(intersec);

        Point reflection = getReflection(r , normal);

        Point refraction = getRefraction(r , normal);

        for(int i=0; i<3; i++)
            current_color[i] = color[i] * co_efficients[0];

        for(int i=0; i<lights.size();i++){

            double dirx = lights[i].x-intersec.x;
            double diry = lights[i].y-intersec.y;
            double dirz = lights[i].z-intersec.z;
            Point direction(dirx, diry, dirz);

            double len = sqrt(direction.x*direction.x +
                               direction.y*direction.y +
                                direction.z*direction.z);

            direction = normalize(direction);

            double startx = intersec.x+direction.x*1.0;
            double starty = intersec.y+direction.y*1.0;
            double startz = intersec.z+direction.z*1.0;
            Point start(startx, starty, startz);

            Ray L(start, direction);

            bool obstacleFlag = false;

            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);
                if( t <= len && t>0 ){
                    obstacleFlag = true;
                    break;
                }
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = dot(reflection,r->dir);
                phong = pow(phong, shine);
                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int j=0; j<3; j++){
                    current_color[j] += lambert*source_factor*color[j]*co_efficients[1];
                    current_color[j] += phong*source_factor*color[j]*co_efficients[2];
                }
            }

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double reflected_color[3];
                double t_min = 9999999;

                for(int j=0; j<objects.size();j++){
                    double t = objects[j]->intersect(&reflectionRay, reflected_color, 0);
                    if(t<=0)continue;
                    else if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }
                if(nearest!=-1){
                    objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += reflected_color[j] * co_efficients[3];
                    }
                }

                start.x = intersec.x + refraction.x * 1.0;
                start.y = intersec.y + refraction.y * 1.0;
                start.z = intersec.z + refraction.z * 1.0;

                Ray refractionRay(start, refraction);
                nearest = -1;
                double refracted_color[3];
                t_min = 9999999;

                for(int j=0; j<objects.size();j++){

                    double t = objects[j]->intersect(&refractionRay, refracted_color, 0);

                    if(t<=0)continue;

                    else if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }

                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&refractionRay,refracted_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += refracted_color[j] * refractive_index;
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
    bitmap_image texture;
    Floor(double FloorWidth, double TileWidth){
        reference_point.x = -FloorWidth/2;
        reference_point.y = -FloorWidth/2;
        reference_point.z = 0;
        texture = bitmap_image("3.bmp");
        length=TileWidth;
        width = FloorWidth;
        cout<<texture.width()<<"texture width"<<endl;
    }
    void draw(){
        int tiles = reference_point.x*2;
        tiles /= 2;
        tiles = abs(tiles);
        for(int i=0;i<tiles;i++){
            for(int j=0; j<tiles;j++){
                if((i+j) % 2 == 0){
                    glColor3f(1,1,1);
                }
                else{
                    glColor3f(0,0,0);
                }
                glPushMatrix();
                {
                    glTranslatef(reference_point.x+length*i,reference_point.y+length*j,reference_point.z);
                    drawTile(length);
                }
                glPopMatrix();
            }
        }
    }

    void setVar(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J){}

    double calculateT(Ray *r){
        Point normal = getNormal(reference_point);
        double nRo = dot(normal,r->start);
        double nRd = dot(normal, r->dir);
        nRo = (-1) * nRo;
        double t = nRo  / nRd;
        return t;
    }

    Point getNormal(Point intersectionPoint){
        Point normal(0,0,1);
        return normal;
    }

    double intersect(Ray *r, double current_color[3], int level){

        double t = calculateT(r);

        double intersecx = r->start.x+r->dir.x*t;
        double intersecy = r->start.y+r->dir.y*t;
        double intersecz = r->start.z+r->dir.z*t;
        Point intersec(intersecx, intersecy, intersecz);

        if(intersec.x<reference_point.x || intersec.y<reference_point.y || intersec.x>((-1) * reference_point.x) || intersec.y>((-1) * reference_point.y)){
            return -1;
        }

        int x = (intersec.x-reference_point.x)/length;
        int y = (intersec.y-reference_point.y)/length;

        unsigned int xx = intersec.x + (width / 2);
        xx = xx * (texture.width() / width );
        unsigned int yy = intersec.y + (width / 2);
        yy = yy * (texture.height() / width);

        if((x + y)%2){
            for(int i=0; i<3; i++)color[i] = 0;
        }
        else {
            for(int i=0; i<3; i++)color[i] = 255;
        }

        unsigned char rr, g, b;
//        xx = xx * texture.width()/1000.0;
//        yy = yy * texture.height()/1000.0;

        texture.get_pixel(xx, yy, rr, g, b);

        double c1 = rr;
        double c2 = g;
        double c3 = b;

        for(int i=0; i<3; i++){
            double cc;
            if(i == 0)cc = c1 / 255;
            else if(i == 1)cc = c2 / 255;
            else if(i == 2)cc = c3 / 255;
            current_color[i] = color[i] * co_efficients[0] * cc;
        }

        Point normal = getNormal(intersec);
        Point reflection = getReflection(r , normal);
        Point refraction = getRefraction(r , normal);

        for(int i=0; i<lights.size();i++){

            double dirx = lights[i].x-intersec.x;
            double diry = lights[i].y-intersec.y;
            double dirz = lights[i].z-intersec.z;
            Point direction(dirx, diry, dirz);

            double len = sqrt(direction.x*direction.x +
                               direction.y*direction.y +
                                direction.z*direction.z);

            direction = normalize(direction);

            double startx = intersec.x+direction.x*1.0;
            double starty = intersec.y+direction.y*1.0;
            double startz = intersec.z+direction.z*1.0;
            Point start(startx, starty, startz);

            Ray L(start, direction);

            bool obstacleFlag = false;

            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);
                if( t <= len && t>0 ){
                    obstacleFlag = true;
                    break;
                }
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = pow(dot(reflection,r->dir),shine);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int j=0; j<3; j++){
                    current_color[j] += lambert*source_factor*color[j]*co_efficients[1];
                    current_color[j] += phong*source_factor*color[j]*co_efficients[2];
                }
            }

            if(level == 0) return t ;

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double reflected_color[3];
                double t_min = 9999999;

                for(int j=0; j<objects.size();j++){

                    double t = objects[j]->intersect(&reflectionRay, reflected_color, 0);

                    if(t<=0)continue;

                    if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }
                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += reflected_color[j] * co_efficients[3];
                    }
                }

//                start.x = intersec.x + refraction.x * 1.0;
//                start.y = intersec.y + refraction.y * 1.0;
//                start.z = intersec.z + refraction.z * 1.0;
//
//                Ray refractionRay(start, refraction);
//                nearest = -1;
//                double refracted_color[3];
//                t_min = 9999999;
//
//                for(int j=0; j<objects.size();j++){
//
//                    double t = objects[j]->calculateT(&refractionRay);
//
//                    if(t<=0)continue;
//
//                    else if(t<t_min){
//                        t_min = t;
//                        nearest = j;
//                    }
//                }
//
//                if(nearest!=-1){
//                    double t = objects[nearest]->intersect(&refractionRay,refracted_color,level+1);
//
//                    for (int j=0; j<3; j++){
//                        current_color[j] += refracted_color[j] * refractive_index;
//                    }
//                }

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
    Point a,b,c;
    Point p[3];

    Triangle(Point a[]){
        for(int i=0; i<3; i++){
            p[i]=a[i];
            cout<<p[i].x<<" "<<p[i].y<<" "<<p[i].z<<endl;
        }
    }

    void draw(){
        double r = color[0];
        double g = color[1];
        double b = color[2];
        glColor3f(r, g, b);
        glBegin(GL_TRIANGLES);
        {
            for(int i=0; i<3; i++){
                glVertex3f(p[i].x,p[i].y,p[i].z);
            }
        }
        glEnd();
    }

    void setVar(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J){}

    Point getNormal(Point intersectionPoint){

        double qx = p[1].x-p[0].x;
        double qy = p[1].y-p[0].y;
        double qz = p[1].z-p[0].z;
        Point q(qx, qy, qz);

        double px = p[2].x-p[0].x;
        double py = p[2].y-p[0].y;
        double pz = p[2].z-p[0].z;
        Point pp(px, py, pz);

        Point normal = cross(q,pp);

        normalize(normal);

        return normal;
    }


    double calculateT(Ray *r){
        double t;
        double epsilon = 0.0000001;
        double a,f,u,v;

        double qx = p[1].x-p[0].x;
        double qy = p[1].y-p[0].y;
        double qz = p[1].z-p[0].z;
        Point edge1(qx, qy, qz);

        double px = p[2].x-p[0].x;
        double py = p[2].y-p[0].y;
        double pz = p[2].z-p[0].z;
        Point edge2(px, py, pz);

        Point h = cross(r->dir,edge2);
        a = dot(edge1,h);
        if (a > -epsilon && a < epsilon)return -1;
        f = 1.0 / a;

        double sx = r->start.x-p[0].x;
        double sy = r->start.y-p[0].y;
        double sz = r->start.z-p[0].z;
        Point s(sx, sy, sz);

        u = f * dot(s,h);
        if (u < 0 || u > 1)return -1;

        Point q = cross(s,edge1);
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


        double intersecx = r->start.x+r->dir.x*t;
        double intersecy = r->start.y+r->dir.y*t;
        double intersecz = r->start.z+r->dir.z*t;
        Point intersec(intersecx, intersecy, intersecz);

        for(int i=0; i<lights.size();i++){

            Point normal = getNormal(intersec);

            double dirx = lights[i].x-intersec.x;
            double diry = lights[i].y-intersec.y;
            double dirz = lights[i].z-intersec.z;
            Point direction(dirx, diry, dirz);

            double len = sqrt(direction.x*direction.x +
                               direction.y*direction.y +
                                direction.z*direction.z);

            direction = normalize(direction);

            double d = dot(direction, normal);

            if(d>0){
                normal.x = normal.x * (-1);
                normal.y = normal.y * (-1);
                normal.z = normal.z * (-1);
            }

            Point reflection = getReflection(r , normal);

            double startx = intersec.x+direction.x*1.0;
            double starty = intersec.y+direction.y*1.0;
            double startz = intersec.z+direction.z*1.0;
            Point start(startx, starty, startz);

            Ray L(start, direction);

            bool obstacleFlag = false;


            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);
                if( t <= len && t>0 ){
                    obstacleFlag = true;
                    break;
                }
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = pow(dot(reflection,r->dir),shine);

                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int j=0; j<3; j++){
                    current_color[j] += lambert*source_factor*color[j]*co_efficients[1];
                    current_color[j] += phong*source_factor*color[j]*co_efficients[2];
                }
            }

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;

                Ray reflectionRay(start, reflection);
                int nearest = -1;
                double reflected_color[3];
                double t_min = 9999999;

                for(int j=0; j<objects.size();j++){

                    double t = objects[j]->intersect(&reflectionRay, reflected_color, 0);

                    if(t<=0)continue;

                    if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }

                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += reflected_color[j] * co_efficients[3];
                    }
                }
                Point refraction = getRefraction(r , normal);
                start.x = intersec.x + refraction.x * 1.0;
                start.y = intersec.y + refraction.y * 1.0;
                start.z = intersec.z + refraction.z * 1.0;

                Ray refractionRay(start, refraction);
                nearest = -1;
                double refracted_color[3];
                t_min = 9999999;

                for(int j=0; j<objects.size();j++){

                    double t = objects[nearest]->intersect(&refractionRay,refracted_color,0);

                    if(t<=0)continue;

                    else if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }

                if(nearest!=-1){
                    objects[nearest]->intersect(&refractionRay,refracted_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += refracted_color[j] * refractive_index;
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
    genQuad(double a, double b, double c, double length,double width,double height){
        Point r(a,b,c);
        this->reference_point = r;
        this->length = length;
        this->width = width;
        this->height = height;
    }

    void setVar(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J){
        this->A =A ;this->B =B ;this->C =C ;this->D =D ;this->E =E ;
        this->F =F ;this->G =G ;this->H =H ;this->I =I ;this->J =J ;
    }

    void draw(){}

    Point getNormal(Point intersec){
        Point normal(2*A*intersec.x + D*intersec.y + F*intersec.z + G,
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
        double b = 2*A*X0*X1 + 2*B*Y0*Y1 + 2*C*Z0*Z1 + D*X0*Y1 + D*X1*Y0 + E*Y0*Z1 + E*Y1*Z0 + F*Z0*X1 + F*Z1*X0 + G*X0 + H*Y0 + I*Z0;
        double c = A*X0*X0 + B*Y0*Y0 + C*Z0*Z0 + D*X0*Y0 + E*Y0*Z0 + F*Z0*X0 + G*X0 + H*Y0 + I*Z0 + J;

        double d = b*b - 4*a*c;

        if(d<0)return -1;

        double t1 = (-b + sqrt(d)) / (2.0*a);
        double t2 = (-b - sqrt(d)) / (2.0*a);

        double p1x = X0+X1*t1;
        double p1y = Y0+Y1*t1;
        double p1z = Z0+Z1*t1;

        double p2x = X0+X1*t2;
        double p2y = Y0+Y1*t2;
        double p2z = Z0+Z1*t2;

        Point p1(p1x, p1y, p1z);
        Point p2(p2x, p2y, p2z);

        Point rp = reference_point;

        bool fx1 = length>0 && (rp.x > p1.x || rp.x+length < p1.x);
        bool fy1 = width>0 && (rp.y > p1.y || rp.y+width < p1.y);
        bool fz1 = height>0 && (rp.z > p1.z || rp.z+height < p1.z);

        bool fx2 = length>0 && (rp.x > p2.x || rp.x+length < p2.x);
        bool fy2 = width>0 && (rp.y > p2.y || rp.y+width < p2.y);
        bool fz2 = height>0 && (rp.z > p2.z || rp.z+height < p2.z);

        bool f1 = fx1 || fy1 || fz1;
        bool f2 = fx2 || fy2 || fz2;

        if(f1 && f2)return -1;
        if(f2)return t1;
        if(f1)return t2;
        if(t2>t1)return t1;
        else return t2;

    }


    double intersect(Ray *r, double current_color[3], int level){

        double t = calculateT(r);

        if(t<=0)return -1;

        if(level == 0)return t;

        double intersecx = r->start.x+r->dir.x*t;
        double intersecy = r->start.y+r->dir.y*t;
        double intersecz = r->start.z+r->dir.z*t;
        Point intersec(intersecx, intersecy, intersecz);

        Point normal = getNormal(intersec);

        Point reflection = getReflection(r , normal);


        for(int i=0; i<3; i++)
            current_color[i] = color[i] * co_efficients[0];

        for(int i=0; i<lights.size();i++){

            double dirx = lights[i].x-intersec.x;
            double diry = lights[i].y-intersec.y;
            double dirz = lights[i].z-intersec.z;
            Point direction(dirx, diry, dirz);

            double len = sqrt(direction.x*direction.x +
                               direction.y*direction.y +
                                direction.z*direction.z);

            direction = normalize(direction);

            double startx = intersec.x+direction.x*1.0;
            double starty = intersec.y+direction.y*1.0;
            double startz = intersec.z+direction.z*1.0;
            Point start(startx, starty, startz);

            Ray L(start, direction);

            bool obstacleFlag = false;

            for(int j=0; j<objects.size(); j++){
                double t = objects[j]->calculateT(&L);
                if( t <= len && t>0 ){
                    obstacleFlag = true;
                    break;
                }
            }

            if(!obstacleFlag){
                double lambert = dot(L.dir,normal);
                double phong = dot(reflection,r->dir);
                phong = pow(phong, shine);
                if(lambert < 0)lambert = 0;
                if(phong < 0)phong = 0;

                for (int j=0; j<3; j++){
                    current_color[j] += lambert*source_factor*color[j]*co_efficients[1];
                    current_color[j] += phong*source_factor*color[j]*co_efficients[2];
                }
            }

            if(level < recursion_level){

                start.x = intersec.x + reflection.x * 1.0;
                start.y = intersec.y + reflection.y * 1.0;
                start.z = intersec.z + reflection.z * 1.0;
                Ray reflectionRay(start, reflection);

                int nearest = -1;
                double reflected_color[3];
                double t_min = 9999999;

                for(int j=0; j<objects.size();j++){
                    double t = objects[j]->intersect(&reflectionRay, reflected_color, 0);
                    if(t<=0)continue;
                    else if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }
                if(nearest!=-1){
                    objects[nearest]->intersect(&reflectionRay,reflected_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += reflected_color[j] * co_efficients[3];
                    }
                }

                Point refraction = getRefraction(r , normal);
                start.x = intersec.x + refraction.x * 1.0;
                start.y = intersec.y + refraction.y * 1.0;
                start.z = intersec.z + refraction.z * 1.0;

                Ray refractionRay(start, refraction);
                nearest = -1;
                double refracted_color[3];
                t_min = 9999999;

                for(int j=0; j<objects.size();j++){

                    double t = objects[j]->intersect(&refractionRay, refracted_color, 0);

                    if(t<=0)continue;

                    else if(t<t_min){
                        t_min = t;
                        nearest = j;
                    }
                }

                if(nearest!=-1){
                    double t = objects[nearest]->intersect(&refractionRay,refracted_color,level+1);

                    for (int j=0; j<3; j++){
                        current_color[j] += refracted_color[j] * refractive_index;
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

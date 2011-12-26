#include "../../config.h"

#include <stdio.h>
#include <FL/Fl.H>
#include <FL/Fl_Gl_Window.H>
#include <FL/gl.h>
#include <vector>
#include <cmath>
#include <complex>

#ifdef HAS_GD
#include <gd.h>
#endif /* HAS_GD */

using std::vector;


class MyGLWindow : public Fl_Gl_Window {
public:
    MyGLWindow(int X, int Y, int W, int H, const char* L=0);

    void set_xbounds();
    void set_xbounds(double x1, double y1, double x2, double y2);
    void unset_xbounds();

    void set_skeleton(bool is_skeleton_);
    void set_colormap(int which_colormap_);
    void set_cbounds();
    void set_cbounds(double c1, double c2);
    void unset_cbounds();

    void colormap(double mag);
    void colorvertex(int p, bool colored = true);
    void draw_quad(int p1, int p2, int p3, int p4, bool colored = true);
    void set_orth2d();
    void drawmesh();

    void draw();
    void resize(int X, int Y, int W, int H);

    void read_mesh();
    void read_bmesh();

    void read_cdata();
    void read_bcdata();

    void read_mode_data();
    void read_mode_bdata();

    void read_mode_cdata();
    void read_mode_bcdata();

    void first_frame();
    void next_frame();

    void set_animated(bool animated) { is_animated_ = animated; }
    bool is_animated() { return is_animated_; }

    void save_png(const char* name);

private:
    typedef std::complex<double> dcomplex;

    int which_colormap;
    bool is_skeleton;
    bool frozen_xbounds;
    bool frozen_cbounds;
    int numnp;
    int nshape;
    double minx[2];
    double maxx[2];
    double cmin, cmax;

    std::vector<double> color_;
    std::vector<double> x_;
    std::vector<int> ix_;

    std::vector<dcomplex> mode_color_;
    std::vector<dcomplex> mode_u_;
    std::vector<double>   mode_x_;
    int mode_phase_; // 2 seconds == 20 frames
    bool is_animated_;

    double min(double x, double y) { return (x < y) ? x : y; }
    double max(double x, double y) { return (x > y) ? x : y; }
    double abs(double x) { return (x < 0) ? x : -x; }

    void modified_spectral_colormap(double x, double& r, double& g, double& b);
    void blue_grey_colormap(double x, double& r, double& g, double& b);
};


int nwin;                  // Number of display windows
int currentwin;            // Index of current window
vector<MyGLWindow*> meshw; // Array of mesh windows

#define ME MyGLWindow


ME::ME(int X, int Y, int W, int H, const char* L) : 
    Fl_Gl_Window(X,Y,W,H,L), 
    which_colormap(0), is_skeleton(0), 
    frozen_xbounds(0), frozen_cbounds(0),
    mode_phase_(0), is_animated_(false)
{ 
    cmin = cmax = 0;
    minx[0] = minx[1] = 0;
    maxx[1] = maxx[1] = 1;
}


void ME::set_skeleton(bool is_skeleton_)
{
    is_skeleton = is_skeleton_;
}


void ME::set_colormap(int which_colormap_)
{
    which_colormap = which_colormap_;
}


void ME::set_xbounds()
{
    if (x_.size()) {
        minx[0] = maxx[0] = x_[0];
        minx[1] = maxx[1] = x_[1];
    }
    for (int ipt = 0; ipt < x_.size(); ipt += 2) {
        if (x_[ipt+0] < minx[0]) minx[0] = x_[ipt+0];
        if (x_[ipt+1] < minx[1]) minx[1] = x_[ipt+1];
        if (x_[ipt+0] > maxx[0]) maxx[0] = x_[ipt+0];
        if (x_[ipt+1] > maxx[1]) maxx[1] = x_[ipt+1];
    }
}


void ME::set_xbounds(double x1, double y1, double x2, double y2)
{
    minx[0] = x1;  maxx[0] = x2;
    minx[1] = y1;  maxx[1] = y2;
    frozen_xbounds = 1;
}


void ME::unset_xbounds()
{
    frozen_xbounds = 0;
}


void ME::set_cbounds()
{
    if (color_.size())
        cmin = cmax = color_[0]; 
    for (int ipt = 0; ipt < color_.size(); ++ipt) {
        if (color_[ipt] < cmin) cmin = color_[ipt];
        if (color_[ipt] > cmax) cmax = color_[ipt];
    }
}


void ME::set_cbounds(double c1, double c2)
{
    cmin = c1;
    cmax = c2;
    frozen_cbounds = 1;
}


void ME::unset_cbounds()
{
    frozen_cbounds = 0;
}


void ME::modified_spectral_colormap(double x, double& r, double& g, double& b)
{
    // http://geography.uoregon.edu/datagraphics/color_scales.htm
    // Modified 11-step spectral scheme
    //
    const double data[] = {
        0.650,   0.000,   0.130,
        0.850,   0.150,   0.196,
        0.970,   0.430,   0.370,
        1.000,   0.680,   0.450,
        1.000,   0.880,   0.600,
        1.000,   1.000,   0.750,
        0.880,   1.000,   1.000,
        0.670,   0.970,   1.000,
        0.450,   0.850,   1.000,
        0.250,   0.630,   1.000,
        0.150,   0.300,   1.000,
        0.150,   0.300,   1.000,}; // Repeat 12th step for ease of indexing

    x = 10*min(max(x,0.),1.);
    int n = (int) x;
    double y = x-n;
    r = (1-y)*data[3*n+0] + y*data[3*n+3];
    g = (1-y)*data[3*n+1] + y*data[3*n+4];
    b = (1-y)*data[3*n+2] + y*data[3*n+5];
}


void ME::blue_grey_colormap(double x, double& r, double& g, double& b)
{
    x = max(x,0);
    r = x;
    g = x;
    b = 0.5+0.5*x;
}


void ME::colormap(double mag)
{
    if (cmin == cmax) {
        glColor3f(1.0, 1.0, 1.0);
        return;
    }
    double x = (mag-cmin)/(cmax-cmin);
    double r, g, b;

    if (which_colormap == 1)
        blue_grey_colormap(x, r, g, b);
    else
        modified_spectral_colormap(x, r, g, b);

    r = min(max(r,0.),1.);
    g = min(max(g,0.),1.);
    b = min(max(b,0.),1.);

    glColor3f(r,g,b);
}


void ME::colorvertex(int p, bool colored)
{
    if (colored) 
        colormap(color_[p]);
    glVertex2f(x_[2*p], x_[2*p+1]);
}


void ME::draw_quad(int p1, int p2, int p3, int p4, bool colored)
{
    colorvertex(p1, colored);
    colorvertex(p2, colored);
    colorvertex(p3, colored);
    colorvertex(p4, colored);
}


void ME::set_orth2d()
{
    // Compute scale factor
    double xs = 1.8*w()/(maxx[0]-minx[0]);
    double ys = 1.8*h()/(maxx[1]-minx[1]);
    double ss = ((xs < ys) ? xs : ys);
    
    // Range midpoints
    double xmr = (maxx[0]+minx[0])/2;
    double ymr = (maxx[1]+minx[1])/2;
    
    glLoadIdentity();
    glOrtho(-w()/ss+xmr,w()/ss+xmr, -h()/ss+ymr,h()/ss+ymr, -1,1);
}


void ME::drawmesh()
{
    glColor3f(0.0, 0.0, 1.0);
    for (int ielt = 0; ielt < ix_.size(); ielt += nshape) {
        if (is_skeleton)
            glBegin(GL_LINE_LOOP);
        else
            glBegin(GL_QUADS);
        if (nshape == 4) {
            draw_quad(ix_[ielt+0], ix_[ielt+1], ix_[ielt+2], ix_[ielt+3]);
        } else {
            draw_quad(ix_[ielt+0], ix_[ielt+1], ix_[ielt+8], ix_[ielt+7]);
            draw_quad(ix_[ielt+1], ix_[ielt+2], ix_[ielt+3], ix_[ielt+8]);
            draw_quad(ix_[ielt+8], ix_[ielt+3], ix_[ielt+4], ix_[ielt+5]);
            draw_quad(ix_[ielt+7], ix_[ielt+8], ix_[ielt+5], ix_[ielt+6]);
        }
        glEnd();
    }
}


void ME::draw() 
{
    if (!valid()) {
        valid(1);
        glViewport(0,0,w(),h());
    }
    
    set_orth2d();
    glClear(GL_COLOR_BUFFER_BIT);
    drawmesh();
}


void ME::resize(int X, int Y, int W, int H) 
{
    Fl_Gl_Window::resize(X,Y,W,H);
    glLoadIdentity();
    glViewport(0,0,W,H);
    set_orth2d();
}


void ME::read_mesh()
{
    set_animated(false);

    int numelt;
    char buf[256];
    fgets(buf, 256, stdin);
    sscanf(buf, "%d %d %d", &numnp, &numelt, &nshape);

    color_.resize(numnp);
    x_.resize(2*numnp);
    ix_.resize(nshape*numelt);

    for (int i = 0; i < numnp; ++i) {
        float xx, yy;
        scanf("%f %f", &xx, &yy);
        x_[2*i+0] = xx;
        x_[2*i+1] = yy;
        color_[i] = cmin;
    }
    
    for (int i = 0; i < nshape*numelt; ++i) 
        scanf("%d", &(ix_[i]));

    if (!frozen_xbounds)
        set_xbounds();
}


void ME::read_bmesh()
{
    set_animated(false);

    int numelt;
    char buf[256];
    fgets(buf, 256, stdin);
    sscanf(buf, "%d %d %d", &numnp, &numelt, &nshape);

    color_.resize(numnp);
    x_.resize(2*numnp);
    ix_.resize(nshape*numelt);

    fread(&(x_[0]),  sizeof(double), 2*numnp,       stdin);
    fread(&(ix_[0]), sizeof(int),    nshape*numelt, stdin);
    memset(&(color_[0]), 0, numnp*sizeof(double));

    if (!frozen_xbounds)
        set_xbounds();
}


void ME::read_cdata()
{
    set_animated(false);
    for (int i = 0; i < numnp; ++i) {
        float c;
        scanf("%f", &c);
        color_[i] = c;
    }
    if (!frozen_cbounds)
        set_cbounds();
}


void ME::read_bcdata()
{
    set_animated(false);
    fread(&(color_[0]), sizeof(double), numnp, stdin);
    if (!frozen_cbounds)
        set_cbounds();
}


void ME::read_mode_data()
{
    set_animated(false);
    mode_color_.resize(numnp);
    mode_x_.resize(2*numnp);
    mode_u_.resize(2*numnp);
    for (int i = 0; i < numnp; ++i) {
        float xx, yy;
        scanf("%f %f", &xx, &yy);
        mode_x_[2*i+0] = xx;
        mode_x_[2*i+1] = yy;
        mode_color_[i] = 0;
    }
    for (int i = 0; i < numnp; ++i) {
        float xx, yy, xxi, yyi;
        scanf("%f %f %f %f", &xx, &xxi, &yy, &yyi);
        mode_u_[2*i+0] = dcomplex(xx,xxi);
        mode_u_[2*i+1] = dcomplex(yy,yyi);
    }    
}


void ME::read_mode_bdata()
{
    set_animated(false);
    mode_color_.resize(numnp);
    mode_x_.resize(2*numnp);
    mode_u_.resize(2*numnp);
    fread(&(mode_x_[0]), sizeof(double),   2*numnp, stdin);
    fread(&(mode_u_[0]), sizeof(dcomplex), 2*numnp, stdin);
}


void ME::read_mode_cdata()
{
    set_animated(false);
    for (int i = 0; i < numnp; ++i) {
        float c, ci;
        scanf("%f", &c, &ci);
        mode_color_[i] = dcomplex(c,ci);
    }
}


void ME::read_mode_bcdata()
{
    set_animated(false);
    mode_color_.resize(numnp);
    fread(&(mode_color_[0]), sizeof(dcomplex), numnp, stdin);
}


void ME::first_frame()
{
    mode_phase_ = 0;
}


void ME::next_frame()
{
    ++mode_phase_;
    double theta = ( (double) mode_phase_ )/10;
    dcomplex z = dcomplex(std::cos(theta), std::sin(theta));
    for (int i = 0; i < 2*numnp; ++i) 
        x_[i] = mode_x_[i] + std::real(z*mode_u_[i]);
    for (int i = 0; i < numnp; ++i)
        color_[i] = std::real(z*mode_color_[i]);
}


void ME::save_png(const char* name)
{
#ifdef HAS_GD
    set_animated(false);
    gdImagePtr im;
    int width = w();
    int height = h();

    im = gdImageCreateTrueColor(width, height);

    vector<unsigned char> src(3*width*height);
    make_current();
    glReadPixels(0,0, width,height, GL_RGB,GL_UNSIGNED_BYTE, &(src[0]));

    for(int y=0; y < height; y++) {
        for (int x=0; x < width; x++) {
            int k = (height-1-y)*width+x;
            unsigned char r = src[3*k+0];
            unsigned char g = src[3*k+1];
            unsigned char b = src[3*k+2];
            gdImageSetPixel(im, x, y, 
                            r*0x00010000 + 
                            g*0x00000100 +
                            b);
        }
    }

    FILE* fp = fopen(name, "wb");
    gdImagePng(im,fp);
    fclose(fp);
#endif /* HAS_GD */
}


#undef ME


void show_next_frame_cb(void*) 
{
    for (int i = 0; i < nwin; ++i) {
        MyGLWindow* draw_x = meshw[i];
        if (draw_x && draw_x->is_animated()) {
            draw_x->next_frame();
            draw_x->redraw();
        }
    }
    Fl::repeat_timeout(0.1, show_next_frame_cb);
}


void handle_stdin(int fd, void*) 
{
    MyGLWindow* draw_x = meshw[currentwin];
    char s[256];
    if (fgets(s, 256, stdin)) {
        if (strncmp(s, "quit", 4) == 0) {
            exit(0);
        } else if (strncmp(s, "current", 4) == 0) {
            fgets(s, 256, stdin);
            sscanf(s, "%d", &currentwin);
            if (currentwin < 0)
                currentwin = 0;
            else if (currentwin >= nwin)
                currentwin = nwin-1;
        } else if (strncmp(s, "mesh", 4) == 0) {
            draw_x->read_mesh();
        } else if (strncmp(s, "bmesh", 5) == 0) {
            draw_x->read_bmesh();
        } else if (strncmp(s, "cdat", 4) == 0) {
            draw_x->read_cdata();
        } else if (strncmp(s, "bcdat", 5) == 0) {
            draw_x->read_bcdata();
        } else if (strncmp(s, "mdat", 4) == 0) {
            draw_x->read_mode_data();
        } else if (strncmp(s, "bmdat", 5) == 0) {
            draw_x->read_mode_bdata();
        } else if (strncmp(s, "mcdat", 5) == 0) {
            draw_x->read_mode_cdata();
        } else if (strncmp(s, "mbcdat", 6) == 0) {
            draw_x->read_mode_bcdata();
        } else if (strncmp(s, "skeleton", 4) == 0) {
            int flag;
            fgets(s, 256, stdin);
            sscanf(s, "%d", &flag);
            draw_x->set_skeleton(flag);
        } else if (strncmp(s, "cmap", 4) == 0) {
            int which_cmap;
            fgets(s, 256, stdin);
            sscanf(s, "%d", &which_cmap);
            draw_x->set_colormap(which_cmap);
        } else if (strncmp(s, "box", 3) == 0) {
            float x1, y1, x2, y2;
            fgets(s, 256, stdin);
            sscanf(s, "%f %f %f %f", &x1, &y1, &x2, &y2);
            draw_x->set_xbounds(x1, y1, x2, y2);
        } else if (strncmp(s, "cbou", 4) == 0) {
            float x1, x2;
            fgets(s, 256, stdin);
            sscanf(s, "%f %f", &x1, &x2);
            draw_x->set_cbounds(x1, x2);
        } else if (strncmp(s, "anim", 4) == 0) {
            int flag;
            fgets(s, 256, stdin);
            sscanf(s, "%d", &flag);
            draw_x->set_animated(flag);
            draw_x->first_frame();
        } else if (strncmp(s, "unbox", 5) == 0) {
            draw_x->unset_xbounds();
        } else if (strncmp(s, "uncbou", 6) == 0) {
            draw_x->unset_cbounds();
        } else if (strncmp(s, "plot", 4) == 0) {
            draw_x->redraw();
        } else if (strncmp(s, "save", 4) == 0) {
            char* fname = s + 5;
            for (char* p = fname; *p; ++p)
                if (isspace(*p))
                    *p = 0;
            draw_x->save_png(fname);
        }
    }
}


int main(int argc, char** argv) 
{
    currentwin = 0;
    if (argc > 1)
        nwin = atoi(argv[1]);
    else
        nwin = 1;

    setbuf(stdin, NULL);
    Fl_Window win(400,400,"FLTK Mesh Viewer");

    meshw.resize(nwin);
    int subh = ( win.h()-10 )/nwin;
    for (int i = 0; i < nwin; ++i) {
        meshw[i] = new MyGLWindow(10, 10+i*subh, win.w()-20, subh-10);
        meshw[i]->end();
    }

    Fl::add_fd(fileno(stdin), handle_stdin, NULL);
    Fl::add_timeout(0.1, show_next_frame_cb);
    Fl::gl_visual(FL_RGB);
    win.resizable(win);
    win.end();
    win.show();
    return(Fl::run());
}

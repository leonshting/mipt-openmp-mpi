#include "lifeviewer.h"

// Cube
void LifeViewer::primCube(const coord center, const float a)
{
    glBegin(GL_QUADS);
    glVertex3f(center.X+a,center.Y+a,center.Z+a);
    glVertex3f(center.X+a,center.Y-a,center.Z+a);
    glVertex3f(center.X-a,center.Y-a,center.Z+a);
    glVertex3f(center.X-a,center.Y+a,center.Z+a);

    glVertex3f(center.X+a,center.Y+a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y+a,center.Z-a);

    glVertex3f(center.X+a,center.Y+a,center.Z+a);
    glVertex3f(center.X+a,center.Y+a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z+a);

    glVertex3f(center.X-a,center.Y+a,center.Z+a);
    glVertex3f(center.X-a,center.Y+a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z+a);

    glVertex3f(center.X+a,center.Y+a,center.Z+a);
    glVertex3f(center.X+a,center.Y+a,center.Z-a);
    glVertex3f(center.X-a,center.Y+a,center.Z-a);
    glVertex3f(center.X-a,center.Y+a,center.Z+a);

    glVertex3f(center.X+a,center.Y-a,center.Z+a);
    glVertex3f(center.X+a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z+a);
    glEnd();
}

// Cube wireframe
void LifeViewer::primCubeWF(const coord center, const float a)
{
    glBegin(GL_LINE_LOOP);
    glVertex3f(center.X+a,center.Y+a,center.Z+a);
    glVertex3f(center.X+a,center.Y+a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z+a);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex3f(center.X-a,center.Y+a,center.Z+a);
    glVertex3f(center.X-a,center.Y+a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z+a);
    glEnd();
    glBegin(GL_LINES);
    glVertex3f(center.X+a,center.Y+a,center.Z+a);
    glVertex3f(center.X-a,center.Y+a,center.Z+a);
    glVertex3f(center.X+a,center.Y+a,center.Z-a);
    glVertex3f(center.X-a,center.Y+a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z-a);
    glVertex3f(center.X-a,center.Y-a,center.Z-a);
    glVertex3f(center.X+a,center.Y-a,center.Z+a);
    glVertex3f(center.X-a,center.Y-a,center.Z+a);
    glEnd();
}

void LifeViewer::primCubeWF(const coord p1, const coord p2)
{
    glBegin(GL_LINE_LOOP);
    glVertex3i(p1.X,p1.Y,p1.Z);
    glVertex3i(p2.X,p1.Y,p1.Z);
    glVertex3i(p2.X,p2.Y,p1.Z);
    glVertex3i(p1.X,p2.Y,p1.Z);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex3i(p1.X,p1.Y,p2.Z);
    glVertex3i(p2.X,p1.Y,p2.Z);
    glVertex3i(p2.X,p2.Y,p2.Z);
    glVertex3i(p1.X,p2.Y,p2.Z);
    glEnd();
    glBegin(GL_LINES);
    glVertex3i(p1.X,p1.Y,p1.Z);
    glVertex3i(p1.X,p1.Y,p2.Z);
    glVertex3i(p2.X,p1.Y,p1.Z);
    glVertex3i(p2.X,p1.Y,p2.Z);
    glVertex3i(p2.X,p2.Y,p1.Z);
    glVertex3i(p2.X,p2.Y,p2.Z);
    glVertex3i(p1.X,p2.Y,p1.Z);
    glVertex3i(p1.X,p2.Y,p2.Z);
    glEnd();
}

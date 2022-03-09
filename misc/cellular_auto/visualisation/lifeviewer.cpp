#include <QFileDialog>

#include "lifeviewer.h"
#include "readerxyz.h"
#include "readerwxyz.h"

#ifndef min
#define min std::min
#endif
#ifndef max
#define max std::max
#endif

void LifeViewer::gotoFrame(int n)
{
    if (n<0)
        n = data_.size()-1;
    if (n>static_cast<int>(data_.size()-1))
        n = 0;
    curFrame_ = n;
    emit frameChanged(n);
    QGLViewer::repaint();
}

void LifeViewer::init()
{
    //restoreStateFromFile();
    glDisable(GL_LIGHTING);
}

void LifeViewer::draw()
{
    if (static_cast<int>(data_.size())<=curFrame_ || curFrame_<0) return;
    glLineWidth(1.f);
    glColor3f(.3f,.3f,.3f);
    primCubeWF(cmin,cmax);
    frame::iterator it;
    QColor color;
    switch(drawingMode_)
    {
        case WCUBES:
            glLineWidth(3.f);
            for(it = data_[curFrame_].begin(); it!=data_[curFrame_].end(); ++it)
            {
                color = colors_.value(it->W,color0_);
                glColor3f(color.redF(),color.greenF(),color.blueF());
                primCubeWF(*it,.4f);
            }
            glLineWidth(1.f);
            break;
        case BCUBES:
            glLineWidth(2.f);
            for(it = data_[curFrame_].begin(); it!=data_[curFrame_].end(); ++it)
            {
                color = colors_.value(it->W,color0_);
                glColor3f(color.redF()/2,color.greenF()/2,color.blueF()/2);
                primCubeWF(*it,.4f);
            }
            glLineWidth(1.f);
            //fall thru
        case CUBES:
            for(it = data_[curFrame_].begin(); it!=data_[curFrame_].end(); ++it)
            {
                color = colors_.value(it->W,color0_);
                glColor3f(color.redF(),color.greenF(),color.blueF());
                primCube(*it,.4f);
            }
            break;
        case POINTS:
        default:
            glPointSize(3.0);
            glBegin(GL_POINTS);
            for(it = data_[curFrame_].begin(); it!=data_[curFrame_].end(); ++it)
            {
                color = colors_.value(it->W,color0_);
                glColor3f(color.redF(),color.greenF(),color.blueF());
                glVertex3i(it->X,it->Y,it->Z);
            }
            glEnd();
            break;
    }
    if (textShown_)
    {
        glColor3f(1.f,1.f,1.f);
        drawText(1, this->height()-1,QString("\"")+filename_+"\" ["+QString::number(curFrame_)+"]",QFont("Arial",10));
    }
}

void LifeViewer::animate()
{
    nextFrame();
}

QString LifeViewer::helpString() const
{
    QString text("<h2>Visualisator for 3D version of Conway's Game of Life</h2>\n");
    text += "It might even work";
    return text;
}

void LifeViewer::loadData()
{
    loadFromFile( QFileDialog::getOpenFileName(parent_,QString("Open file..."),"",QString("All files (*.*);;WXYZ files (*.wxyz);;XYZ files (*.xyz)")) );
    setRadiusAndCenter();
    gotoFrame(0);
    emit frameMaxChanged(data_.size()-1);
    showEntireScene();
	//std::cout<<"Loading finished!"<<std::endl;
}

void LifeViewer::loadFromFile(QString filename)
{
    if (filename=="") return;
    Reader *r;
    if (filename.endsWith(ReaderXYZ::filetype()))
        r = new ReaderXYZ(filename);
    else if (filename.endsWith(ReaderWXYZ::filetype()))
        r = new ReaderWXYZ(filename);
	else // default
		r = new ReaderWXYZ(filename);

    data_ = r->read();
    delete r;
    filename_ = filename;
}

void LifeViewer::setRadiusAndCenter()
{
	if (data_.size() <= 0) return;
	bool set = false;
	for(int i=0; i!=data_.size(); ++i)
	{
		if (!set && data_[i].size()>0)
		{
			cmax.X = cmin.X =data_[i][0].X;	cmax.Y = cmin.Y =data_[i][0].Y; cmax.Z = cmin.Z =data_[i][0].Z;
			set = true;
		}
		for(int j=0; j!=data_[i].size(); ++j)
        {
			cmin.X = min(cmin.X,data_[i][j].X); cmin.Y = min(cmin.Y,data_[i][j].Y); cmin.Z = min(cmin.Z,data_[i][j].Z);
            cmax.X = max(cmax.X,data_[i][j].X); cmax.Y = max(cmax.Y,data_[i][j].Y); cmax.Z = max(cmax.Z,data_[i][j].Z);
        }
	}
    setSceneBoundingBox(qglviewer::Vec(cmin.X,cmin.Y,cmin.Z),qglviewer::Vec(cmax.X,cmax.Y,cmax.Z));
}


#include "readerwxyz.h"

ReaderWXYZ::ReaderWXYZ(QString filename)
{
    file = new QFile(filename);
    file->open(QIODevice::ReadOnly | QIODevice::Text);
}
ReaderWXYZ::~ReaderWXYZ()
{
    delete file;
}

frames ReaderWXYZ::read()
{
    QTextStream stream(file);
    frames a;
    int w,x,y,z;
    QString t;
    while (!stream.atEnd())
    {
        a.push_back(frame());
        for(stream>>w>>x>>y>>z;w!=0;stream>>w>>x>>y>>z)
        {
            a[a.size()-1].push_back( coord(std::pair<int,int>(w,x), std::pair<int,int>(y,z)) );
        }
    }
    if (a[a.size()-1].size() == 0)
        a.erase(--a.end());
    return a;
}

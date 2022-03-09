#include "readerxyz.h"

ReaderXYZ::ReaderXYZ(QString filename)
{
    file = new QFile(filename);
    file->open(QIODevice::ReadOnly | QIODevice::Text);
}
ReaderXYZ::~ReaderXYZ()
{
    delete file;
}

frames ReaderXYZ::read()
{
    QTextStream stream(file);
    frames a;
    int l,w,x,y,z;
    QString t;
    while (!stream.atEnd())
    {
        stream >> l;
        if (l!=0)
        {
            a.push_back(frame(l));
            for(int i=0;i!=l;++i)
            {
                stream >> w >> x >> y >> z;
                a[a.size()-1][i] = coord(std::pair<int,int>(w,x), std::pair<int,int>(y,z));
            }
        }
    }
    return a;
}

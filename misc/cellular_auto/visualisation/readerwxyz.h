#ifndef READERWXYZ_H
#define READERWXYZ_H

#include <QFile>
#include <QTextStream>
#include "reader.h"

class ReaderWXYZ : public Reader
{
public:
    ReaderWXYZ(QString filename);
    virtual frames read();
    ~ReaderWXYZ();
    static const QString filetype() {return ".wxyz";}
private:
    QFile *file;
};

#endif // READERWXYZ_H

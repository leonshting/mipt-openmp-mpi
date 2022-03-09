#ifndef READERXYZ_H
#define READERXYZ_H

#include <QFile>
#include <QTextStream>
#include "reader.h"

class ReaderXYZ : public Reader
{
public:
    ReaderXYZ(QString filename);
    virtual frames read();
    ~ReaderXYZ();
    static const QString filetype() {return ".xyz";}
private:
    QFile *file;
};

#endif // READERXYZ_H

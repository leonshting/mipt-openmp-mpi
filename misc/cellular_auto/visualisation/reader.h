#ifndef READER_H
#define READER_H

#include <QFile>
#include "common.h"

class Reader
{
public:
    virtual frames read() = 0;
    static const QString filetype() {return ".";}
};

#endif // READER_H

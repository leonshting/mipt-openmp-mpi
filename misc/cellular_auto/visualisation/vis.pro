# -------------------------------------------------
# Project created by QtCreator 2010-07-09T17:06:30
# -------------------------------------------------
QT += opengl \
    xml
TARGET = vis
TEMPLATE = app
SOURCES += main.cpp \
    mainwindow.cpp \
    lifeviewer.cpp \
    readerxyz.cpp \
    readerwxyz.cpp \
    GLPrimitives.cpp
HEADERS += mainwindow.h \
    lifeviewer.h \
    readerxyz.h \
    common.h \
    readerwxyz.h \
    reader.h
FORMS += mainwindow.ui
unix:LIBS += -lqglviewer-qt4
win32:LIBS += -lQGLViewer2 -L"C:\tempgrav\libQGLViewer-2.3.6\QGLViewer\release"
win32:INCLUDEPATH += "C:\Program Files\libQGLViewer\"
win32:DEFINES += WIN32
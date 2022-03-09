#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QColor>
#include <QColorDialog>
#include <QSettings>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

protected:
    void changeEvent(QEvent *e);

private:
    Ui::MainWindow *ui;
    QSettings *settings;
    void loadSettings();
    void saveSettings();

private slots:
    void on_sbMax_valueChanged( int );
    void makeMovie();
    void setColorClicked();
    void setColorWChanged( int );
};

#endif // MAINWINDOW_H

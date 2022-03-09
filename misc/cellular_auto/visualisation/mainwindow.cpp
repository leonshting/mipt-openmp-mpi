#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    ui->cbType->addItem("Points",ui->viewer->POINTS);
    ui->cbType->addItem("Cubes",ui->viewer->CUBES);
    ui->cbType->addItem("BCubes",ui->viewer->BCUBES);
    ui->cbType->addItem("SLines",ui->viewer->WCUBES);

    settings = new QSettings("MIPT","Life3D Visualiser",this);
    loadSettings();

}

MainWindow::~MainWindow()
{
    delete ui;
    saveSettings();
    delete settings;
}

void MainWindow::setColorClicked()
{
    int w = ui->sbColorW->value();
    ui->viewer->setColor(w, QColorDialog::getColor(ui->viewer->getColor(w),this));
    setColorWChanged(w);
}
void MainWindow::setColorWChanged( int w )
{
    QColor c = ui->viewer->getColor(w);
    QString styleSheet = QString("* { background-color: rgb(")+ QString::number(c.red()) +","+ QString::number(c.green()) +","+ QString::number(c.blue()) +") }";
    ui->pColorC->setStyleSheet(styleSheet);
}

void MainWindow::changeEvent(QEvent *e)
{
    QMainWindow::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        ui->retranslateUi(this);
        break;
    default:
        break;
    }
}

void MainWindow::makeMovie()
{
    for(int i=0;i<=ui->viewer->frameMax();++i)
    {
        ui->viewer->gotoFrame(i);
        ui->viewer->saveSnapshot();
    }
}

void MainWindow::on_sbMax_valueChanged(int m)
{
    ui->sbNum->setMaximum(m);
}

void MainWindow::saveSettings()
{
    settings->setValue("animation/delay",   ui->sbAnimDelay->value());
    settings->setValue("display/fps",       ui->cbFPS->isChecked());
    settings->setValue("display/grid",      ui->cbGrid->isChecked());
    settings->setValue("display/text",      ui->cbText->isChecked());
    settings->setValue("display/type",      ui->cbType->currentIndex());
    QHash<int,QColor> qh = ui->viewer->getColors();
    for (QHash<int,QColor>::iterator it = qh.begin(); it!=qh.end(); ++it)
    {
        settings->setValue(QString("colors/")+QString::number(it.key()),it.value());
    }
    settings->sync();
}

void MainWindow::loadSettings()
{
    ui->sbAnimDelay->setValue( settings->value("animation/delay",1000).toInt() );
    ui->cbFPS->setChecked( settings->value("display/fps",false).toBool() );
    ui->cbGrid->setChecked( settings->value("display/grid",false).toBool() );
    ui->cbText->setChecked( settings->value("display/text",false).toBool() );
    ui->cbType->setCurrentIndex( settings->value("display/type",2).toInt() );

    settings->beginGroup("colors");
    QStringList ws = settings->allKeys();
    for (QStringList::iterator it = ws.begin(); it!=ws.end(); ++it)
    {
        ui->viewer->setColor( (*it).toInt() , settings->value(*it).value<QColor>() );
    }
    settings->endGroup();
    setColorWChanged(ui->sbColorW->value());

    ui->viewer->setAnimationPeriod(ui->sbAnimDelay->value());
    ui->viewer->setFPSIsDisplayed(ui->cbFPS->isChecked());
    ui->viewer->setShowText(ui->cbText->isChecked());
    ui->viewer->setGridIsDrawn(ui->cbGrid->isChecked());
    ui->viewer->setDrawingMode(ui->cbType->currentIndex());
}

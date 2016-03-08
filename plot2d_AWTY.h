#ifndef PLOT2D_AWTY_H
#define PLOT2D_AWTY_H

#include <QDialog>
#include "wstring.h"
#include "qcustomplot.h"
#include <string>

namespace Ui {
class Plot2D;
}

class Plot2D_AWTY : public QDialog
{
    Q_OBJECT

public:
    explicit Plot2D_AWTY(QString xlabel, QString ylabel, double xstart, double xend, double ystart, double yend, QString title, QWidget *parent = 0, bool newlambda = 0);
    void initialization_AWTY(QString xlabel, QString ylabel, double xstart, double xend, double ystart, double yend, QString title, bool newlambda);
    ~Plot2D_AWTY();

    void Plot2DCurves_AWTY(QVector<double> x, QVector<double> y, QString qtname, QColor color);
    void Clear2DCurves_AWTY();

private slots:
    void titleDoubleClick(QMouseEvent *event, QCPPlotTitle *title);
    void axisLabelDoubleClick(QCPAxis* axis, QCPAxis::SelectablePart part);
    void legendDoubleClick(QCPLegend* legend, QCPAbstractLegendItem* item);
    void selectionChanged();
    void mousePress();
    void mouseWheel();
    void removeSelectedGraph();
    void removeAllGraphs();
    void contextMenuRequest(QPoint pos);
    void moveLegend();
    void showPointToolTip(QMouseEvent *event);


private:
    Ui::Plot2D *ui;
};

#endif // PLOT2D_AWTY_H

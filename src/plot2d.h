#ifndef PLOT2D_H
#define PLOT2D_H

#include <QDialog>
#include "wstring.h"
#include "qcustomplot.h"
#include <string>

namespace Ui {
class Plot2D;
}

class Plot2D : public QDialog
{
    Q_OBJECT

public:
//    explicit Plot2D(double **com_info, int length, QWidget *parent = 0);
    explicit Plot2D(QString xlabel, QString ylabel, double xstart, double xend, double ystart, double yend, QString title, QWidget *parent = 0, bool newlambda = true);
    void initialization(QString xlabel, QString ylabel, double xstart, double xend, double ystart, double yend, QString title, bool newlambda = true);
    ~Plot2D();

    void Plot2DCurves(QVector<double> x, QVector<double> y, QString qtname, QColor color);
    void Clear2DCurves();

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
//    void graphClicked(QCPAbstractPlottable *plottable);

    void on_pushButton_clicked();

signals:
    void sendnewlambda(double, Plot2D *, int);
private:
    Ui::Plot2D *ui;
    int plottype; // 1: community cova, 2: community affi
};

#endif // PLOT2D_H

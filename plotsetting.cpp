#include "plotsetting.h"
#include "ui_plotsetting.h"

plotsetting::plotsetting(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::plotsetting)
{
    ui->setupUi(this);
}

plotsetting::~plotsetting()
{
    delete ui;
}

#ifndef PLOTSETTING_H
#define PLOTSETTING_H

#include <QDialog>

namespace Ui {
class plotsetting;
}

class plotsetting : public QDialog
{
    Q_OBJECT

public:
    explicit plotsetting(QWidget *parent = 0);
    ~plotsetting();

private:
    Ui::plotsetting *ui;
};

#endif // PLOTSETTING_H

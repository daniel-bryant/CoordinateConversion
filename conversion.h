#ifndef CONVERSION_H
#define CONVERSION_H

#include <QStringList>

class Conversion
{
public:
    Conversion();
    static QStringList GeogToUTM(double, double);
    static QStringList UTMtoGeog(double, double, double);
};

#endif // CONVERSION_H

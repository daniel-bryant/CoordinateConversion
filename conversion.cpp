#include "conversion.h"

Conversion::Conversion()
{
}

QStringList Conversion::GeogToUTM(double latIn, double lonIn)
{
    //Convert Latitude and Longitude to UTM

    //Declarations
    const double PI = 4.0*atan(1.0);
    //Symbols as used in USGS PP 1395: Map Projections - A Working Manual
    double DatumEqRad[14] = {6378137.0,6378137.0,6378137.0,6378135.0,6378160.0,6378245.0,6378206.4,
        6378388.0,6378388.0,6378249.1,6378206.4,6377563.4,6377397.2,6377276.3};
    double DatumFlat[14] = {298.2572236, 298.2572236, 298.2572215,	298.2597208, 298.2497323, 298.2997381, 294.9786982,
        296.9993621, 296.9993621, 293.4660167, 294.9786982, 299.3247788, 299.1527052, 300.8021499};
    int Item = 0;//Default
    double k0 = 0.9996;//scale on central meridian
    double a = DatumEqRad[Item];//equatorial radius, meters.
    double f = 1/DatumFlat[Item];//polar flattening.
    double b = a*(1-f);//polar axis.
    double e = sqrt(1 - b*b/a*a);//eccentricity
    double drad = PI/180;//Convert degrees to radians)
    double latd = 0;//latitude in degrees
    double phi = 0;//latitude (north +, south -), but uses phi in reference
    double N = 0;
    double T = pow(tan(phi),2);
    double C = pow(e*cos(phi),2);
    double lngd = 0;//longitude in degrees
    double M = 0;//M requires calculation
    double x = 0;//x coordinate
    double y = 0;//y coordinate
    double utmz = 30;//utm zone
    double zcm = 0;//zone central meridian
    //end of Declarations

    k0 = 0.9996;//scale on central meridian
    b = a*(1-f);//polar axis.
    e = sqrt(1 - (b/a)*(b/a));//eccentricity
    //Input Geographic Coordinates
    //Decimal Degree Option
    double latd0 = latIn;
    double lngd0 = lonIn;

    lngd=lngd0;
    latd=latd0;

    if(isnan(latd)|| isnan(lngd)){
        qDebug("Non-Numeric Input Value");
    }
    if(latd <-90 || latd> 90){
        qDebug("Latitude must be between -90 and 90");
    }
    if(lngd <-180 || lngd > 180){
        qDebug("Latitude must be between -180 and 180");
    }

    phi = latd*drad;//Convert latitude to radians
    utmz = 1 + floor((lngd+180)/6);//calculate utm zone

    zcm = 3 + 6*(utmz-1) - 180;//Central meridian of zone
    //Calculate Intermediate Terms
    double esq = (1 - (b/a)*(b/a));//e squared for use in expansions
    double e0sq = e*e/(1-e*e);// e0 squared - always even powers
    N = a/sqrt(1-pow(e*sin(phi),2));
    T = pow(tan(phi),2);
    C = e0sq*pow(cos(phi),2);
    double A = (lngd-zcm)*drad*cos(phi);
    //Calculate M
    M = phi*(1 - esq*(1/4 + esq*(3/64 + 5*esq/256)));
    M = M - sin(2*phi)*(esq*(3/8 + esq*(3/32 + 45*esq/1024)));
    M = M + sin(4*phi)*(esq*esq*(15/256 + esq*45/1024));
    M = M - sin(6*phi)*(esq*esq*esq*(35/3072));
    M = M*a;//Arc length along standard meridian
    double M0 = 0;//M0 is M for some origin latitude other than zero. Not needed for standard UTM
    //Calculate UTM Values
    x = k0*N*A*(1 + A*A*((1-T+C)/6 + A*A*(5 - 18*T + T*T + 72*C -58*e0sq)/120));//Easting relative to CM
    x=x+500000;//Easting standard
    y = k0*(M - M0 + N*tan(phi)*(A*A*(1/2 + A*A*((5 - T + 9*C + 4*C*C)/24 + A*A*(61 - 58*T + T*T + 600*C - 330*e0sq)/720))));//Northing from equator
    if (y < 0){y = 10000000+y;}
    QStringList utmOut;
    utmOut.append(QString::number(utmz)); //zone
    utmOut.append(QString::number(y)); //northing
    utmOut.append(QString::number(x)); //easting

    return utmOut;

}//close Geog to UTM

QStringList Conversion::UTMtoGeog(double zoneIn, double northingIn, double eastingIn){
    //Convert UTM Coordinates to Geographic

    //Declarations
    const double PI = 4.0*atan(1.0);
    //Symbols as used in USGS PP 1395: Map Projections - A Working Manual
    double DatumEqRad[14] = {6378137.0,6378137.0,6378137.0,6378135.0,6378160.0,6378245.0,6378206.4,
        6378388.0,6378388.0,6378249.1,6378206.4,6377563.4,6377397.2,6377276.3};
    double DatumFlat[14] = {298.2572236, 298.2572236, 298.2572215,	298.2597208, 298.2497323, 298.2997381, 294.9786982,
        296.9993621, 296.9993621, 293.4660167, 294.9786982, 299.3247788, 299.1527052, 300.8021499};
    int Item = 0;//Default
    double k0 = 0.9996;//scale on central meridian
    double a = DatumEqRad[Item];//equatorial radius, meters.
    double f = 1/DatumFlat[Item];//polar flattening.
    double b = a*(1-f);//polar axis.
    double e = sqrt(1 - b*b/a*a);//eccentricity
    double drad = PI/180;//Convert degrees to radians)
    double phi = 0;//latitude (north +, south -), but uses phi in reference
    double lng = 0;//Longitude (e = +, w = -) - can't use long - reserved word
    double lngd = 0;//longitude in degrees
    double M = 0;//M requires calculation
    double x = 0;//x coordinate
    double y = 0;//y coordinate
    double utmz = 30;//utm zone
    double zcm = 0;//zone central meridian
    //end of Declarations

    k0 = 0.9996;//scale on central meridian
    b = a*(1-f);//polar axis.
    e = sqrt(1 - (b/a)*(b/a));//eccentricity
    double esq = (1 - (b/a)*(b/a));//e squared for use in expansions
    double e0sq = e*e/(1-e*e);// e0 squared - always even powers
    x = eastingIn;
    if (x<160000 || x>840000){qDebug("Outside permissible range of easting values \n Results may be unreliable \n Use with caution");}
    y = northingIn;
    if (y<0){qDebug("Negative values not allowed \n Results may be unreliable \n Use with caution");}
    if (y>10000000){qDebug("Northing may not exceed 10,000,000 \n Results may be unreliable \n Use with caution");}
    utmz = zoneIn;
    zcm = 3 + 6*(utmz-1) - 180;//Central meridian of zone
    double e1 = (1 - sqrt(1 - e*e))/(1 + sqrt(1 - e*e));//Called e1 in USGS PP 1395 also
    double M0 = 0;//In case origin other than zero lat - not needed for standard UTM
    M = M0 + y/k0;//Arc length along standard meridian.
    double mu = M/(a*(1 - esq*(1/4 + esq*(3/64 + 5*esq/256))));
    double phi1 = mu + e1*(3/2 - 27*e1*e1/32)*sin(2*mu) + e1*e1*(21/16 -55*e1*e1/32)*sin(4*mu);//Footprint Latitude
    phi1 = phi1 + e1*e1*e1*(sin(6*mu)*151/96 + e1*sin(8*mu)*1097/512);
    double C1 = e0sq*pow(cos(phi1),2);
    double T1 = pow(tan(phi1),2);
    double N1 = a/sqrt(1-pow(e*sin(phi1),2));
    double R1 = N1*(1-e*e)/(1-pow(e*sin(phi1),2));
    double D = (x-500000)/(N1*k0);
    phi = (D*D)*(1/2 - D*D*(5 + 3*T1 + 10*C1 - 4*C1*C1 - 9*e0sq)/24);
    phi = phi + pow(D,6)*(61 + 90*T1 + 298*C1 + 45*T1*T1 -252*e0sq - 3*C1*C1)/720;
    phi = phi1 - (N1*tan(phi1)/R1)*phi;

    QStringList geoOut;
    //Output Latitude
    geoOut.append(QString::number( (1000000*phi/drad)/1000000 ));

    //Longitude
    lng = D*(1 + D*D*((-1 -2*T1 -C1)/6 + D*D*(5 - 2*C1 + 28*T1 - 3*C1*C1 +8*e0sq + 24*T1*T1)/120))/cos(phi1);
    lngd = zcm+lng/drad;

    //Output Longitude
    geoOut.append(QString::number( floor(1000000*lngd)/1000000 ));

    return geoOut;

}//End UTM to Geog

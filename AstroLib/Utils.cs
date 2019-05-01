using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


namespace AstroLib {

    public static class Math2 {
        public static double Sec(double rad) {
            return 1 / Math.Cos(rad);
        }

        public static double Cosec(double rad) {
            return 1 / Math.Sin(rad);
        }

        public static double Cot(double rad) {
            return 1 / Math.Tan(rad);
        }

    }
   

    public static class MyExtensions {
        public static double ToRad(this Double deg) {
            return (deg * Math.PI / 180.0);
        }

        public static double ToDeg(this Double rad) {
            return (rad / Math.PI * 180.0);
        }

         }

    public static class Utils {
                                       
        //Return the current julian date
        public static double jd() {     
            
            TimeSpan _TimeSpan = (DateTime.UtcNow - new DateTime(1970, 1, 1, 0, 0, 0));
            
            return ((_TimeSpan.TotalSeconds) / 86400) + 2440587.5;
        }

        public static double lst(double JD, double T, Location loc) {

            double theta = (280.46061837 + (360.98564736629 * (JD - 2451545)) + (0.000387933 * T * T) - ((T * T * T) / 38710000)).ToRad();
            
            return theta- loc.lon;
        }

        //Given an angle in decimal degrees, return the degree component with +/- sign and leading zeros
        public static string dd(double deg) {

            string degs;

            if (deg >= 0) {
                degs = Math.Floor(deg).ToString("+00");
            } else if (deg > -1) {
                degs = Math.Ceiling(deg).ToString("-00");
            } else {
                degs = Math.Ceiling(deg).ToString("00");
            }

            return degs;
        }

        //Same as above but no positive sign
        public static string dd1(double deg) {
            string degs;

            if (deg >= 0) {
                degs = Math.Floor(deg).ToString("00");
            } else if (deg > -1) {
                degs = Math.Ceiling(deg).ToString("-00");
            } else {
                degs = Math.Ceiling(deg).ToString("00");
            }

            return degs;
        }

        //Same above but no sign at all (for hours)
        public static string hh(double deg) {

            return Math.Floor(deg).ToString("00");
        }

        //Given an angle in decimal degrees, return the minutes component with leading zeros
        public static string dm(double deg) {

            deg = Math.Abs(deg);

            return Math.Floor((deg - Math.Floor(deg)) * 60).ToString("00");
        }

        //Given an angle in decimal degrees, return the seconds component with leading zeros
        public static string ds(double deg) {

            deg = Math.Abs(deg);

            double mini = ((deg - Math.Floor(deg)) * 60);

            return Math.Round((mini - Math.Floor(mini)) * 60).ToString("00");
        }

        public static CoordsString EquatorialToDMS(EquatorialCoords coords) {

            CoordsString StringOut = new CoordsString();

            double RAH = rad2deg(coords.ra) / 15;
            double DecD = rad2deg(coords.dec);

            StringOut.ra = hh(RAH) + "h " + dm(RAH) + "m " + ds(RAH) + "s";
            StringOut.dec = dd(DecD) + (Char)176 + " " + dm(DecD) + (Char)39 + " " + ds(DecD) + (Char)34;

            return StringOut;
        }      

        //return a value in the range 0<x<360
        public static double quad(double inval) {

            if (inval < 0) {
                return inval + ((0 - (Math.Floor(inval / 360))) * 360);
            } else {
                if (inval > 360) {
                    return inval - ((Math.Floor(inval / 360)) * 360);
                } else {
                    return inval;
                }
            }
        }

        //Convert degrees to radians
        public static  double deg2rad(double deg) {

            return (deg * Math.PI / 180.0);
        }

        //Convert radians to degrees
        public static  double rad2deg(double rad) {

            return (rad / Math.PI * 180.0);
        }

        // Calculate the secant of angle, in radians.
        public static double Sec(double angle) {

            return 1.0 / Math.Cos(angle);
        }

        //Convert Alt/Az to RA/Dec. Input and output in rads
        public static EquatorialCoords HorizontalToEquatorial(HorizontalCoords coords, Location loc, double T) {

            EquatorialCoords coordsout = new EquatorialCoords();

            double HA = Math.Atan2(Math.Sin(coords.az + Math.PI), Math.Cos(coords.az + Math.PI) * Math.Sin(loc.lat) + Math.Tan(coords.alt) * Math.Cos(loc.lat));

            double theta = deg2rad(280.46061837 + (360.98564736629 * T * 36525) + (0.000387933 * T * T) - ((T * T * T) / 38710000));

            coordsout.ra = theta - loc.lon - HA;

            coordsout.dec = Math.Asin(Math.Sin(loc.lat) * Math.Sin(coords.alt) - Math.Cos(loc.lat) * Math.Cos(coords.alt) * Math.Cos(coords.az + Math.PI));

            return coordsout;
        }

        public static HorizontalCoords EquatorialToHorizontal(EquatorialCoords coords, Location loc, double T) {

            HorizontalCoords coordsout = new HorizontalCoords();

            double theta = deg2rad(280.46061837 + (360.98564736629 * T * 36525) + (0.000387933 * T * T) - ((T * T * T) / 38710000));

            double HA = theta - loc.lon - coords.ra;

            coordsout.az = Math.Atan2(Math.Sin(HA), Math.Cos(HA) * Math.Sin(loc.lat) - Math.Tan(coords.dec) * Math.Cos(loc.lat)) + Math.PI;

            coordsout.alt = Math.Asin(Math.Sin(loc.lat) * Math.Sin(coords.dec) + Math.Cos(loc.lat) * Math.Cos(coords.dec) * Math.Cos(HA));

            return coordsout;
        }


        //Precess coords. Input and output in rads. JD0 is Julian date of starting epoch, and JD final epoch
        public static EquatorialCoords Precess(EquatorialCoords coords, double JD0, double JD) {

            EquatorialCoords coordsout = new EquatorialCoords();

            double T = (JD0 - 2451545.0) / 36525;
            double t = (JD - JD0) / 36525;

            double eta = deg2rad(((2306.2181 + 1.39656 * T - 0.000139 * T * T) * t + (0.30188 - 0.000344 * T) * t * t + 0.017998 * t * t * t) / 3600);
            double zeta = deg2rad(((2306.2181 + 1.39656 * T - 0.000139 * T * T) * t + (1.09468 + 0.000066 * T) * t * t + 0.018203 * t * t * t) / 3600);
            double theta = deg2rad(((2004.3109 - 0.8533 * T - 0.000217 * T * T) * t - (0.42665 + 0.000217 * T) * t * t - 0.041833 * t * t * t) / 3600);

            double A = Math.Cos(coords.dec) * Math.Sin(coords.ra + eta);
            double B = Math.Cos(theta) * Math.Cos(coords.dec) * Math.Cos(coords.ra + eta) - Math.Sin(theta) * Math.Sin(coords.dec);
            double C = Math.Sin(theta) * Math.Cos(coords.dec) * Math.Cos(coords.ra + eta) + Math.Cos(theta) * Math.Sin(coords.dec);

            coordsout.ra = Math.Atan2(A, B) + zeta;

            if (coords.dec > 1.5) { //~86 degres
                coordsout.dec = Math.Acos(Math.Sqrt(A * A + B * B));
            } else {
                coordsout.dec = Math.Asin(C);
            }
            return coordsout;
        }


        public static EquatorialCoords J2000ToB1950(EquatorialCoords coords) {

            double JD0 = 2451545.0; //j2000
            double JD = 2433282.4235; //b1950

            return Precess(coords, JD0, JD);
        }



        public static EquatorialCoords NowToJ2000(EquatorialCoords coords, double JD0) {

            double JD = 2451545.0; //j2000

            return Precess(coords, JD0, JD);
        }



        public static EquatorialCoords J2000ToNow(EquatorialCoords coords, double JD) {

            double JD0 = 2451545.0;

            return Precess(coords, JD0, JD);
        }



        public static GalacticCoords B1950ToGalactic(EquatorialCoords coords) {

            GalacticCoords coordsout = new GalacticCoords();

            double x = Math.Atan2(Math.Sin(deg2rad(192.25) - coords.ra), Math.Cos(deg2rad(192.25) - coords.ra) * Math.Sin(deg2rad(27.4)) - Math.Tan(coords.dec) * Math.Cos(deg2rad(27.4)));
            coordsout.l = deg2rad(quad(303 - rad2deg(x)));

            double sinB = Math.Sin(coords.dec) * Math.Sin(deg2rad(27.4)) + Math.Cos(coords.dec) * Math.Cos(deg2rad(27.4)) * Math.Cos(deg2rad(192.25) - coords.ra);
            coordsout.b = Math.Asin(sinB);

            return coordsout;
        }


        public static EquatorialCoords GalacticToB1950(GalacticCoords coords) {

            EquatorialCoords coordsout = new EquatorialCoords();

            double y = Math.Atan2(Math.Sin(coords.l - deg2rad(123)), Math.Cos(coords.l - deg2rad(123)) * Math.Sin(deg2rad(27.4)) - Math.Tan(coords.b) * Math.Cos(deg2rad(27.4)));
            coordsout.ra = deg2rad(quad(12.25 + rad2deg(y)));

            double sind = Math.Sin(coords.b) * Math.Sin(deg2rad(27.4)) + Math.Cos(coords.b) * Math.Cos(deg2rad(27.4)) * Math.Cos(coords.l - deg2rad(123));
            coordsout.dec = Math.Asin(sind);

            return coordsout;
        }

        //Convert Julian date to t= Julian millenia since J2000.0
        public static double JDtot(double JD) {
            return (JD - 2451545) / 365250;
        }

        //Convert Julian date to T= Julian centuries since J2000.0
        public static double JDtoT(double JD) {
            return (JD - 2451545) / 36525;
        }

        //Convert a dateTime object to Julian date
        public static double DateTimetoJD(this DateTime date) {
            return date.ToOADate() + 2415018.5;
        }

        //Heliocentric ecliptical to Geocentric ecliptical
        public static EclipticalCoords HelToGeo(EclipticalCoords Earth, EclipticalCoords Planet) {

            EclipticalCoords geo = new EclipticalCoords();

            double x = Planet.r * Math.Cos(Planet.b) * Math.Cos(Planet.l) - Earth.r * Math.Cos(Earth.b) * Math.Cos(Earth.l);
            double y = Planet.r * Math.Cos(Planet.b) * Math.Sin(Planet.l) - Earth.r * Math.Cos(Earth.b) * Math.Sin(Earth.l);
            double z = Planet.r * Math.Sin(Planet.b) - Earth.r * Math.Sin(Earth.b);

            geo.l = Math.Atan2(y, x);
            geo.b = Math.Atan2(z, Math.Sqrt(x * x + y * y));
            geo.r = Math.Sqrt(x * x + y * y + z * z);

            return geo;
        }

        //Ecliptical to Celestial coordinate conversion. Epsilon (obliquity of the ecliptic) in rads
        public static EquatorialCoords EclipticalToCelestial(EclipticalCoords eclipt, NutObl nutobl) {

            EquatorialCoords coords = new EquatorialCoords();

            coords.ra = Math.Atan2((Math.Sin(eclipt.l) * Math.Cos(nutobl.epsilon)) - (Math.Tan(eclipt.b) * Math.Sin(nutobl.epsilon)), Math.Cos(eclipt.l));
            coords.dec = Math.Asin((Math.Sin(eclipt.b) * Math.Cos(nutobl.epsilon)) + (Math.Cos(eclipt.b) * Math.Sin(nutobl.epsilon) * Math.Sin(eclipt.l)));

            return coords;
        }

        public static EclipticalCoords EquatorialToEcliptical(EquatorialCoords celest, NutObl nutobl) {

            EclipticalCoords coords = new EclipticalCoords();

            coords.l = Math.Atan2((Math.Sin(celest.ra) * Math.Cos(nutobl.epsilon)) + (Math.Tan(celest.dec) * Math.Sin(nutobl.epsilon)), Math.Cos(celest.ra));
            coords.b = Math.Asin((Math.Sin(celest.dec) * Math.Cos(nutobl.epsilon)) - (Math.Cos(celest.dec) * Math.Sin(nutobl.epsilon) * Math.Sin(celest.ra)));

            return coords;
        }

        //Calculate the components  of nutation in longitude and obliquity
        public static NutObl calcNutObl(Double T) {
            NutObl nutobl = new NutObl();

            double L = deg2rad(280.4665 + (36000.7698 * T));
            double Ld = deg2rad(218.3165 + (481267.8813 * T));
            double omega = deg2rad(125.0445479 - (1934.1362891 * T) + (0.0020754 * T * T) + ((T * T * T) / 467441) - ((T * T * T * T) / 60616000));

            nutobl.epsilon0 = deg2rad(23.4392911111 - (0.013004166667 * T) - ((1.63889 / 10000000) * T * T) + (T * T * T * (5.03611 / 10000000)));

            nutobl.deltaepsilon = deg2rad(((9.2 / 3600) * (Math.Cos(omega))) + ((0.57 / 3600) * (Math.Cos(2 * L))) + ((0.1 / 3600) * (Math.Cos(2 * Ld))) - ((0.09 / 3600) * (Math.Cos(2 * omega))));
            nutobl.deltapsi = deg2rad(((-17.2 / 3600) * (Math.Sin(omega))) - ((1.32 / 3600) * (Math.Sin(2 * L))) - ((2.23 / 3600) * (Math.Sin(2 * Ld))) + ((0.21 / 3600) * (Math.Sin(2 * omega))));

            nutobl.epsilon = nutobl.epsilon0 + nutobl.deltaepsilon;

            return nutobl;

        }

        //Calculate the corrections in RA/Dec due to nutation
        public static EquatorialCoords Nutation(EquatorialCoords celest, NutObl nutobl) {

            EquatorialCoords delta = new EquatorialCoords();

            delta.ra = nutobl.deltapsi * (Math.Cos(nutobl.epsilon) + Math.Sin(nutobl.epsilon) * Math.Sin(celest.ra) * Math.Tan(celest.dec)) - nutobl.deltaepsilon * Math.Cos(celest.ra) * Math.Tan(celest.dec);
            delta.dec = nutobl.deltapsi * Math.Sin(nutobl.epsilon) * Math.Cos(celest.ra) + nutobl.deltaepsilon * Math.Sin(celest.ra);

            return delta;
        }

        //Calculate the corrections in RA/Dec due to aberation
        public static EquatorialCoords Aberation(EquatorialCoords celest, EclipticalCoords sun, NutObl nutobl, double T) {

            EquatorialCoords delta = new EquatorialCoords();

            double k = deg2rad(20.49552 / 3600);
            double ee = deg2rad(0.016708634 - 4.2037E-05 * T - 1.267E-07 * T * T);
            double pir = deg2rad(102.93735 + 1.71946 * T + 0.00046 * T * T);

            delta.ra = (-k / Math.Cos(celest.ra)) * (Math.Cos(celest.ra) * Math.Cos(sun.l) * Math.Cos(nutobl.epsilon) + Math.Sin(celest.ra) * Math.Sin(sun.l));
            delta.ra += (ee * k / Math.Cos(celest.dec)) * (Math.Cos(celest.ra) * Math.Cos(pir) * Math.Cos(nutobl.epsilon) + Math.Sin(celest.ra) * Math.Sin(pir));
            delta.dec = ee * k * (Math.Cos(pir) * Math.Cos(nutobl.epsilon) * (Math.Tan(nutobl.epsilon) * Math.Cos(celest.dec) - Math.Sin(celest.ra) * Math.Sin(celest.dec)) + Math.Cos(celest.ra) * Math.Sin(celest.dec) * Math.Sin(pir));
            delta.dec -= k * (Math.Cos(sun.l) * Math.Cos(nutobl.epsilon) * (Math.Tan(nutobl.epsilon) * Math.Cos(celest.dec) - Math.Sin(celest.ra) * Math.Sin(celest.dec)) + Math.Cos(celest.ra) * Math.Sin(celest.dec) * Math.Sin(sun.l));

            return delta;
        }

    }


}

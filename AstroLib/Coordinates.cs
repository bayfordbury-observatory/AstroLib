using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace AstroLib {
    
    //Everything must be in radians

    public class EquatorialCoords {

        public double ra { get; set; }
        public double dec { get; set; }

        public EquatorialCoords() {
            ra = 0;
            dec = 0;
        }

        public EquatorialCoords(double new_ra, double new_dec) {
            ra = new_ra;
            dec = new_dec;
        }

        public HorizontalCoords ToHorizontal(Location loc, double T) {

            HorizontalCoords coordsout = new HorizontalCoords();

            double theta = Utils.deg2rad(280.46061837 + (360.98564736629 * T * 36525) + (0.000387933 * T * T) - ((T * T * T) / 38710000));

            double HA = theta - loc.lon - this.ra;

            coordsout.az = Math.Atan2(Math.Sin(HA), Math.Cos(HA) * Math.Sin(loc.lat) - Math.Tan(this.dec) * Math.Cos(loc.lat)) + Math.PI;

            coordsout.alt = Math.Asin(Math.Sin(loc.lat) * Math.Sin(this.dec) + Math.Cos(loc.lat) * Math.Cos(this.dec) * Math.Cos(HA));

            return coordsout;

        }

        public CoordsString ToDMS() {

            CoordsString StringOut = new CoordsString();

            double RAH = Utils.quad(this.ra.ToDeg()) / 15;
            double DecD = (this.dec.ToDeg());

            StringOut.ra = Utils.hh(RAH) + "h " + Utils.dm(RAH) + "m " + Utils.ds(RAH) + "s";
            StringOut.dec = Utils.dd(DecD) + (Char)176 + " " + Utils.dm(DecD) + (Char)39 + " " + Utils.ds(DecD) + (Char)34;

            return StringOut;
        }

        public EquatorialCoords Precess(double JD0, double JD) {

            EquatorialCoords coordsout = new EquatorialCoords();

            double T = (JD0 - 2451545.0) / 36525;
            double t = (JD - JD0) / 36525;

            double eta = (((2306.2181 + 1.39656 * T - 0.000139 * T * T) * t + (0.30188 - 0.000344 * T) * t * t + 0.017998 * t * t * t) / 3600) * Math.PI / 180.0;
            double zeta = (((2306.2181 + 1.39656 * T - 0.000139 * T * T) * t + (1.09468 + 0.000066 * T) * t * t + 0.018203 * t * t * t) / 3600) * Math.PI / 180.0;
            double theta = (((2004.3109 - 0.8533 * T - 0.000217 * T * T) * t - (0.42665 + 0.000217 * T) * t * t - 0.041833 * t * t * t) / 3600) * Math.PI / 180.0;

            double A = Math.Cos(this.dec) * Math.Sin(this.ra + eta);
            double B = Math.Cos(theta) * Math.Cos(this.dec) * Math.Cos(this.ra + eta) - Math.Sin(theta) * Math.Sin(this.dec);
            double C = Math.Sin(theta) * Math.Cos(this.dec) * Math.Cos(this.ra + eta) + Math.Cos(theta) * Math.Sin(this.dec);

            coordsout.ra = Math.Atan2(A, B) + zeta;

            if (this.dec > 1.5) { //~86 degres
                coordsout.dec = Math.Acos(Math.Sqrt(A * A + B * B));
            } else {
                coordsout.dec = Math.Asin(C);
            }

            ra = coordsout.ra;
            dec = coordsout.dec;

            return coordsout;
        }

        //Calculate the corrections in RA/Dec due to nutation
        public EquatorialCoords Nutation(NutObl nutobl) {

            EquatorialCoords delta = new EquatorialCoords();

            delta.ra = nutobl.deltapsi * (Math.Cos(nutobl.epsilon) + Math.Sin(nutobl.epsilon) * Math.Sin(this.ra) * Math.Tan(this.dec)) - nutobl.deltaepsilon * Math.Cos(this.ra) * Math.Tan(this.dec);
            delta.dec = nutobl.deltapsi * Math.Sin(nutobl.epsilon) * Math.Cos(this.ra) + nutobl.deltaepsilon * Math.Sin(this.ra);

            ra = ra + delta.ra;
            dec = dec + delta.dec;

            return delta;
        }

        //Calculate the corrections in RA/Dec due to aberation
        public EquatorialCoords Aberation(EclipticalCoords sun, NutObl nutobl, double T) {

            EquatorialCoords delta = new EquatorialCoords();

            double k = (20.49552 / 3600) * Math.PI / 180.0;
            double ee = (0.016708634 - 4.2037E-05 * T - 1.267E-07 * T * T) * Math.PI / 180.0;
            double pir = (102.93735 + 1.71946 * T + 0.00046 * T * T) * Math.PI / 180.0;

            delta.ra = (-k / Math.Cos(this.ra)) * (Math.Cos(this.ra) * Math.Cos(sun.l) * Math.Cos(nutobl.epsilon) + Math.Sin(this.ra) * Math.Sin(sun.l));
            delta.ra += (ee * k / Math.Cos(this.dec)) * (Math.Cos(this.ra) * Math.Cos(pir) * Math.Cos(nutobl.epsilon) + Math.Sin(this.ra) * Math.Sin(pir));
            delta.dec = ee * k * (Math.Cos(pir) * Math.Cos(nutobl.epsilon) * (Math.Tan(nutobl.epsilon) * Math.Cos(this.dec) - Math.Sin(this.ra) * Math.Sin(this.dec)) + Math.Cos(this.ra) * Math.Sin(this.dec) * Math.Sin(pir));
            delta.dec -= k * (Math.Cos(sun.l) * Math.Cos(nutobl.epsilon) * (Math.Tan(nutobl.epsilon) * Math.Cos(this.dec) - Math.Sin(this.ra) * Math.Sin(this.dec)) + Math.Cos(this.ra) * Math.Sin(this.dec) * Math.Sin(sun.l));

            ra = ra + delta.ra;
            dec = dec + delta.dec;

            return delta;
        }

        public EclipticalCoords ToEcliptical(NutObl nutobl) {

            EclipticalCoords coords = new EclipticalCoords();

            coords.l = Math.Atan2((Math.Sin(this.ra) * Math.Cos(nutobl.epsilon)) + (Math.Tan(this.dec) * Math.Sin(nutobl.epsilon)), Math.Cos(this.ra));
            coords.b = Math.Asin((Math.Sin(this.dec) * Math.Cos(nutobl.epsilon)) - (Math.Cos(this.dec) * Math.Sin(nutobl.epsilon) * Math.Sin(this.ra)));

            return coords;
        }

        public EquatorialCoords J2000ToB1950() {

            double JD0 = 2451545.0; //j2000
            double JD = 2433282.4235; //b1950

            return Precess(JD0, JD);
        }

        public EquatorialCoords B1950ToJ2000() {

            double JD0 = 2433282.4235; //b1950
            double JD = 2451545.0; //j2000

            return Precess(JD0, JD);
        }

        public EquatorialCoords B1950ToNow(double JD) {

            double JD0 = 2433282.4235; //b1950

            return Precess(JD0, JD);
        }

        public EquatorialCoords NowToJ2000(double JD0) {

            double JD = 2451545.0; //j2000

            return Precess(JD0, JD);
        }

        public EquatorialCoords J2000ToNow(double JD) {

            double JD0 = 2451545.0;

            return Precess(JD0, JD);
        }

        public GalacticCoords B1950ToGalactic() {

            GalacticCoords coordsout = new GalacticCoords();

            double x = Math.Atan2(Math.Sin(Utils.deg2rad(192.25) - this.ra), Math.Cos(Utils.deg2rad(192.25) - this.ra) * Math.Sin(Utils.deg2rad(27.4)) - Math.Tan(this.dec) * Math.Cos(Utils.deg2rad(27.4)));
            coordsout.l = Utils.deg2rad(Utils.quad(303 - Utils.rad2deg(x)));

            double sinB = Math.Sin(this.dec) * Math.Sin(Utils.deg2rad(27.4)) + Math.Cos(this.dec) * Math.Cos(Utils.deg2rad(27.4)) * Math.Cos(Utils.deg2rad(192.25) - this.ra);
            coordsout.b = Math.Asin(sinB);

            return coordsout;
        }

     }

    public class CoordsString {

        public string ra { get; set; }
        public string dec { get; set; }
        public string l { get; set; }
        public string b { get; set; }
        public string r{ get; set; }


    }

    public class HorizontalCoords {

        public double alt { get; set; }
        public double az { get; set; }

        public HorizontalCoords() {
            alt = 0;
            az = 0;
        }

        public HorizontalCoords(double new_alt, double new_az) {
            alt = new_alt;
            az = new_az;
        }

        public EquatorialCoords ToEquatorial(Location loc, double T) {

            EquatorialCoords coordsout = new EquatorialCoords();

            double HA = Math.Atan2(Math.Sin(this.az + Math.PI), Math.Cos(this.az + Math.PI) * Math.Sin(loc.lat) + Math.Tan(this.alt) * Math.Cos(loc.lat));

            double theta = (280.46061837 + (360.98564736629 * T * 36525) + (0.000387933 * T * T) - ((T * T * T) / 38710000)) * Math.PI / 180.0;

            coordsout.ra = theta - loc.lon - HA;

            coordsout.dec = Math.Asin(Math.Sin(loc.lat) * Math.Sin(this.alt) - Math.Cos(loc.lat) * Math.Cos(this.alt) * Math.Cos(this.az + Math.PI));

            return coordsout;
        }    

    }

    public class GalacticCoords {

        public double l { get; set; }
        public double b { get; set; }

        public GalacticCoords() {
            l = 0;
            b = 0;
        }

        public GalacticCoords(double new_l, double new_b) {
            l = new_l;
            b = new_b;
        }

        public EquatorialCoords ToB1950() {

            EquatorialCoords coordsout = new EquatorialCoords();

            double y = Math.Atan2(Math.Sin(this.l - Utils.deg2rad(123)), Math.Cos(this.l - Utils.deg2rad(123)) * Math.Sin(Utils.deg2rad(27.4)) - Math.Tan(this.b) * Math.Cos(Utils.deg2rad(27.4)));
            coordsout.ra = Utils.deg2rad(Utils.quad(12.25 + Utils.rad2deg(y)));

            double sind = Math.Sin(this.b) * Math.Sin(Utils.deg2rad(27.4)) + Math.Cos(this.b) * Math.Cos(Utils.deg2rad(27.4)) * Math.Cos(this.l - Utils.deg2rad(123));
            coordsout.dec = Math.Asin(sind);

            return coordsout;
        }

        public CoordsString ToDMS() {

            CoordsString StringOut = new CoordsString();

            double l = Utils.quad(this.l.ToDeg()) ;
            double b = (this.b.ToDeg());

            StringOut.l= Utils.dd(l) + (Char)176 + " " + Utils.dm(l) + (Char)39 + " " + Utils.ds(l) + (Char)34;
            StringOut.b = Utils.dd(b) + (Char)176 + " " + Utils.dm(b) + (Char)39 + " " + Utils.ds(b) + (Char)34;

            return StringOut;
        }


    }

    public class Location {

        public double lat { get; set; }
        public double lon { get; set; } //Longitude is POSITIVE to the WEST
        public double alt { get; set; } //metres      

        public Location() {
            lat = 0;
            lon = 0;
            alt = 0;
        }

        public Location(double new_lat, double new_lon, double new_alt) {
            lat = new_lat;
            lon = new_lon;
            alt = new_alt;
        }

        public Location(double new_lat, double new_lon) {
            lat = new_lat;
            lon = new_lon;
            alt = 0;
        }

    }

    public class EclipticalCoords {

        public double l { get; set; }
        public double b { get; set; }
        public double r { get; set; }

        public EclipticalCoords() {
            l = 0;
            b = 0;
            r = 0;
        }

        public EclipticalCoords(double new_l, double new_b, double new_r) {
            l = new_l;
            b = new_b;
            r = new_r;
        }

        public EclipticalCoords(double new_l, double new_b) {
            l = new_l;
            b = new_b;
            r = 0;
        }

        public EclipticalCoords HelToGeo(EclipticalCoords Earth) {

            EclipticalCoords geo = new EclipticalCoords();

            double x = this.r * Math.Cos(this.b) * Math.Cos(this.l) - Earth.r * Math.Cos(Earth.b) * Math.Cos(Earth.l);
            double y = this.r * Math.Cos(this.b) * Math.Sin(this.l) - Earth.r * Math.Cos(Earth.b) * Math.Sin(Earth.l);
            double z = this.r * Math.Sin(this.b) - Earth.r * Math.Sin(Earth.b);

            geo.l = Math.Atan2(y, x);
            geo.b = Math.Atan2(z, Math.Sqrt(x * x + y * y));
            geo.r = Math.Sqrt(x * x + y * y + z * z);

            l = geo.l;
            b = geo.b;
            r = geo.r;

            return geo;
        }


        public EquatorialCoords EclipticalToEquatorial(NutObl nutobl) {

            EquatorialCoords coords = new EquatorialCoords();

            coords.ra = Math.Atan2((Math.Sin(this.l) * Math.Cos(nutobl.epsilon)) - (Math.Tan(this.b) * Math.Sin(nutobl.epsilon)), Math.Cos(this.l));
            coords.dec = Math.Asin((Math.Sin(this.b) * Math.Cos(nutobl.epsilon)) + (Math.Cos(this.b) * Math.Sin(nutobl.epsilon) * Math.Sin(this.l)));

            return coords;
        }

        public CoordsString ToDMS() {

            CoordsString StringOut = new CoordsString();

            double ld = Utils.quad(this.l.ToDeg()) ;
            double bd= (this.b.ToDeg());

            StringOut.l = Utils.dd(ld) + (Char)176 + " " + Utils.dm(ld) + (Char)39 + " " + Utils.ds(ld) + (Char)34;
            StringOut.b = Utils.dd(bd) + (Char)176 + " " + Utils.dm(bd) + (Char)39 + " " + Utils.ds(bd) + (Char)34;

            return StringOut;
        }

        public EquatorialCoords HelToEquatorial(EclipticalCoords Earth, NutObl nutobl){

            EclipticalCoords geo = Utils.HelToGeo(Earth, this);

            geo.r += nutobl.deltapsi;

            return geo.EclipticalToEquatorial(nutobl);

        }

    }

    public class NutObl {

        public double epsilon0 { get; set; }
        public double deltaepsilon { get; set; }
        public double epsilon { get; set; }
        public double deltapsi { get; set; }

        public NutObl() {

            epsilon0 = 0;
            deltaepsilon = 0;
            epsilon = 0;
            deltapsi = 0;

        }

        public NutObl(Double T) {

            double L = Utils.deg2rad(280.4665 + (36000.7698 * T));
            double Ld = Utils.deg2rad(218.3165 + (481267.8813 * T));

            double omega = Utils.deg2rad(125.04452 - (1934.136261 * T) + (0.0020708 * T * T) + ((T * T * T) / 450000));

            epsilon0 = (23.4392911111 - (0.01300416666 * T) - (0.00000016388888 * T * T) + (T * T * T * 0.00000050361111)) * Math.PI / 180.0;

            deltapsi = (((-17.2 * Math.Sin(omega)) - (1.32 * Math.Sin(2 * L)) - (2.23 * Math.Sin(2 * Ld)) + (0.21 * Math.Sin(2 * omega))) / 3600) * Math.PI / 180.0;
            deltaepsilon = (((9.2 * Math.Cos(omega)) + (0.57 * Math.Cos(2 * L)) + (0.1 * Math.Cos(2 * Ld)) - (0.09 * Math.Cos(2 * omega))) / 3600) * Math.PI / 180.0;

            epsilon = epsilon0 + deltaepsilon;

        }

    }

    public class MoonSun {

        public double Ld { get; set; } //Moon's mean longitude
        public double D { get; set; }  //Mean elongation of the Moon
        public double M { get; set; }  //Sun's mean anomoly
        public double Md { get; set; } //Moon's mean anomoly
        public double F { get; set; }  //Moon's argument of latitude
        public double omega { get; set; } //Mean longitude of the ascending node of the lunar orbit
        public double E { get; set; } //Eccentricity of the earth

        public MoonSun() {

            Ld = 0;
            D = 0;
            M = 0;
            Md = 0;
            F = 0;
            E = 0; ;
            omega = 0;

        }

        public MoonSun(double T) {

            Ld = (218.3164477 + (481267.88123421 * T) - (0.0015786 * T * T) + ((T * T * T) / 538841) - ((T * T * T * T) / 65194000)) * Math.PI / 180.0;
            D = (297.8501921 + (445267.1114034 * T) - (0.0018819 * T * T) + ((T * T * T) / 545868) - ((T * T * T * T) / 113065000)) * Math.PI / 180.0;
            M = (357.5291092 + (35999.0502909 * T) - (0.0001536 * T * T) + ((T * T * T) / 24490000)) * Math.PI / 180.0;
            Md = (134.9633964 + (477198.8675055 * T) + (0.0087414 * T * T) + ((T * T * T) / 69699) - ((T * T * T * T) / 14712000)) * Math.PI / 180.0;
            F = (93.272095 + (483202.0175233 * T) - (0.0036539 * T * T) + ((T * T * T) / 3526000) - ((T * T * T * T) / 863310000)) * Math.PI / 180.0;
            E = 1 - (0.002516 * T) - (7.4E-06 * T * T);
            omega = (125.0445479 - (1934.1362891 * T) + (0.0020754 * T * T) + ((T * T * T) / 467441) - ((T * T * T * T) / 60616000)) * Math.PI / 180.0;

        }
    }

 


}

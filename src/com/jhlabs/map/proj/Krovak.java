package com.jhlabs.map.proj;

import static java.lang.Math.atan2;
import static java.lang.Math.sin;
import static java.lang.Math.asin;
import static java.lang.Math.cos;
import static java.lang.Math.tan;
import static java.lang.Math.atan;
import static java.lang.Math.sqrt;
import static java.lang.Math.pow;
import static java.lang.Math.abs;
import static java.lang.Math.toRadians;
import static java.lang.Math.toDegrees;

import com.jhlabs.geom.Point2D;
import com.jhlabs.geom.Point3D;

public class Krovak extends Projection {

	static double SEC_TO_RAD = 4.84813681109535993589914102357e-6;

	public boolean tczech = false;

	// towgs84 parameters
	double dx = 485.021;
	double dy = 169.465;
	double dz = 483.839;

	double rx = 7.786342 * SEC_TO_RAD;
	double ry = 4.397554 * SEC_TO_RAD;
	double rz = 4.102655 * SEC_TO_RAD;

	double m = 0.0;

	private Point3D wgs_xyz = new Point3D();
	private Point3D xyz_wgs = new Point3D();
	
	/**
	 * 
	 * @param lp
	 *            x - longitude (lambda) y - latitude (phi)
	 * @param xy
	 */
	public final Point2D transform(Point2D lp, Point2D xy) {
		//long t1 = System.currentTimeMillis();

		double s45, s90, e2, e, alfa, u0, g, k, n0, ro0, ad, a, s0, n;
		double gfi, u, fi0, v, s, d, eps, ro;
		double phi, lam;
		
		phi = toRadians(lp.y);
		lam = toRadians(lp.x);

		BLH_xyz(new Point3D(lam, phi, 0), wgs_xyz);
		transform(wgs_xyz);
		xyz_BLH(wgs_xyz, xyz_wgs);
		phi = xyz_wgs.y;
		lam = xyz_wgs.x;

		s45 = 0.785398163397448;
		s90 = 2 * s45;
		fi0 = projectionLatitude;//phi0;

		a = 6377397.155;
		e2 = 0.006674372230614;
		e = sqrt(e2);

		alfa = sqrt(1. + (e2 * pow(cos(fi0), 4)) / (1. - e2)); // B

		u0 = asin(sin(fi0) / alfa); // gama0
		g = pow((1. + e * sin(fi0)) / (1. - e * sin(fi0)), alfa * e / 2.); // in

		k = tan(u0 / 2. + s45) / pow(tan(fi0 / 2. + s45), alfa) * g; // t0

		n0 = a * sqrt(1. - e2) / (1. - e2 * pow(sin(fi0), 2)); // A
		s0 = 1.37008346281555; // phi1
		n = sin(s0); // n
		ro0 = scaleFactor * n0 / tan(s0);
		ad = s90 - 1.04216856380474; //(true) azimuth of initial line passing through the projection center
		gfi = pow(((1. + e * sin(phi)) / (1. - e * sin(phi))), (alfa * e / 2.));

		u = 2. * (atan(k * pow(tan(phi / 2. + s45), alfa) / gfi) - s45); // U
		v = (projectionLongitude - lam) * alfa;

		s = asin(cos(ad) * sin(u) + sin(ad) * cos(u) * cos(v));
		d = asin(cos(u) * sin(v) / cos(s));
		eps = n * d;
		ro = ro0 * pow(tan(s0 / 2. + s45), n) / pow(tan(s / 2. + s45), n);

		xy.y = ro * cos(eps);
		xy.x = ro * sin(eps);
		//xy.y = falseNorthing + ro * cos(eps);
		//xy.x = falseEasting + ro * sin(eps);
		
		if (!tczech) {
			xy.y *= -1.0;
			xy.x *= -1.0;
		}
		//System.out.println("projecting time: "+(System.currentTimeMillis()-t1));
		return xy;
	}

	// vypocet pravouhlych souradnic z geodetickych souradnic
	private void BLH_xyz(Point3D lph, Point3D xyz) {
		double a = 6378137.0;
		double f1 = 298.257223563;

		double ro, e2;
		e2 = 1.0 - pow(1.0 - 1.0 / f1, 2);
		ro = a / sqrt(1.0 - e2 * pow(sin(lph.y), 2));
		xyz.x = (ro + lph.z) * cos(lph.y) * cos(lph.x);
		xyz.y = (ro + lph.z) * cos(lph.y) * sin(lph.x);
		xyz.z = ((1.0 - e2) * ro + lph.z) * sin(lph.y);
	}

	private void transform(Point3D xyz) {
		double s = (-m / 1000000.0) + 1.0;
		double x2 = -dx + s * ( xyz.x    + rz*xyz.y - ry*xyz.z);
		double y2 = -dy + s * (-rz*xyz.x + xyz.y    + rx*xyz.z);
		double z2 = -dz + s * ( ry*xyz.x - rx*xyz.y + xyz.z);
		xyz.set(x2, y2, z2);
	}
	
	// vypocet geodetickych souradnic z pravouhlych souradnic
	private void xyz_BLH(Point3D xyz, Point3D lph) {
		double a = 6377397.15508; // parametry Besselova elipsoidu
		double f1 = 299.152812853;
		double ab, e2, th, st, ct, p, t;
		ab = f1 / (f1 - 1.0);
		p = sqrt(pow(xyz.x, 2) + pow(xyz.y, 2));
		e2 = 1.0 - pow(1.0 - 1.0 / f1, 2);
		th = atan(xyz.z * ab / p);
		st = sin(th);
		ct = cos(th);
		t = (xyz.z + e2 * ab * a * pow(st, 3)) / (p - e2 * a * pow(ct, 3));

		lph.y = atan(t);
		lph.z = sqrt(1 + t * t) * (p - a / sqrt(1 + (1 - e2) * t * t));
		lph.x = 2 * atan(xyz.y / (p + xyz.x));
	}

	public final Point2D inverseTransform(Point2D xy, Point2D lp) {
		double s45, s90, fi0, e2, e, alfa, u0, g, k, n0, ro0, ad, a, s0, n;
		double u, deltav, s, d, eps, ro, fi1, xy0;

		s45 = 0.785398163397448;
		s90 = 2 * s45;
		fi0 = projectionLatitude;

		a = 6377397.15508;
		e2 = 0.006674372230614;
		e = sqrt(e2);

		alfa = sqrt(1. + (e2 * pow(cos(fi0), 4)) / (1. - e2));
		u0 = asin(sin(fi0) / alfa);
		g = pow((1. + e * sin(fi0)) / (1. - e * sin(fi0)), alfa * e / 2.);

		k = tan(u0 / 2. + s45) / pow(tan(fi0 / 2. + s45), alfa) * g;

		n0 = a * sqrt(1. - e2) / (1. - e2 * pow(sin(fi0), 2));
		s0 = 1.37008346281555;
		n = sin(s0);
		ro0 = scaleFactor * n0 / tan(s0);
		ad = s90 - 1.04216856380474;

		xy0 = xy.x;
		xy.x = xy.y;
		xy.y = xy0;

		if (!tczech) {
			xy.x *= -1.0;
			xy.y *= -1.0;
		}

		ro = sqrt(xy.x * xy.x + xy.y * xy.y);
		eps = atan2(xy.y, xy.x);
		d = eps / sin(s0);
		s = 2. * (atan(pow(ro0 / ro, 1. / n) * tan(s0 / 2. + s45)) - s45);

		u = asin(cos(ad) * sin(s) - sin(ad) * cos(s) * cos(d));
		deltav = asin(cos(s) * sin(d) / cos(u));

		lp.x = projectionLongitude - deltav / alfa;

		fi1 = u;
		int ok = 0;
		do {
			lp.y = 2.0 * (atan(pow(k, -1.0 / alfa) * pow(tan(u / 2.0 + s45), 1.0 / alfa)
					* pow((1.0 + e * sin(fi1)) / (1.0 - e * sin(fi1)), e / 2.0)) - s45);

			if (abs(fi1 - lp.y) < 0.000000000000001) {
				ok = 1;
			}
			fi1 = lp.y;
		} while (ok == 0);

		double f1, x1, y1, z1, x2, y2, z2;
		a = 6377397.15508;
		f1 = 299.152812853;
		e2 = 1 - pow(1. - 1. / f1, 2);
		ro = a / sqrt(1.0 - e2 * pow(sin(lp.y), 2));
		double H = 0;
		x1 = (ro+H) * cos(lp.y) * cos(lp.x);  
		y1 = (ro+H) * cos(lp.y) * sin(lp.x);  
		z1 = ((1-e2) * ro+H) * sin(lp.y);
		
		double scale = (m / 1000000.0) + 1.0;
		
		x2 = dx + scale * ( x1    - rz*y1 + ry*z1);
		y2 = dy + scale * ( rz*x1 + y1    - rx*z1);
		z2 = dz + scale * (-ry*x1 + rx*y1 + z1);
		
		a = 6378137.0;
		f1 = 298.257223563;
		double a_b = f1/(f1-1.0);
		double p = sqrt(x2*x2+y2*y2);
		e2 = 1.0 - pow(1.0 - 1.0 / f1, 2);
		double theta = atan(z2 * a_b / p);
		double st = sin(theta);
		double ct = cos(theta);
		double t = (z2+e2*a_b*a*pow(st, 3))/(p-e2*a*pow(ct, 3));
		lp.y = atan(t);
		lp.x = 2*atan(y2/(p+x2));
		//H = sqrt(1+t*t)*(p-a/sqrt(1+(1-e2)*t*t));
		
		lp.x = toDegrees(lp.x);
		lp.y = toDegrees(lp.y);
		return lp;
	}

	public static void main(String[] args) {
		Krovak krov = new Krovak();
		krov.setScaleFactor(0.9999);
		krov.setProjectionLatitude(0.86393797973719311);
		krov.setProjectionLongitude(0.43342343091192509);
		
		Point2D wgsPos = new Point2D(21.23886386, 49.00096926);
		Point2D out = new Point2D();
		krov.transform(wgsPos, out);
		//krov.transform(new Point2D(16.849771965029522, 50.209011566079397), out);
		System.out.println(wgsPos);
		System.out.println(out);
		
		Point2D xy = new Point2D();
		krov.inverseTransform(out, xy);
		System.out.println(xy);
	}
}

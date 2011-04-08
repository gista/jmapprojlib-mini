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

import com.jhlabs.geom.Point2D;
import com.jhlabs.geom.Point3D;
import com.jhlabs.map.Ellipsoid;

public class KrovakSk extends Projection {

	public static double SEC_TO_RAD = 4.84813681109535993589914102357e-6;
	public static final double s45 = 0.785398163397448;
	public static final double s90 = 1.5707963267948959;

	// towgs84 parameters
	double dx = 485.021;
	double dy = 169.465;
	double dz = 483.839;

	double rx = 7.786342 * SEC_TO_RAD;
	double ry = 4.397554 * SEC_TO_RAD;
	double rz = 4.102655 * SEC_TO_RAD;

	double m = 0.0;

	private double n0;   // A
	private double alfa; // B
	private double u0;   // gama0
	private double k;    // t0
	private double s0 = 1.37008346281555; // phi1 (Latitude of pseudo standard parallel 78° 30'00" N)
	private double n = sin(s0); // n
	private double ad = s90 - 1.04216856380474; //(true) azimuth of initial line passing through the projection center

	private Point3D lph = new Point3D();
	private Point3D xyz = new Point3D();

	public KrovakSk() {
		setEllipsoid(Ellipsoid.BESSEL);
	}

	@Override
	public void initialize() {
		super.initialize();
		// projectionLatitude is fi0
		n0 = sqrt(1.0 - es) / (1.0 - es * pow(sin(projectionLatitude), 2)); // A (a=1)
		alfa = sqrt(1.0 + (es * pow(cos(projectionLatitude), 4)) / (1.0 - es)); // B
		u0 = asin(sin(projectionLatitude) / alfa); // gama0
		
		double g = pow((1.0 + e * sin(projectionLatitude)) / (1.0 - e * sin(projectionLatitude)), alfa * e / 2.0);
		k = tan(u0 / 2.0 + s45) / pow(tan(projectionLatitude / 2.0 + s45), alfa) * g; // t0
	}

	/**
	 * 
	 * @param lp
	 *            x - longitude (lambda)
	 *            y - latitude (phi)
	 * @param xy
	 */
	public final Point2D project(double lam, double phi, Point2D xy) {
		double ro0, gfi, u, v, s, d, eps, ro;

		// apply towgs84 datum shift
		lph.set(lam + projectionLongitude, phi, 0);
		Ellipsoid.WGS_1984.toGeocentric(lph, xyz);
		shiftDatum(xyz);
		ellipsoid.toGeodetic(xyz, lph);
		lam = lph.x - projectionLongitude;
		phi = lph.y;
		
		ro0 = scaleFactor * n0 / tan(s0);
		gfi = pow(((1.0 + e * sin(phi)) / (1.0 - e * sin(phi))), (alfa * e / 2.0));
		u = 2.0 * (atan(k * pow(tan(phi / 2.0 + s45), alfa) / gfi) - s45);
		v = - lam * alfa;
		s = asin(cos(ad) * sin(u) + sin(ad) * cos(u) * cos(v));
		d = asin(cos(u) * sin(v) / cos(s));
		eps = n * d;
		ro = ro0 * pow(tan(s0 / 2.0 + s45), n) / pow(tan(s / 2.0 + s45), n);

		xy.y = ro * cos(eps);
		xy.x = ro * sin(eps);
		
		xy.y *= -1.0;
		xy.x *= -1.0;
		return xy;
	}

	public final Point2D projectInverse(double x, double y, Point2D lp) {
		double ro0, u, deltav, s, d, eps, ro, fi1;

		x = -x;
		y = -y;
		
		ro0 = scaleFactor * n0 / tan(s0);
		ro = sqrt(x * x + y * y);
		eps = atan2(x, y);
		d = eps / sin(s0);
		s = 2.0 * (atan(pow(ro0 / ro, 1.0 / n) * tan(s0 / 2.0 + s45)) - s45);
		u = asin(cos(ad) * sin(s) - sin(ad) * cos(s) * cos(d));
		deltav = asin(cos(s) * sin(d) / cos(u));

		lp.x = projectionLongitude - deltav / alfa;

		fi1 = u;
		boolean ok = false;
		do {
			lp.y = 2.0 * (atan(pow(k, -1.0 / alfa) * pow(tan(u / 2.0 + s45), 1.0 / alfa)
					* pow((1.0 + e * sin(fi1)) / (1.0 - e * sin(fi1)), e / 2.0)) - s45);

			if (abs(fi1 - lp.y) < 0.000000000000001) {
				ok = true;
			}
			fi1 = lp.y;
		} while (!ok);

		// apply towgs84 datum shift
		lph.set(lp.x, lp.y, 0);
		ellipsoid.toGeocentric(lph, xyz);
		inverseShiftDatum(xyz);
		Ellipsoid.WGS_1984.toGeodetic(xyz, lph);
		lp.x = lph.x - projectionLongitude;
		lp.y = lph.y;
		return lp;
	}

	private void shiftDatum(Point3D xyz) {
		double s = (-m / 1000000.0) + 1.0;
		double x2 = -dx + s * ( xyz.x    + rz*xyz.y - ry*xyz.z);
		double y2 = -dy + s * (-rz*xyz.x + xyz.y    + rx*xyz.z);
		double z2 = -dz + s * ( ry*xyz.x - rx*xyz.y + xyz.z);
		xyz.set(x2, y2, z2);
	}

	private void inverseShiftDatum(Point3D xyz) {
		double scale = (m / 1000000.0) + 1.0;
		
		double x2 = dx + scale * ( xyz.x    - rz*xyz.y + ry*xyz.z);
		double y2 = dy + scale * ( rz*xyz.x + xyz.y    - rx*xyz.z);
		double z2 = dz + scale * (-ry*xyz.x + rx*xyz.y + xyz.z);
		xyz.set(x2, y2, z2);
	}

	public boolean hasInverse() {
		return true;
	}

	public boolean isConformal() {
		return true;
	}

	public static void main(String[] args) {
		KrovakSk krov = new KrovakSk();
		krov.setScaleFactor(0.9999);
		krov.setProjectionLatitude(0.86393797973719311);
		krov.setProjectionLongitude(0.43342343091192509);
		krov.initialize();

		Point2D wgsPos = new Point2D(21.23886386, 49.00096926);
		Point2D out = new Point2D();
		krov.transform(wgsPos, out);
		//krov.transform(new Point2D(16.849771965029522, 50.209011566079397), out);
		System.out.println(wgsPos);
		System.out.println("EXPECTED: [-262731.724126, -1208353.150846]");
		System.out.println(out);
		
		Point2D xy = new Point2D();
		krov.inverseTransform(out, xy);
		System.out.println(xy);
		
		// test2
		System.out.println();
		wgsPos = new Point2D(21.127741, 49.0156973);
		krov.transform(wgsPos, out);
		System.out.println(wgsPos);
		System.out.println("EXPECTED: [-270773.835347, -1206330.287303]");
		System.out.println(out);
		krov.inverseTransform(out, xy);
		System.out.println(xy);
	}
}

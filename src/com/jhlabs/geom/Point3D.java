package com.jhlabs.geom;

public final class Point3D {

	public double x;
	public double y;
	public double z;
	
	public Point3D() {
	}
	
	public Point3D(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public void set(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	public void set(Point3D p) {
		this.x = p.x;
		this.y = p.y;
		this.z = p.z;
	}

	@Override
	public final String toString() {
		return String.format("[%f, %f, %f]", x, y, z);
	}
}

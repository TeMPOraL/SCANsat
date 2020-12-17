#region license
/*
 * [Scientific Committee on Advanced Navigation]
 * 			S.C.A.N. Satellite
 * 
 * SCANmapProjection - SCANsat projection types
 *
 * Copyright (c)2013 damny;
 * Copyright (c)2014 technogeeky <technogeeky@gmail.com>;
 * Copyright (c)2014 (Your Name Here) <your email here>; see LICENSE.txt for licensing details.
*/
#endregion

using System;
using UnityEngine;

namespace SCANsat.SCAN_Map
{
	/// <summary>Available map projections.</summary>
	public enum MapProjection
	{
		/// <summary>
		/// A <a href="https://en.wikipedia.org/wiki/Equirectangular_projection"><i>plate carrée</i> equirectangular projection</a>.
		/// <para>Interprets longitudes and latitudes as <c>x</c> and <c>y</c> coordinates.</para>
		/// <para>Uses fixed coorinates of origin (0° N, 0° W).</para>
		/// </summary>
		Rectangular = 0,
		/// <summary>
		/// A general-purpose, low-distortion cylindrical projection with flat poles, projecting pole lines to be half the length of the equator.
		/// <para>See <a href="https://en.wikipedia.org/wiki/Kavrayskiy_VII_projection">general summary</a>.</para>
		/// <para>Uses fixed coorinates of origin (0° N, 0° W).</para>
		/// </summary>
		KavrayskiyVII = 1,

		/// <summary>
		/// A compound projection made by applying <see cref="Orthographic"/> to southern and northern hemispheres separately.
		/// <para>Applies offsets to projected coordinates to form two images - one of each hemisphere - on the same projection plane.</para>
		/// <para>Uses fixed coorinates of origin -90° N, 0° W and 90° N, 0° W.</para>
		/// </summary>
		Polar = 2,

		/// <summary>
		/// A projection of the globe as seen from infinity (approximating deep space), above a given geographic position of origin (longitude/latitude).
		/// <para>Projects at most a hemisphere's worth of data. Has low distortion near the centre, making it suitable for zoomed map.
		/// See <a href="http://wiki.gis.com/wiki/index.php/Orthographic_projection_(cartography)">general summary</a>.</para>
		/// </summary>
		Orthographic = 3,
	}

	public class SCANmapProjection
	{
		public static string[] projectionNames = getProjectionNames();

		private static string[] getProjectionNames()
		{
			MapProjection[] v = (MapProjection[])Enum.GetValues(typeof(MapProjection));
			string[] r = new string[v.Length];
			for (int i = 0; i < v.Length; ++i)
				r[i] = v[i].ToString();
			return r;
		}

		/// <summary>A fixed multiplier applied to coordinates projected via <see cref="MapProjection.Polar"/>.</summary>
		/// <remarks>This seems to serve only to make the polar projection map bigger.</remarks>
		private const double polarProjectionScale = 1.3;
		/// <summary>
		/// A fixed multiplier applied to coordinates projected via <see cref="MapProjection.Orthographic"/>.
		/// </summary>
		private const double orthographicProjectionScale = 1.5;

		/// <summary>
		/// Project geographical coordinates to cartesian coordinates using a given <see cref="MapProjection"/>.
		/// <para>Some projections do not support recentering by <c>centerLon</c> and <c>centerLan</c>. In those cases,
		/// these parameters will be ignored.</para>
		/// </summary>
		/// <param name="projection">Projection to use.</param>
		/// <param name="lon">Longitude in degrees.</param>
		/// <param name="lat">Latitude in degrees.</param>
		/// <param name="easting">Cartesian <c>x</c> coordinate. Generally expected to be in range of [-180, 180]; can be -200 on projection error.</param>
		/// <param name="northing">Cartesian <c>y</c> coordinate. Generally expected to be in range of [-90, 90]; can be -200 on projection error.</param>
		/// <param name="centerLon">Projection longitude of origin, in degrees.</param>
		/// <param name="centerLat">Projection latitude of origin, in degrees.</param>
		/// <remarks>Calculations seem to hardcode some assumptions about the projection plane - particularly, that it'll be 360 grid cells wide, and 180 grid cells high.</remarks>
		public static void project(MapProjection projection, double lon, double lat, out double easting, out double northing, double centerLon = 0.0, double centerLat = 0.0)
		{
			//FIXME understand this. Seems to correct for... something, probably having to do with raster-space coords.
			//NOTE: SCANutil.cs has a bunch of functions doing similar things in there, albeit with 360/180 vs. 3600/1800.
			lon = (lon + 3600 + 180) % 360 - 180; // = (lon + 360*(10 + 0.5)) % 360 
			lat = (lat + 1800 + 90) % 180 - 90;

			double lonR = Mathf.Deg2Rad * lon;
			double latR = Mathf.Deg2Rad * lat;

			switch (projection)
			{
				case MapProjection.KavrayskiyVII:
					easting = Mathf.Rad2Deg * (3.0f * lonR / 2.0f / Math.PI) * Math.Sqrt(Math.PI * Math.PI / 3.0f - latR * latR);
					northing = lat;
					return;
				case MapProjection.Polar:
					if (latR < 0)
					{
						easting = Mathf.Rad2Deg * (polarProjectionScale * Math.Cos(latR) * Math.Sin(lonR) - Math.PI / 2);
						northing = Mathf.Rad2Deg * (polarProjectionScale * Math.Cos(latR) * Math.Cos(lonR));
					}
					else
					{
						easting = Mathf.Rad2Deg * (polarProjectionScale * Math.Cos(latR) * Math.Sin(lonR) + Math.PI / 2);
						northing = Mathf.Rad2Deg * (-polarProjectionScale * Math.Cos(latR) * Math.Cos(lonR));
					}
					return;
				case MapProjection.Orthographic:
					double centerLonR = Mathf.Deg2Rad * centerLon;
					double centerLatR = Mathf.Deg2Rad * centerLat;

					if (Math.Sin(centerLatR) * Math.Sin(latR) + Math.Cos(centerLatR) * Math.Cos(latR) * Math.Cos(lonR - centerLonR) < 0)
					{
						//XXX some kind of weird error handler?
						easting = -200;
						northing = -200;
						return;
					}

					easting = Mathf.Rad2Deg * (orthographicProjectionScale * Math.Cos(latR) * Math.Sin(lonR - centerLonR));
					northing = Mathf.Rad2Deg * (orthographicProjectionScale * (Math.Cos(centerLatR) * Math.Sin(latR) - Math.Sin(centerLatR) * Math.Cos(latR) * Math.Cos(lonR - centerLonR)));
					return;
				case MapProjection.Rectangular:
				default:
					easting = lon;
					northing = lat;
					return;
			}
		}

		/// <summary>
		/// Project cartesian coordinates to geographical coordinates (longitude, latitude) using a given <see cref="MapProjection"/>.
		/// <para>Some projections do not support recentering by <c>centerLon</c> and <c>centerLan</c>. In those cases,
		/// these parameters will be ignored.</para>
		/// </summary>
		/// <param name="projection">Projection to use.</param>
		/// <param name="easting">Cartesian <c>x</c> coordinate.</param>
		/// <param name="northing">Cartesian <c>y</c> coordinate.</param>
		/// <param name="lon">Longitude in degrees. Usually between -180° and 180° (inclusive), but can be = 300 in case of projection error.</param>
		/// <param name="lat">Latitude in degrees. Usually between -90° and 90° (inclusive), but can be = 300 in case of projection error.</param>
		/// <param name="centerLon">Projection longitude of origin, in degrees.</param>
		/// <param name="centerLat">Projection latitude of origin, in degrees.</param>
		/// <remarks>Calculations seem to hardcode some assumptions about the projection plane - particularly, that it'll be 360 grid cells wide, and 180 grid cells high.</remarks>
		public static void unproject(MapProjection projection, double easting, double northing, out double lon, out double lat, double centerLon = 0.0, double centerLat = 0.0)
        {
			//FIXME Again, magic correction block. Need to understand the meaning of this.
			//This has definitely something to do with raster-space.
            if (northing > 90)
            {
                northing = 180 - northing;
                easting += 180;
            }
            else if (northing < -90)
            {
                northing = -180 - northing;
                easting += 180;
            }
            easting = (easting + 3600 + 180) % 360 - 180;
            northing = (northing + 1800 + 90) % 180 - 90;

			double eastingR = Mathf.Deg2Rad * easting;
			double northingR= Mathf.Deg2Rad * northing;

			switch (projection)
            {
                case MapProjection.KavrayskiyVII:
                    lon = Mathf.Rad2Deg * (eastingR / Math.Sqrt(Mathf.PI * Math.PI / 3.0f - northingR * northingR) * 2.0f * Math.PI / 3.0f);
					lat = northing;
					return;
                case MapProjection.Polar:
                    double lat0R = Math.PI / 2;
                    if (eastingR < 0)
                    {
                        eastingR += Math.PI / 2;
                        lat0R = -Math.PI / 2;
                    }
                    else
                    {
                        eastingR -= Math.PI / 2;
                    }
                    eastingR /= polarProjectionScale;
                    northingR /= polarProjectionScale;

                    double p = Math.Sqrt(eastingR * eastingR + northingR * northingR);
                    double c = Math.Asin(p);

                    eastingR = Math.Atan2((eastingR * Math.Sin(c)), (p * Math.Cos(lat0R) * Math.Cos(c) - northingR * Math.Sin(lat0R) * Math.Sin(c)));
                    lon = (Mathf.Rad2Deg * eastingR + 180) % 360 - 180;
					if (lon <= -180)
					{
						lon = -180;
					}
					lat = Mathf.Rad2Deg * Math.Asin(Math.Cos(c) * Math.Sin(lat0R) + (northingR * Math.Sin(c) * Math.Cos(lat0R)) / (p));
					return;
                case MapProjection.Orthographic:
                    double centerLonR = Mathf.Deg2Rad * centerLon;
                    double centerLatR = Mathf.Deg2Rad * centerLat;

                    double p2 = Math.Sqrt(eastingR * eastingR + northingR * northingR);
                    double c2 = Math.Asin(p2 / orthographicProjectionScale);

					if (Math.Cos(c2) < 0)
					{
						//XXX some weird error handler.
						lon = 300;
						lat = 300;
						return;
					}

                    eastingR = centerLonR + Math.Atan2(eastingR * Math.Sin(c2), p2 * Math.Cos(c2) * Math.Cos(centerLatR) - northingR * Math.Sin(c2) * Math.Sin(centerLatR));

                    lon = (Mathf.Rad2Deg * eastingR + 180) % 360 - 180;
					if (lon <= -180)
					{
						lon += 360;
					}
					lat = Mathf.Rad2Deg * Math.Asin(Math.Cos(c2) * Math.Sin(centerLatR) + (northingR * Math.Sin(c2) * Math.Cos(centerLatR)) / p2);
					return;
				case MapProjection.Rectangular:
                default:
					lon = easting;
					lat = northing;
					return;
            }
		}
    }
}

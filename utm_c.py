
"""
UTM - Universal Transverse Mercator python Class
Transformasi utm --> latlon and latlon --> utm
written by : Eko Rahmadi

"""
from math import floor, pi, degrees, radians, sqrt, sin, cos, tan, asin, acos, atan, atan2

class utm_:
    def __init__(self, ellipsoid):
        if ellipsoid.lower()=='wgs84':
            self.a =6378137.0
            self.f = 1 / 298.257223563
            self.e2 = 1 - (1 - self.f) * (1 - self.f)
        elif ellipsoid.lower()=='grs80':
            self.a =6378137.0
            self.f = 1 / 298.257222101
            self.e2 = 1 - (1 - self.f) * (1 - self.f)
    def ll2zone(self, lat, lon):

        if 56 <= lat < 64 and 3 <= lon < 12:
            return 32

        if 72 <= lat <= 84 and lon >= 0:
            if lon < 9:
                return 31
            elif lon < 21:
                return 33
            elif lon < 33:
                return 35
            elif lon < 42:
                return 37

        return int((lon + 180) / 6) + 1


    def central_meredian(self, zone):
        return (zone - 1) * 6 - 180 + 3


    def ll2utm(self, lat, lon):


        K0=0.9996
        E2 =self.e2 * self.e2
        E3 = E2 * self.e2
        
        E_P2 =self.e2 / (1.0 - self.e2)
        
        M1 = (1 - self.e2 / 4 - 3 * E2 / 64 - 5 * E3 / 256)
        M2 = (3 * self.e2 / 8 + 3 * E2 / 32 + 45 * E3 / 1024)
        M3 = (15 * E2 / 256 + 45 * E3 / 1024)
        M4 = (35 * E3 / 3072)
        
        lat_rad = radians(lat)
        lat_sin = sin(lat_rad)
        lat_cos = cos(lat_rad)

        lat_tan = lat_sin / lat_cos
        lat_tan2 = lat_tan * lat_tan
        lat_tan4 = lat_tan2 * lat_tan2

    
        zone =self.ll2zone(lat, lon)
    
        lon_rad = radians(lon)
        central_lon = self.central_meredian(zone)
        central_lon_rad = radians(central_lon)

        n = self.a / sqrt(1 - self.e2 * lat_sin**2)
        c = E_P2 * lat_cos**2

        a = lat_cos * (lon_rad - central_lon_rad)
        a2 = a * a
        a3 = a2 * a
        a4 = a3 * a
        a5 = a4 * a
        a6 = a5 * a

        m = self.a * (M1 * lat_rad -
             M2 * sin(2 * lat_rad) +
             M3 * sin(4 * lat_rad) -
             M4 * sin(6 * lat_rad))

        easting = K0 * n * (a +
                        a3 / 6 * (1 - lat_tan2 + c) +
                        a5 / 120 * (5 - 18 * lat_tan2 + lat_tan4 + 72 * c - 58 * E_P2)) + 500000

        northing = K0 * (m + n * lat_tan * (a2 / 2 +
                                        a4 / 24 * (5 - lat_tan2 + 9 * c + 4 * c**2) +
                                        a6 / 720 * (61 - 58 * lat_tan2 + lat_tan4 + 600 * c - 330 * E_P2)))
        hemi='N'
        if lat<0:
            northing += 10000000
            hemi='S'

        return easting, northing, zone, hemi

    def utm2ll(self, easting, northing, zone, hemi):
        K0=0.9996
        E2 =self.e2 * self.e2
        E3 = E2 * self.e2
        
        E_P2 =self.e2 / (1.0 - self.e2)
        
        _E = (1 - sqrt(1 - self.e2)) / (1 + sqrt(1 - self.e2))
        _E2 = _E * _E
        _E3 = _E2 * _E
        _E4 = _E3 * _E
        _E5 = _E4 * _E
        M1 = (1 - self.e2 / 4 - 3 * E2 / 64 - 5 * E3 / 256)
        P2 = (3. / 2 * _E - 27. / 32 * _E3 + 269. / 512 * _E5)
        P3 = (21. / 16 * _E2 - 55. / 32 * _E4)
        P4 = (151. / 96 * _E3 - 417. / 128 * _E5)
        P5 = (1097. / 512 * _E4)

        x = easting - 500000
        y = northing
        if hemi.upper()=='S':
            y -= 10000000

        m = y / K0
        mu = m / (self.a * M1)

        p_rad = (mu +
             P2 * sin(2 * mu) +
             P3 * sin(4 * mu) +
             P4 * sin(6 * mu) +
             P5 * sin(8 * mu))

        p_sin = sin(p_rad)
        p_sin2 = p_sin * p_sin

        p_cos = cos(p_rad)

        p_tan = p_sin / p_cos
        p_tan2 = p_tan * p_tan
        p_tan4 = p_tan2 * p_tan2

        ep_sin = 1 - self.e2 * p_sin2
        ep_sin_sqrt = sqrt(1 - self.e2 * p_sin2)

        n = self.a / ep_sin_sqrt
        r = (1 - self.e2) / ep_sin

        c = _E * p_cos**2
        c2 = c * c

        d = x / (n * K0)
        d2 = d * d
        d3 = d2 * d
        d4 = d3 * d
        d5 = d4 * d
        d6 = d5 * d

        lat = (p_rad - (p_tan / r) *
                (d2 / 2 -
                 d4 / 24 * (5 + 3 * p_tan2 + 10 * c - 4 * c2 - 9 * E_P2)) +
                 d6 / 720 * (61 + 90 * p_tan2 + 298 * c + 45 * p_tan4 - 252 * E_P2 - 3 * c2))

        lon = (d -
                 d3 / 6 * (1 + 2 * p_tan2 + c) +
                 d5 / 120 * (5 - 2 * c + 28 * p_tan2 - 3 * c2 + 8 * E_P2 + 24 * p_tan4)) / p_cos

        return (degrees(lat), degrees(lon)+self.central_meredian(zone))

def main():
    '''
    #Testing
    gc=utm_('wgs84')
    print(gc.ll2utm(-6.4910555327, 106.8489137899))
    print(gc.utm2ll(704462.4938, 9282139.5482, 48,'S'))
    print("y = ", gc.ll2utm(-6.4910555327, 106.8489137899)[1])
    '''
if __name__ == "__main__":
    main()

    

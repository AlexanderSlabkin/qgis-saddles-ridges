import math
import time


# rlayer = QgsRasterLayer("D:/GEOINF/Elevation GMTED2010/GMTED2010N10E030_300/10n030e_20101117_gmted_mea300.tif", "Elevation") #africa

# rlayer = QgsRasterLayer("D:/GEOINF/Elevation GMTED2010/Austria/GMTED2010N30E000_075/30n000e_20101117_gmted_mea075.tif", "Elevation") #austria
# vlayer = QgsVectorLayer("D:/GEOINF/lines.shp", "lines",'ogr')
rlayer = QgsRasterLayer("D:/GEOINF/Elevation GMTED2010/Austria/GMTED2010N30E000_300/30n000e_20101117_gmted_mea300.tif",
                        "Elevation")  # austria
vlayer = QgsVectorLayer("D:/GEOINF/lines.shp", "lines", 'ogr')
if not vlayer.isValid():
    print
    "vlayer failed to load!"
if not rlayer.isValid():
    print
    "rlayer failed to load!"

if QgsMapLayerRegistry.instance().count() is 0:
    crs = QgsCoordinateReferenceSystem(4326, QgsCoordinateReferenceSystem.EpsgCrsId)
    rlayer.setCrs(crs)
    QgsMapLayerRegistry.instance().addMapLayer(rlayer)
    QgsMapLayerRegistry.instance().addMapLayer(vlayer)

xmin = rlayer.extent().xMinimum()
xmax = rlayer.extent().xMaximum()
ymin = rlayer.extent().yMinimum()
ymax = rlayer.extent().yMaximum()

xstep = rlayer.rasterUnitsPerPixelX()
ystep = rlayer.rasterUnitsPerPixelY()
print xstep


def elev(point):
    ident = rlayer.dataProvider().identify(point, QgsRaster.IdentifyFormatValue)
    if ident.isValid():
        return ident.results()[1]
    else:
        print
        'elev: error'


a = [QgsVector(0, ystep), QgsVector(xstep, ystep), QgsVector(xstep, 0), QgsVector(xstep, -ystep), QgsVector(0, -ystep),
     QgsVector(-xstep, -ystep), QgsVector(-xstep, 0), QgsVector(-xstep, ystep)]


def is_loc_max(point):
    max = elev(point)
    for i in range(8):
        if (elev(point + a[i]) >= max):
            return False
    return True


def is_loc_min(point):
    min = elev(point)
    for i in range(8):
        if (elev(point + a[i]) <= min):
            return False
    return True


def find_loc_maxs(point1, vector1, vector2):
    if not (abs((abs(vector1.angle(vector2)) - math.pi / 2)) < 0.001 or abs(
            (abs(vector1.angle(vector2)) - math.pi * 3 / 2)) < 0.001):
        print "find_loc_maxs error: non perp"
        return 1
    point = QgsPoint(point1)
    pointx = point.x()
    pointy = point.y()
    alpha = vector1.angle()
    betha = vector2.angle()
    if abs(alpha - math.pi) < 0.1 or abs(alpha - math.pi / 2) < 0.1 or abs(alpha - math.pi * 3 / 2) < 0.1 or abs(
            alpha) < 0.1:
        step1 = vector1.normalized() * xstep
        step2 = vector2.normalized() * ystep
    elif abs(alpha - math.pi / 4) < 0.1 or abs(alpha - math.pi * 3 / 4) < 0.1 or abs(
            alpha - math.pi * 5 / 4) < 0.1 or abs(alpha - math.pi * 7 / 4) < 0.1:
        step1 = vector1.normalized() * xstep / math.cos(math.pi / 4)
        step2 = vector1.normalized() * ystep / math.cos(math.pi / 4)
    elif alpha > math.pi / 4 and alpha < math.pi * 3 / 4 or alpha > 5 * math.pi / 4 and alpha < math.pi * 7 / 4:
        alpha = alpha if alpha < math.pi else alpha - math.pi
        betha = betha if betha < math.pi else betha - math.pi
        step1 = vector1.normalized() * (xstep / math.cos(math.pi / 2 - alpha))
        step2 = vector2.normalized() * (ystep / abs(math.cos(betha)))
    else:
        alpha = alpha if alpha < math.pi else alpha - math.pi
        betha = betha if betha < math.pi else betha - math.pi
        step1 = vector1.normalized() * (xstep / abs(math.cos(alpha)))
        step2 = vector2.normalized() * (ystep / math.cos(math.pi / 2 - betha))
    Nstep1 = int(vector1.length() / step1.length())
    Nstep2 = int(vector2.length() / step2.length())
    maxs = []
    for i in range(Nstep2):
        for j in range(Nstep1):
            if is_loc_max(point):
                maxs.append(QgsPoint(point))
            point += step1
        point = QgsPoint(pointx, pointy) + step2 * (i + 1)
    return maxs


def find_loc_mins(xmin, xmax, ymin, ymax):
    x = xmin
    y = ymin
    mins = []
    while (ymax - y) > ystep / 2:
        while (xmax - x) > xstep / 2:
            point = QgsPoint(x, y)
            if is_loc_min(point):
                mins.append(point)
            x += xstep
        x = xmin
        y += ystep
    return mins


def find_loc_maxs_and_mins(point1, vector1, vector2):
    if not (abs((abs(vector1.angle(vector2)) - math.pi / 2)) < 0.001 or abs(
            (abs(vector1.angle(vector2)) - math.pi * 3 / 2)) < 0.001):
        print
        "find_loc_maxs error: non perp"
        return 1
    point = QgsPoint(point1)
    pointx = point.x()
    pointy = point.y()
    alpha = vector1.angle()
    betha = vector2.angle()
    if abs(alpha - math.pi) < 0.1 or abs(alpha - math.pi / 2) < 0.1 or abs(alpha - math.pi * 3 / 2) < 0.1 or abs(
            alpha) < 0.1:
        step1 = vector1.normalized() * xstep
        step2 = vector2.normalized() * ystep
    elif abs(alpha - math.pi / 4) < 0.1 or abs(alpha - math.pi * 3 / 4) < 0.1 or abs(
            alpha - math.pi * 5 / 4) < 0.1 or abs(alpha - math.pi * 7 / 4) < 0.1:
        step1 = vector1.normalized() * xstep / math.cos(math.pi / 4)
        step2 = vector1.normalized() * ystep / math.cos(math.pi / 4)
    elif alpha > math.pi / 4 and alpha < math.pi * 3 / 4 or alpha > 5 * math.pi / 4 and alpha < math.pi * 7 / 4:
        alpha = alpha if alpha < math.pi else alpha - math.pi
        betha = betha if betha < math.pi else betha - math.pi
        step1 = vector1.normalized() * (xstep / math.cos(math.pi / 2 - alpha))
        step2 = vector2.normalized() * (ystep / abs(math.cos(betha)))
    else:
        alpha = alpha if alpha < math.pi else alpha - math.pi
        betha = betha if betha < math.pi else betha - math.pi
        step1 = vector1.normalized() * (xstep / abs(math.cos(alpha)))
        step2 = vector2.normalized() * (ystep / math.cos(math.pi / 2 - betha))
    Nstep1 = int(vector1.length() / step1.length())
    Nstep2 = int(vector2.length() / step2.length())
    maxs = []
    mins = []
    for i in range(Nstep2):
        for j in range(Nstep1):
            if is_loc_max(point):
                maxs.append(QgsPoint(point))
            elif is_loc_min(point):
                mins.append(QgsPoint(point))
            point += step1
        point = QgsPoint(pointx, pointy) + step2 * (i + 1)
    return mins, maxs


left_bottom = QgsPoint(xmin, ymin)


def sortByLen(inputVector):
    return inputVector.length()


class Saddle:
    def __init__(self, min1, min2, max1, max2):
        if min1.distance(left_bottom) <= min2.distance(left_bottom):
            self.min1 = min1
            self.min2 = min2
        else:
            self.min1 = min2
            self.min2 = min1
        if max1.distance(left_bottom) <= max2.distance(left_bottom):
            self.max1 = max1
            self.max2 = max2
        else:
            self.max1 = max2
            self.max2 = max1

    def __str__(self):
        return "min1 = " + str(self.min1) + ", " + str(elev(self.min1)) + \
               "   min2 = " + str(self.min2) + ", " + str(elev(self.min2)) + \
               "   max1 = " + str(self.max1) + ", " + str(elev(self.max1)) + \
               "   max2 = " + str(self.max2) + ", " + str(elev(self.max2))

    def __eq__(self, other):
        if self.min1.compare(other.min1, 0.0001) & \
                self.min2.compare(other.min2, 0.0001) & \
                self.max1.compare(other.max1, 0.0001) & \
                self.max2.compare(other.max2, 0.0001):
            return True
        else:
            return False

    def __sub__(self, other):
        return min((self.min1 - other.min1), (self.min2 - other.min2), (self.max1 - other.max1),
                   (self.max2 - other.max2), key=sortByLen)

    def draw(self):
        caps = vlayer.dataProvider().capabilities()
        if caps & QgsVectorDataProvider.AddFeatures:
            feat1 = QgsFeature()
            feat2 = QgsFeature()
            feat3 = QgsFeature()
            feat4 = QgsFeature()
            feat1.setGeometry(QgsGeometry.fromPolyline([self.min1, self.min2]))
            feat2.setGeometry(QgsGeometry.fromPolyline([self.max2, self.min2]))
            feat3.setGeometry(QgsGeometry.fromPolyline([self.min1, self.max1]))

            (res, outFeats) = vlayer.dataProvider().addFeatures([feat1, feat2, feat3, feat4])
        return outFeats


def sort_by_distance_to_point(list, point):
    def sortByDistance(inputPoint):
        return point.distance(inputPoint)

    list.sort(key=sortByDistance)
    return list


def find_Saddles(xmin, xmax, ymin, ymax):
    mins1 = find_loc_mins(xmin, xmax, ymin, ymax)
    print len(mins1)
    saddles = []
    mins = sort_by_distance_to_point(mins1, left_bottom)
    for point in mins:
        b = sort_by_distance_to_point(mins, point)
        k = 1
        distance = point.distance(b[k])
        while (distance < 0.1):
            if (abs(elev(point) - elev(b[k])) < 1000):
                vec1 = (point - b[k]) * 0.5
                vec2 = vec1.perpVector()
                q = b[k] + vec1
                z1 = q + vec1 / 2 + vec2 / 2
                z2 = q - vec1 / 2 - vec2 / 2
                (center_mins, center_maxs) = find_loc_maxs_and_mins(point - vec2 / 2, -vec1 * 2, vec2)
                right = find_loc_maxs(z1, -vec1, vec2)
                left = find_loc_maxs(z2, vec1, -vec2)
                if len(right) > 0 and len(left) > 0 and (
                        len(center_mins) is 0 or elev(min(center_mins, key=elev)) > elev(b[k]) and elev(
                        min(center_mins, key=elev)) > elev(point)):
                    max_right = max(right, key=elev)
                    max_left = max(left, key=elev)
                    if len(center_maxs) > 0:
                        max_center = max(center_maxs, key=elev)
                        if elev(max_right) > elev(max_center) and elev(max_left) > elev(max_center):
                            saddles.append(Saddle(point, b[k], max_left, max_right))


                    else:
                        saddles.append(Saddle(point, b[k], max_left, max_right))
            k += 1
            distance = point.distance(b[k])
        mins.remove(point)

    return saddles


ii = []
for i in range(1000):
    ii.append(i)
res = vlayer.dataProvider().deleteFeatures(ii)
# saddles1 = find_Saddles(6, 6.2, 46, 46.2)
##saddles1 = find_Saddles(36, 37, 11, 13)
saddles1 = find_Saddles(7, 9, 45, 47)
print
len(saddles1)
for i in saddles1:
    # print i
    i.draw()
# duplicates = []
# for i in range(len(saddles1)):
#    for j in range(len(saddles1)):
#        if not(i is  j) and abs((saddles1[i] - saddles1[j]).length()) < xstep and not ((j,i) in duplicates):
#            duplicates.append((i,j))
# print duplicates
# for i in duplicates:
#    saddles1[i[0]].draw()
#    saddles1[i[1]].draw()
#
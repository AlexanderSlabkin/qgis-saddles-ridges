import math
import time


# rlayer = QgsRasterLayer("D:/GEOINF/Elevation GMTED2010/Austria/GMTED2010N30E000_300/30n000e_20101117_gmted_mea300.tif", "Elevation") #austria
rlayer = QgsRasterLayer("D:/GEOINF/Elevation GMTED2010/Austria/GMTED2010N30E000_075/30n000e_20101117_gmted_mea075.tif",
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
    for i in a:
        if elev(point + i) > max:
            return False
    return True


def is_loc_min(point):
    min = elev(point)
    for i in a:
        if elev(point + i) <= min:
            return False
    return True


def is_same_value(point):
    value = elev(point)
    for i in a:
        if elev(point + i) == value:
            return True
    return False


step1 = QgsVector(xstep, 0)
step2 = QgsVector(0, ystep)


def find_loc_maxs(startpoint, stoppoint):
    xmin = startpoint.x()
    point = QgsPoint(startpoint)
    maxs = []
    N1 = int((stoppoint.x() - startpoint.x()) / xstep)
    N2 = int((stoppoint.y() - startpoint.y()) / ystep)
    for i in range(N2):
        for j in range(N1):
            if is_loc_max(point):
                maxs.append(QgsPoint(point))
            point += step1
        point.setX(xmin)
        point += step2
    return maxs


def draw_ridge(points):
    for i in range(len(points) - 1):
        feat = QgsFeature()
        feat.setGeometry(QgsGeometry.fromPolyline([points[i], points[i + 1]]))
        (res, outFeats) = vlayer.dataProvider().addFeatures([feat])


def mean(startpoint, stoppoint):
    xmin = startpoint.x()
    point = QgsPoint(startpoint)
    string_mean = 0
    N1 = int((stoppoint.x() - startpoint.x()) / xstep)
    N2 = int((stoppoint.y() - startpoint.y()) / ystep)
    for i in range(N2):
        sum = 0
        for j in range(N1):
            sum += elev(point)
            point += step1
        string_mean += sum / N1
        point.setX(xmin)
        point += step2
    return string_mean / N2


def deviation(mean, startpoint, stoppoint):
    xmin = startpoint.x()
    point = QgsPoint(startpoint)
    string_mean = 0
    N1 = int((stoppoint.x() - startpoint.x()) / xstep)
    N2 = int((stoppoint.y() - startpoint.y()) / ystep)
    for i in range(N2):
        sum = 0
        for j in range(N1):
            sum += abs(elev(point) - mean)
            point += step1
        string_mean += sum / N1
        point.setX(xmin)
        point += step2
    return string_mean / N2


def gradient(point, index, eps):
    max = elev(point + a[index])
    z = index
    indexes = [(index - 1) % 8, index, (index + 1) % 8]
    for i in indexes:  # , index - 2, (index + 2)%8]:
        m = elev(point + a[i])
        if m > max:
            max = m
    if max - elev(point) < eps:
        result = []
        ind = 0
        indexes.append(index)
        for i in indexes:
            if elev(point + a[i]) == max:
                result.append(i)
        if len(result) > 1:
            max = elev(point + a[result[0]] * 2)
            for j in range(1, len(result)):
                if elev(point + a[result[j]] * 2) > max:
                    max = elev(point + a[result[j]] * 2)
                    ind = j
        return result[ind]
    else:
        return False


def find_ridge(point, eps):
    supermax = elev(point)
    max = elev(point + a[0])
    indexes = []
    way_n = 0
    passed_points = []
    for i in range(1, 8):
        m = elev(point + a[i])
        if m > max:
            max = m
    for i in range(8):
        if elev(point + a[i]) == max:
            indexes.append(i)
    for i in indexes:
        k = point + a[i]
        passed_points.append([point, QgsPoint(k)])
        index = i
        number_of_steps = 0
        # eps = 200
        while (True):
            number_of_steps += 1
            next_point = gradient(k, index, eps / 2)
            if next_point and number_of_steps < 100:
                k += a[next_point]
                index = next_point
                passed_points[way_n].append(QgsPoint(k))
                #                if is_loc_max(k):
                #                    break
                if supermax - elev(k) > eps:
                    break
            else:
                break
        way_n += 1

    return passed_points


def is_point_in_list(points_list, point_to_find):
    for point in points_list:
        if point_to_find.compare(point, 0.0001):
            return True
    return False


def divide_to_small_squares(point1, point2):
    points = []
    if point1.distance(point2) <= 0.141421356237:
        return point1, point2
    else:
        x1 = point1.x()
        x2 = point2.x()
        y1 = point1.y()
        y2 = point2.y()
        n = int((x2 - x1) / 0.1)
        m = int((y2 - y1) / 0.1)
        for j in range(m):
            for i in range(n):
                points.append((QgsPoint(x1 + i * 0.1, y1 + j * 0.1), QgsPoint(x1 + (i + 1) * 0.1, y1 + (j + 1) * 0.1)))
        return points


def find_Ridges(point1, point2):
    squares = divide_to_small_squares(point1, point2)
    for j in squares:
        m = mean(*j)
        d = deviation(m, *j)
        if d < 70:
            print "flat surface", j
            return 0

    m = mean(point1, point2)
    d = deviation(m, point1, point2)
    if d < 150:
        print
        "flat surface"
        return 0
    maxs = find_loc_maxs(point1, point2)
    ridges = []
    for point in maxs:
        if elev(point) > m - d:
            ridges.append(find_ridge(point, d))
    return ridges


def distance_between(ridge1, ridge2):
    return max(min(ridge1[0].distance(ridge2[0]), ridge1[0].distance(ridge2[-1])),
               min(ridge1[-1].distance(ridge2[0]), ridge1[-1].distance(ridge2[-1])))


def is_same_points(list1, list2):
    n = 0
    for i in list1:
        if is_point_in_list(list2, i):
            n += 1
    return n


def mean_distance(ridge1, ridge2):
    mean = 0
    if (len(ridge1) > len(ridge2)):
        list1 = list(ridge1)
        list2 = list(ridge2)
    else:
        list1 = list(ridge2)
        list2 = list(ridge1)

    for i in list1:
        min = i.distance(list2[0])
        point = list2[0]
        for j in list2:
            if i.distance(j) < min:
                min = i.distance(j)
                point = j
        mean += min
        list2.remove(point)
        if len(list2) == 0:
            break
        mo = len(ridge2)
        mo = len(ridge1) if len(ridge1) < len(ridge2) else len(ridge2)
        return mean / mo


def find_Saddles(point1, point2, distance):
    ridges = find_Ridges(point1, point2)
    if len(ridges) < 2:
        print
        'No saddles'
        return False
    saddles = []
    for k in ridges:
        for q in ridges:
            for i in k:
                for j in q:
                    if not (j, i) in saddles and distance_between(i, j) < distance and not q is k and len(
                            i) > 10 and len(j) > 10 \
                            and not is_same_points(i, j) and mean_distance(i, j) < 0.001:
                        angle = (i[0] - i[-1]).angle(j[0] - j[-1])
                        angle = abs(angle)
                        if (math.pi * 3 / 4 < angle < math.pi * 5 / 4) or (
                                math.pi / 4 > angle > math.pi * 7 / 4):
                            saddles.append((i, j))
    return saddles


#
# point3 = QgsPoint(9, 45)
# point4 = QgsPoint(10, 46)
point3 = QgsPoint(6, 44)
point4 = QgsPoint(7, 45)
# point3 = QgsPoint(5.5, 45)
# point4 = QgsPoint(6.5, 46)
ii = []
for i in range(100000):
    ii.append(i)
res = vlayer.dataProvider().deleteFeatures(ii)
# saddles = find_Saddles(7, 9, 45, 47)

saddles = find_Saddles(point3, point4, 0.05)
print
len(saddles)
if saddles:
    for i in saddles:
        for j in i:
            draw_ridge(j)

# ridges = find_Ridges(point3, point4)
##print len(saddles), len(ridges)
# for i in ridges:
#    for j in i:
#        if len(j) > 5:
#            #print j
#            draw_ridge(j)


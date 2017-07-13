from math import *
from decimal import Decimal, getcontext

getcontext().prec = 13


class Vector(object):
    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)
        except ValueError:
            raise ValueError('The coordinates must be nonempty')
        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return str(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        new_coordinates = [x + y for x, y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def minus(self, v):
        new_coordinates = [x - y for x, y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def times_scalar(self, c):
        new_coordinates = [Decimal(str(c)) * x for x in self.coordinates]
        return Vector(new_coordinates)

    def magnitude(self):
        coordinates_squared = [x ** 2 for x in self.coordinates]
        return sqrt(sum(coordinates_squared))

    def normalized(self):
        try:
            magnitude = self.magnitude()
            return self.times_scalar(Decimal('1.0') / Decimal(str(magnitude)))
        except ZeroDivisionError:
            raise Exception('Cannot normalize the zero vector')

    def dot(self, v):
        return sum([x * y for x, y in zip(self.coordinates, v.coordinates)])

    def arc(self, v):
        u1 = self.normalized()
        u2 = v.normalized()
        try:
            return acos(u1.dot(u2))
        except ZeroDivisionError:
            raise Exception('cannot compute an angle with the zero vector')

    def degree(self, v):
        return self.arc(v) * 180 / pi

    def project(self, v):
        try:
            normv = v.normalized()
            return normv.times_scalar(self.dot(normv))
        except Exception as exception:
            raise exception

    def orth(self, v):
        try:
            return self.minus(self.project(v))
        except Exception as exception:
            raise exception

    def cross(self, v):
        new_coordinates = []
        for i in range(1, self.dimension + 1):
            first = i % self.dimension
            second = (i + 1) % self.dimension
            new_coordinates.append(
                self.coordinates[first] * v.coordinates[second] - self.coordinates[second] * v.coordinates[first])
        return Vector(new_coordinates)

    def parallelogram(self, v):
        return self.cross(v).magnitude()

    def triangle(self, v):
        return self.parallelogram(v) / Decimal('2.0')

# -*- coding:utf-8 -*-

from decimal import Decimal, getcontext
from math import *
from vector import Vector

getcontext().prec = 12


class Plane(object):
    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'

    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

        if not normal_vector:
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    # 一个变量不为0,其他变量取0,可得到平面上一点作为基点(直线参数化后,可得到基点)
    def set_basepoint(self):
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0'] * self.dimension

            initial_index = Plane.first_nonzero_index(n.coordinates)
            initial_coefficient = n.coordinates[initial_index]

            basepoint_coords[initial_index] = c / initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __str__(self):

        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '

            if abs(coefficient) != 1:
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Plane.first_nonzero_index(n.coordinates)
            terms = [write_coefficient(n.coordinates[i], is_initial_term=(i == initial_index)) + 'x_{}'.format(i + 1)
                     for i in range(self.dimension) if round(n.coordinates[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e

        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += ' = {}'.format(constant)

        return output

    @staticmethod
    def first_nonzero_index(iterable):
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():
                return k
        raise Exception(Plane.NO_NONZERO_ELTS_FOUND_MSG)

    def parallel(self, p):
        norm1 = self.normal_vector.normalized()
        norm2 = p.normal_vector.normalized()
        if norm1.arc(norm2) == 0 or norm1.arc(norm2) == pi:
            return "parallel"
        else:
            return "not parallel"

    def __eq__(self, other):
        norm = self.normal_vector
        base1 = self.basepoint
        base2 = other.basepoint
        try:
            return (self.basepoint == other.basepoint and self.constant_term == other.constant_term) or norm.dot(
                base2.minus(base1)) == 0
        except Exception as e:
            return False


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


'''
plane11 = Plane(Vector(['-0.412', '3.806', '0.728']), '-3.46')
plane12 = Plane(Vector(['1.03', '-9.515', '-1.82']), '8.65')
print plane11.parallel(plane12)
print "equal?", (plane11 == plane12)

plane21 = Plane(Vector(['2.611', '5.528', '0.283']), '4.6')
plane22 = Plane(Vector(['7.715', '8.306', '5.342']), '3.76')
print plane21.parallel(plane22)
print "equal?", (plane21 == plane22)

plane31 = Plane(Vector(['-7.926', '8.625', '-7.212']), '-7.952')
plane32 = Plane(Vector(['-2.642', '2.875', '-2.404']), '-2.443')
print plane31.parallel(plane32)
print "equal?", (plane31 == plane32)
'''

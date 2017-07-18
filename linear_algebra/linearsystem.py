# -*- coding:utf-8 -*-

from decimal import Decimal, getcontext

from vector import Vector
from plane import Plane
from copy import deepcopy

getcontext().prec = 30


class LinearSystem(object):
    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        n = self[row].normal_vector
        k = self[row].constant_term
        self[row] = Plane(
            normal_vector=n.times_scalar(coefficient),
            constant_term=k * coefficient)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        n1 = self[row_to_add].normal_vector
        k1 = self[row_to_add].constant_term
        n2 = self[row_to_be_added_to].normal_vector
        k2 = self[row_to_be_added_to].constant_term
        self.planes[row_to_be_added_to] = Plane(
            normal_vector=n1.times_scalar(coefficient).plus(n2),
            constant_term=str(k1 * coefficient + k2))

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def compute_triangular_form(self):
        system = deepcopy(self)

        indices = self.indices_of_first_nonzero_terms_in_each_row()

        for i, v in enumerate(indices):
            # 完全的行简化阶梯型应该满足i(行数)==v(主元所在列数) 例如第0行的主元应该在第0列
            if i < v:
                # 如果i<v,说明此行位置不是标准位置,在整个方程组中寻找主元列位置(v)等于该行索引(i)的方程,进行交换
                for i2, v2 in enumerate(indices):
                    if v2 == i:
                        system.swap_rows(i, i2)
                        # 更新主元位置数据
                        indices = system.indices_of_first_nonzero_terms_in_each_row()
                        break
            elif v < i:
                # 如果i>v,说明应当消去部分变量,例如第2行的第一个不为0的变量位置位置如果是0,那么应该消除位置0,位置1的变量(即2个变量)
                # reduceIndex从当前不为0的位置开始,消除i(当前行索引)个变量
                reduceIndex = v
                while reduceIndex < i:
                    for i2, v2 in enumerate(indices):
                        if v2 == reduceIndex:
                            # 找到主元位置是当前要消去变量位置的方程,计算并消去
                            try:
                                coefficient = system[i].normal_vector.coordinates[reduceIndex] / \
                                              system[i2].normal_vector.coordinates[reduceIndex]
                                system.add_multiple_times_row_to_row(-1 * coefficient, i2, i)
                                # 更新主元位置数据
                                indices = system.indices_of_first_nonzero_terms_in_each_row()
                                break
                            except Exception:
                                continue
                    reduceIndex += 1

        return system

    # 化简成行简化阶梯型.官方教程中从下到上消元的方法更好,只需要两个循环
    def compute_rref(self):
        tf = self.compute_triangular_form()
        rows = len(tf.planes)
        dim = tf.dimension

        # 确定每行主元位置
        indices = tf.indices_of_first_nonzero_terms_in_each_row()

        for i in range(rows):
            # 确认主元位置
            m = indices[i]

            # 主元位置系数不为1,缩放 (m<0说明没有主元)
            if m >= 0 and tf[i].normal_vector.coordinates[m] != 1:
                tf.multiply_coefficient_and_row(Decimal("1") / tf[i].normal_vector.coordinates[m], i)

            for j in range(m + 1, dim):
                # 要化简的位置不为0
                if tf[i].normal_vector.coordinates[j] != 0:
                    # 寻找主元位置在该列的行,与该行进行消元
                    for k, v in enumerate(indices):
                        if v == j:
                            coefficient = tf[i].normal_vector.coordinates[j] / \
                                          tf[k].normal_vector.coordinates[j]
                            tf.add_multiple_times_row_to_row(-1 * coefficient, k, i)

        return tf

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i + 1, p) for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p1 = Plane(normal_vector=Vector(['0.786', '0.786', '0.588']), constant_term='-0.714')
p2 = Plane(normal_vector=Vector(['-0.138', '-0.138', '0.244']), constant_term='0.319')
s = LinearSystem([p1, p2])
ref = s.compute_rref()
print ref
print ref.indices_of_first_nonzero_terms_in_each_row()
print set(range(ref.dimension)) - set(ref.indices_of_first_nonzero_terms_in_each_row())
'''
print "----------------------------------------------------"
p1 = Plane(normal_vector=Vector(['8.631', '5.112', '-1.816']), constant_term='-5.113')
p2 = Plane(normal_vector=Vector(['4.315', '11.132', '-5.27']), constant_term='-6.775')
p3 = Plane(normal_vector=Vector(['-2.158', '3.01', '-1.727']), constant_term='-0.831')
s = LinearSystem([p1, p2, p3])
print s.compute_rref()
print "----------------------------------------------------"
p1 = Plane(normal_vector=Vector(['5.262', '2.739', '-9.878']), constant_term='-3.441')
p2 = Plane(normal_vector=Vector(['5.111', '6.358', '7.638']), constant_term='-2.152')
p3 = Plane(normal_vector=Vector(['2.016', '-9.924', '-1.367']), constant_term='-9.278')
p4 = Plane(normal_vector=Vector(['2.167', '-13.543', '-18.883']), constant_term='-10.567')
s = LinearSystem([p1, p2, p3, p4])
print s.compute_rref()
'''
